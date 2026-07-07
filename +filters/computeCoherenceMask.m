function [norm_trace_map, norm_trace_map_nofilt, norm_regularity_map, ...
    norm_persist_map, norm_persist_map_raw, norm_coherence_map, hot_pixel_accumulator,...
    aperiodic_mask, local_hot_mask, global_hot_mask, secondary_cleaning, filtered_counts_mask] = ...
    computeCoherenceMask(sorted_x, sorted_y, sorted_t, imgSz, ...
    t_interval, coh_params, ...
    frameIndex, norm_trace_map_prev, std_map, mean_map, counts, hot_pixel_accumulator, ...
    genFigure, global_hot_mask)

    % ----------------------------------------------------------------
    % 0. Parse parameters
    % ----------------------------------------------------------------
    r_s                   = coh_params.r_s;
    s_bnd                 = coh_params.s_bnd;
    hpa_decay             = coh_params.hpa_decay;
    hpa_bnd               = coh_params.hpa_bnd;

    % ----------------------------------------------------------------
    % 1. Rule 1 — IEI regularity (regularity map)
    % ----------------------------------------------------------------
    [~, ~, norm_regularity_map] = filters.findRegularity( ...
        sorted_x, sorted_y, std_map, mean_map, imgSz);

    % handle zeros and inf
    norm_regularity_map(isnan(norm_regularity_map)) = 0;
    norm_regularity_map(norm_regularity_map==0)=inf;
    
    % Generate the aperiodic map
    aperiodic_mask = (norm_regularity_map <= s_bnd); 
        
    % Leaky-bucket for persistent tracking
    hot_pixel_accumulator = (hot_pixel_accumulator + aperiodic_mask).*hpa_decay;
    local_hot_mask   = hot_pixel_accumulator >= hpa_bnd;

    if frameIndex == 1
        global_hot_mask = local_hot_mask;
    else
        global_hot_mask = global_hot_mask | local_hot_mask;
    end
    
    % 4. Strip the defective pixels globally
    hot_pixel_idx = find(local_hot_mask);

    filtered_x = sorted_x;
    filtered_y = sorted_y;
    filtered_t = sorted_t;
    
    if ~isempty(hot_pixel_idx)
        % Strip from the coordinate arrays
        sorted_lin_idx = sub2ind(imgSz, filtered_x, filtered_y);
        remove_mask = ismember(sorted_lin_idx, hot_pixel_idx);
        
        filtered_x(remove_mask) = [];
        filtered_y(remove_mask) = [];
        filtered_t(remove_mask) = [];
        
    end

    % ----------------------------------------------------------------
    % 2. Rule 2 - Spatial density (trace map) - Unfiltered
    % ----------------------------------------------------------------

    if genFigure

        % Temporal binning at pixel-equivalent resolution
        Nt_nofilt = max(round(t_interval / (r_s * t_interval)), 1);
        tb_nofilt = min(floor(Nt_nofilt * (sorted_t - min(sorted_t)) / t_interval) + 1, Nt_nofilt);
    
        % 3D event histogram at full spatial resolution
        V_nofilt = accumarray([sorted_x, sorted_y, tb_nofilt], 1, [imgSz(1) imgSz(2) Nt_nofilt]);

        kr_x_nofilt = round(r_s * imgSz(1));
        kr_y_nofilt = round(r_s * imgSz(2));
        kr_t_nofilt = round(r_s * Nt_nofilt);

        % Create a 2D Spatial Kernel (Flat Disk)
        [kx_nofilt, ky_nofilt] = ndgrid(-kr_x_nofilt:kr_x_nofilt, -kr_y_nofilt:kr_y_nofilt);
        K_xy_dist_nofilt = sqrt((kx_nofilt/imgSz(1)).^2 + (ky_nofilt/imgSz(2)).^2);
        K_xy_nofilt = double(K_xy_dist_nofilt <= r_s);

        % Create a 1D Temporal Kernel (Line)
        % Reshaped to 1x1xL to apply strictly along the 3rd dimension (time)
        K_t_nofilt = ones(1, 1, 2 * kr_t_nofilt + 1);

        % Apply sequentially (O(N^2) + O(N) instead of O(N^3))
        density_xy_nofilt = convn(V_nofilt, K_xy_nofilt, 'same');
        density_nofilt = convn(density_xy_nofilt, K_t_nofilt, 'same');
    
        % Sample per-event, then max per-pixel
        sum_exp_dist_nofilt = density_nofilt(sub2ind(size(V_nofilt), sorted_x, sorted_y, tb_nofilt));
    
        sum_exp_dist_map_nofilt = accumarray([sorted_x, sorted_y], sum_exp_dist_nofilt, imgSz, @max, 0);
    
        % Log-normalize to compress the heavy-tailed distribution
        log_trace_map_nofilt = log1p(sum_exp_dist_map_nofilt');
        norm_trace_map_nofilt = log_trace_map_nofilt' ./ max(log_trace_map_nofilt(:));

    else

        norm_trace_map_nofilt = zeros(size(norm_regularity_map));

    end

    % ----------------------------------------------------------------
    % 3. Rule 2 - Spatial density (trace map) - Filtered
    % ----------------------------------------------------------------

    % Temporal binning at pixel-equivalent resolution
    Nt = max(round(t_interval / (r_s * t_interval)), 1);
    tb = min(floor(Nt * (filtered_t - min(filtered_t)) / t_interval) + 1, Nt);
    
    % 3D event histogram at full spatial resolution
    V = accumarray([filtered_x, filtered_y, tb], 1, [imgSz(1) imgSz(2) Nt]);

    kr_x = round(r_s * imgSz(1));
    kr_y = round(r_s * imgSz(2));
    kr_t = round(r_s * Nt);

    % Create a 2D Spatial Kernel (Flat Disk)
    [kx, ky] = ndgrid(-kr_x:kr_x, -kr_y:kr_y);
    K_xy_dist = sqrt((kx/imgSz(1)).^2 + (ky/imgSz(2)).^2);
    K_xy = double(K_xy_dist <= r_s);

    % Create a 1D Temporal Kernel (Line)
    % Reshaped to 1x1xL to apply strictly along the 3rd dimension (time)
    K_t = ones(1, 1, 2 * kr_t + 1);

    % Apply sequentially 
    density_xy = convn(V, K_xy, 'same');
    density = convn(density_xy, K_t, 'same');
    
    % Sample per-event, then max per-pixel
    sum_exp_dist = density(sub2ind(size(V), filtered_x, filtered_y, tb));

    sum_exp_dist_map = accumarray([filtered_x, filtered_y], sum_exp_dist, imgSz, @max, 0);

    % Log-normalize to compress the heavy-tailed distribution
    log_trace_map = log1p(sum_exp_dist_map');
    norm_trace_map = log_trace_map' ./ max(log_trace_map(:));

    % ----------------------------------------------------------------
    % 4. Rule 3a — Temporal persistence
    % ----------------------------------------------------------------
    if frameIndex == 1
        % No previous frame available — use trace map as proxy
        persist_map = norm_trace_map;
        filtered_counts_mask = ones(size(persist_map));
    else
        persist_map = zeros(size(norm_trace_map));

        % KNN search between current and previous trace maps
        [~, ~, minDists, validIdx] = ...
            filters.findPersistenceVectorized(...
            norm_trace_map, norm_trace_map_prev, imgSz);

        if ~isempty(validIdx)
            persist_map(validIdx) = minDists;
        end

        % Remove the hot pixels from the counts mask as well
        filtered_counts_mask = (counts > 0) & ~local_hot_mask;

        % Calculate the exponential decayed persistance to "invert" the meaning
        persist_map = exp(-persist_map / median(persist_map(persist_map > 0))).*(filtered_counts_mask>0);
    end
    
    % Log-normalize the persistence map
    log_persist_map = log1p(persist_map);
    norm_persist_map_raw = log_persist_map ./ max(log_persist_map(:));

    % ----------------------------------------------------------------
    % 4. Rule 3b — Overly persistent single event removal
    % ----------------------------------------------------------------
    if frameIndex == 1
        norm_persist_map = norm_persist_map_raw;
        secondary_cleaning = zeros(size(norm_persist_map));
    else
        secondary_cleaning = (norm_persist_map_raw==1);
        norm_persist_map = norm_persist_map_raw;
        % Apply secondary cleaning to remove overly persistent events
        norm_persist_map(secondary_cleaning==1) = 0;
    end

    % ----------------------------------------------------------------
    % 5. Combine rule maps
    % ----------------------------------------------------------------

    filtered_coherence_map = norm_trace_map ...
        .* norm_persist_map;

    log_coherence_map = log1p(filtered_coherence_map);
    norm_coherence_map = log_coherence_map ./ max(log_coherence_map(:));


end

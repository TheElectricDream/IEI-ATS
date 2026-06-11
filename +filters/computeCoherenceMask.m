function [norm_trace_map, norm_similarity_map, ...
    norm_persist_map, norm_coherence_map, hot_pixel_accumulator,...
    aperiodic_mask, persistent_hot_mask, global_hot_mask] = ...
    computeCoherenceMask(sorted_x, sorted_y, sorted_t, imgSz, ...
    t_interval, coh_params, ...
    frameIndex, norm_trace_map_prev, std_map, mean_map, counts, hot_pixel_accumulator, ...
    plottingFrame, genFigure, plottingType, global_hot_mask)

    % ----------------------------------------------------------------
    % 0. Parse parameters
    % ----------------------------------------------------------------
    r_s                   = coh_params.r_s;

    % ----------------------------------------------------------------
    % 2. Rule 1 — IEI regularity (regularity map)
    % ----------------------------------------------------------------
    [~, ~, norm_similarity_map] = filters.findSimilarities(...
        sorted_x, sorted_y, std_map, mean_map, imgSz);
    norm_similarity_map(isnan(norm_similarity_map)) = 0;

    % ----------------------------------------------------------------
    % Hot Pixel Removal
    % ----------------------------------------------------------------
    % 1. Define the "Poisson Band" 
    % A CV between 0.8 and 1.25 yields a regularity score between ~0.44 and ~0.55.
    poisson_mask = (norm_similarity_map > 0.0) & (norm_similarity_map < 0.8);
    aperiodic_mask = poisson_mask;
    
    % 2. Isolate hyperactivity
    % We only care about the event counts of pixels that lack temporal structure
    poisson_counts = counts(poisson_mask);
    
    if ~isempty(poisson_counts) && length(poisson_counts) > 10
        
        count_th = median(poisson_counts) + 5 * mad(poisson_counts, 1);
        current_hot_mask = poisson_mask & (counts > count_th);

    end
        
    % 2. Update the Leaky Bucket Memory
    % Add current flags, then decay by 10% (multiply by 0.9). 
    % A consistent hot pixel will accumulate up to ~10.0 over time.
    hot_pixel_accumulator = (hot_pixel_accumulator + current_hot_mask) * 0.90;
    
    % 3. Define the Persistent Hot Mask
    % If a pixel's accumulator is > 2.0, it has been flagged in at least 
    % 3 recent frames. It is considered a confirmed hardware defect.
    persistent_hot_mask = hot_pixel_accumulator > 3.0;

    if frameIndex == 1
        global_hot_mask = persistent_hot_mask;
    else
        global_hot_mask = global_hot_mask | persistent_hot_mask;
    end
    
    % 4. Strip the defective pixels globally
    hot_pixel_idx = find(persistent_hot_mask);

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
    % 1. Rule 2 - Spatial density (trace map)
    % ----------------------------------------------------------------

    if frameIndex == plottingFrame && genFigure

        % Temporal binning at pixel-equivalent resolution
        Nt = max(round(t_interval / (r_s * t_interval)), 1);
        tb = min(floor(Nt * (sorted_t - min(sorted_t)) / t_interval) + 1, Nt);

        % 3D event histogram at full spatial resolution
        V = accumarray([sorted_x, sorted_y, tb], 1, [imgSz(1) imgSz(2) Nt]);

        kr_x = round(r_s * imgSz(1));
        kr_y = round(r_s * imgSz(2));
        kr_t = round(r_s * Nt);

        [kx, ky, kt] = ndgrid(-kr_x:kr_x, -kr_y:kr_y, -kr_t:kr_t);

        % Normalize kernel indices back to the same units as the point cloud
        K_dist = sqrt((kx/imgSz(1)).^2 + (ky/imgSz(2)).^2 + (kt/Nt).^2);
        K = double(K_dist <= r_s); % 1 inside the ball, 0 outside

        density = convn(V, K, 'same');

        % Sample per-event, then max per-pixel
        sum_exp_dist = density(sub2ind(size(V), sorted_x, sorted_y, tb));

        sum_exp_dist_map = accumarray([sorted_x, sorted_y], sum_exp_dist, imgSz, @max, 0);

        % Log-normalize to compress the heavy-tailed distribution
        log_trace_map = log1p(sum_exp_dist_map');
        norm_trace_map = log_trace_map' ./ max(log_trace_map(:));

        journal.showScatterPlotOfRuleMaps3D(norm_trace_map,['Norm-Trace-Map-Sample-Density-w-Hot-Pixels-Nominal-Rot' plottingType '-3D.pdf'], false)
    end

    % Temporal binning at pixel-equivalent resolution
    Nt = max(round(t_interval / (r_s * t_interval)), 1);
    tb = min(floor(Nt * (filtered_t - min(filtered_t)) / t_interval) + 1, Nt);
    
    % 3D event histogram at full spatial resolution
    V = accumarray([filtered_x, filtered_y, tb], 1, [imgSz(1) imgSz(2) Nt]);
      
    kr_x = round(r_s * imgSz(1));
    kr_y = round(r_s * imgSz(2));
    kr_t = round(r_s * Nt);
    
    [kx, ky, kt] = ndgrid(-kr_x:kr_x, -kr_y:kr_y, -kr_t:kr_t);
    
    % Normalize kernel indices back to the same units as the point cloud
    K_dist = sqrt((kx/imgSz(1)).^2 + (ky/imgSz(2)).^2 + (kt/Nt).^2);
    K = double(K_dist <= r_s); % 1 inside the ball, 0 outside
    
    density = convn(V, K, 'same');
    
    % Sample per-event, then max per-pixel
    sum_exp_dist = density(sub2ind(size(V), filtered_x, filtered_y, tb));

    sum_exp_dist_map = accumarray([filtered_x, filtered_y], sum_exp_dist, imgSz, @max, 0);

    % Log-normalize to compress the heavy-tailed distribution
    log_trace_map = log1p(sum_exp_dist_map');
    norm_trace_map = log_trace_map' ./ max(log_trace_map(:));

    % ----------------------------------------------------------------
    % 3. Rule 2 — Temporal persistence
    % ----------------------------------------------------------------
    if frameIndex == 1
        % No previous frame available — use trace map as proxy
        persist_map = norm_trace_map;
    else
        persist_map = zeros(size(norm_trace_map));

        % KNN search between current and previous trace maps
        [~, ~, minDists, validIdx] = ...
            filters.findPersistenceVectorized(...
            norm_trace_map, norm_trace_map_prev, imgSz);

        if ~isempty(validIdx)
            persist_map(validIdx) = minDists;
        end

        % Calculate the exponential decayed persistance to "invert" the meaning
        persist_map = exp(-persist_map / median(persist_map(persist_map > 0))).*(counts>0);
    end
    
    % Log-normalize the persistence map
    log_persist_map = log1p(persist_map);
    norm_persist_map = log_persist_map ./ max(log_persist_map(:));

    % ----------------------------------------------------------------
    % 4. Combine rule maps
    % ----------------------------------------------------------------
    % [th_lo_trace, ~] = stats.findElbowThreshold(norm_trace_map,500);
    % [th_lo_persist, ~] = stats.findElbowThreshold(norm_persist_map,500);
    % [th_lo_sim, ~] = stats.findElbowThreshold(norm_similarity_map,500);
    % 
    % norm_trace_map(norm_trace_map < th_lo_trace) = 0;
    % norm_persist_map(norm_persist_map < th_lo_persist) = 0;
    % norm_similarity_map(norm_similarity_map < th_lo_sim) = 0;

    filtered_coherence_map = norm_trace_map ...
        .* norm_persist_map;

    log_coherence_map = log1p(filtered_coherence_map);
    norm_coherence_map = log_coherence_map ./ max(log_coherence_map(:));


end

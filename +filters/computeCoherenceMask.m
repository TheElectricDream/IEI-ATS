function [norm_trace_map, norm_similarity_map, ...
    norm_persist_map, norm_coherence_map, hot_pixel_accumulator] = ...
    computeCoherenceMask(sorted_x, sorted_y, sorted_t, imgSz, ...
    t_interval, unique_idx, pos, group_ends, coh_params, ...
    frameIndex, norm_trace_map_prev, std_map, mean_map, counts, hot_pixel_accumulator)

% COMPUTECOHERENCEMASK  Three-rule coherence filtering pipeline.
%
%   [NORM_TRACE_MAP, NORM_SIMILARITY_MAP, NORM_PERSIST_MAP,
%   NORM_COHERENCE_MAP] = COMPUTECOHERENCEMASK(...) combines three
%   independent coherence rules into a scalar per-pixel score that
%   separates real edge events from noise before temporal surface
%   accumulation.
%
%   Inputs:
%     sorted_x, sorted_y  - [N x 1] Pixel coordinates (row, col),
%                           sorted by linear pixel index.
%     sorted_t             - [N x 1] Timestamps [s], sorted to match.
%     imgSz                - [1 x 2] Image dimensions [nRows, nCols].
%     t_interval           - Scalar frame interval [s].
%     unique_idx           - [K x 1] Linear indices of active pixels.
%     pos                  - [K x 1] Start of each pixel group in the
%                           sorted vectors.
%     group_ends           - [K x 1] End of each pixel group.
%     coh_params           - Struct with fields:
%       .r_s                   - Spatial search radius (normalized).
%     frameIndex           - Current frame number (1-indexed). Rule 2
%                           is skipped on the first frame.
%     norm_trace_map_prev  - [imgSz] Previous frame's trace map for
%                           persistence assessment (Rule 2).
%     std_map              - [imgSz] Per-pixel IEI standard deviation.
%     mean_map             - [imgSz] Per-pixel IEI mean.
%     counts               - [imgSz] Per-pixel event count map; used
%                           to mask persistence scores to active pixels.
%
%   Outputs:
%     norm_trace_map       - [imgSz] Rule 1: spatial density score.
%     norm_similarity_map  - [imgSz] Rule 3: IEI regularity score.
%     norm_persist_map     - [imgSz] Rule 2: temporal persistence.
%     norm_coherence_map   - [imgSz] Combined score (log-normalized
%                           sum of the three rule maps).
%
%   Algorithm:
%     1. Rule 1 — Spatial density: events are binned into a 3D
%        histogram (x, y, temporal-bin) at full spatial resolution.
%        A distance-weighted convolution kernel of radius r_s
%        computes local density. Per-event density is sampled from
%        the volume and aggregated per-pixel via max. The result
%        is log-normalized to compress the heavy-tailed distribution.
%     2. Rule 3 — IEI regularity: local coefficient of variation
%        (CV = sigma/mu) of the IEI map via normalized convolution.
%        High-CV pixels indicate irregular firing and are penalized.
%     3. Rule 2 — Temporal persistence: KNN search between the
%        current and previous trace maps in normalized (row, col,
%        value) space. Close matches receive high scores via
%        exponential decay of distance, masked to active pixels.
%        Skipped on frame 1 (trace map used as proxy).
%     4. Combine: elementwise sum of the three [0, 1] rule maps,
%        followed by log1p normalization of the combined map.
%
%   Notes:
%     - The three rules are independent and additive. Each produces
%       a [0, 1] normalized map. The combined map is thresholded
%       downstream in main.m (not inside this function).
%     - Coordinates: x = row, y = col, sub2ind(imgSz, x, y).
%
%   See also: filters.findSpatialNeighbours,
%             filters.findSimilarities,
%             filters.findPersistenceVectorized

    % ----------------------------------------------------------------
    % 0. Parse parameters
    % ----------------------------------------------------------------
    r_s                   = coh_params.r_s;

    % ----------------------------------------------------------------
    % 2. Rule 3 — IEI regularity (similarity map)
    % ----------------------------------------------------------------
    [~, ~, norm_similarity_map] = filters.findSimilarities(...
        sorted_x, sorted_y, std_map, mean_map, imgSz);
    norm_similarity_map(isnan(norm_similarity_map)) = 0;

    % ----------------------------------------------------------------
    % Bonus Rules: Hot Pixel Removal
    % ----------------------------------------------------------------
    % 1. Define the "Poisson Band" 
    % A CV between 0.8 and 1.25 yields a regularity score between ~0.44 and ~0.55.
    poisson_mask = (norm_similarity_map > 0.2) & (norm_similarity_map < 0.8);
    
    % 2. Isolate hyperactivity
    % We only care about the event counts of pixels that lack temporal structure
    poisson_counts = counts(poisson_mask);
    
    if ~isempty(poisson_counts) && length(poisson_counts) > 10
        
        % (Fallback alternative if elbow fails: robust statistics)
        count_th = median(poisson_counts) + 5 * mad(poisson_counts, 1);
        current_hot_mask = poisson_mask & (counts > count_th);

    end
        
    % 2. Update the Leaky Bucket Memory
    % Add current flags, then decay by 10% (multiply by 0.9). 
    % A consistent hot pixel will accumulate up to ~10.0 over time.
    hot_pixel_accumulator = (hot_pixel_accumulator + current_hot_mask) * 0.9;
    
    % 3. Define the Persistent Hot Mask
    % If a pixel's accumulator is > 2.0, it has been flagged in at least 
    % 3 recent frames. It is considered a confirmed hardware defect.
    persistent_hot_mask = hot_pixel_accumulator > 2.0;
    
    % Combine current transients and historical defects
    final_hot_mask = current_hot_mask | persistent_hot_mask;
    
    % 4. Strip the defective pixels globally
    hot_pixel_idx = find(final_hot_mask);
    
    if ~isempty(hot_pixel_idx)
        % Strip from the coordinate arrays
        sorted_lin_idx = sub2ind(imgSz, sorted_x, sorted_y);
        remove_mask = ismember(sorted_lin_idx, hot_pixel_idx);
        
        sorted_x(remove_mask) = [];
        sorted_y(remove_mask) = [];
        sorted_t(remove_mask) = [];
        
    end

    % ----------------------------------------------------------------
    % 1. Rule 1 — Spatial density (trace map)
    % ----------------------------------------------------------------

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
    temp = norm_persist_map;

    norm_persist_map(norm_persist_map==1)=0;

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

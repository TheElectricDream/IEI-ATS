function [bf_t_mean, bf_t_std, bf_t_mean_diff, bf_t_std_diff, counts, pixel_history] = ...
    backfillIEIStatisticsReference(sorted_t, unique_idx, pos, group_ends, ...
        imgSz, counts, K_buf_backfill, pixel_history, kernel_size)
% FASTIEISTATISTICS Vectorized, stateful implementation of IEI backfill.
%
% pixel_history must be initialized in main.m BEFORE the loop as:
% pixel_history = NaN(prod(imgSz), K_buf_backfill);

    num_pixels = prod(imgSz);
    
    % Initialize outputs
    bf_t_mean      = zeros(imgSz);
    bf_t_std       = zeros(imgSz);
    bf_t_mean_diff = zeros(imgSz);
    bf_t_std_diff  = zeros(imgSz);

    % ---------------------------------------------------------------------
    % Step 1: Update the Rolling History Matrix
    % ---------------------------------------------------------------------
    % For pixels that have new events, push them into their history row
    for k = 1:numel(unique_idx)
        p_idx = unique_idx(k);
        if counts(p_idx) == 0, continue; end % Skipped by mask
        
        % Extract current frame events for this pixel
        new_events = sorted_t(pos(k):group_ends(k));
        n_new = numel(new_events);
        
        if n_new >= K_buf_backfill
            % If we have more new events than the buffer, just take the latest K
            pixel_history(p_idx, :) = new_events(end-K_buf_backfill+1:end)';
        else
            % Shift old history left, insert new events on the right
            pixel_history(p_idx, 1:(end-n_new)) = pixel_history(p_idx, (n_new+1):end);
            pixel_history(p_idx, (end-n_new+1):end) = new_events';
        end
    end

    % ---------------------------------------------------------------------
    % Step 2: Define the Spatial Trigger Zone (Vectorized)
    % ---------------------------------------------------------------------
    active_mask = false(imgSz);
    active_idx = unique_idx(counts(unique_idx) > 0);
    active_mask(active_idx) = true;

    se = strel('square', kernel_size); 
    trigger_zone_mask = imdilate(active_mask, se);
    target_pixels = find(trigger_zone_mask);

    % ---------------------------------------------------------------------
    % Step 3: Vectorized Statistics Calculation
    % ---------------------------------------------------------------------
    if isempty(target_pixels)
        return;
    end

    % Extract the history buffer strictly for the triggered pixels
    % size: [N_target_pixels, K_buf_backfill]
    target_history = pixel_history(target_pixels, :); 

    % Compute absolute time statistics across the rows (dim 2)
    % 'omitnan' gracefully handles pixels that haven't reached K_buf events yet
    means = mean(target_history, 2, 'omitnan');
    stds  = std(target_history, 0, 2, 'omitnan');
    
    % Compute valid counts per row
    valid_counts = sum(~isnan(target_history), 2);

    % Write back absolute time stats
    bf_t_mean(target_pixels) = means;
    bf_t_std(target_pixels)  = stds;
    counts(target_pixels)    = valid_counts; % Overwrite with total history count

    % Compute IEI (diff) across rows
    % diff on NaNs results in NaNs, which 'omitnan' will handle
    iei_matrix = diff(target_history, 1, 2);
    
    mean_diffs = mean(iei_matrix, 2, 'omitnan');
    std_diffs  = std(iei_matrix, 0, 2, 'omitnan');

    % Write back IEI stats
    bf_t_mean_diff(target_pixels) = mean_diffs;
    bf_t_std_diff(target_pixels)  = std_diffs;
    
    % Zero out any NaNs that resulted from rows with < 2 events
    bf_t_mean_diff(isnan(bf_t_mean_diff)) = 0;
    bf_t_std_diff(isnan(bf_t_std_diff)) = 0;
    bf_t_mean(isnan(bf_t_mean)) = 0;
    bf_t_std(isnan(bf_t_std)) = 0;

end
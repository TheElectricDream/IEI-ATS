function [bf_t_mean_diff,bf_t_std_diff, counts] = backfillIEIStatistics(sorted_t, unique_idx, pos, ...
    group_ends, imgSz, counts, K_buf_backfill, event_storage)
    
    % Calculate the window statistics for the current frame
    % [t_mean, t_std, t_max, t_min, t_mean_diff, t_std_diff] = ...
    %     stats.computeNeighborhoodStats(sorted_t, unique_idx, pos, ...
    %     group_ends, imgSz);

    % Based on the counts map, determine which pixels have fewer then K_buf
    % but greater then 1
    lowEventCount = find((counts >= 1) & (counts < K_buf_backfill));
    
    % Event storage contains the last M window slices of events, so we loop
    % through them and grab events at the pixel locations identified in
    % lowCountPixels. We only grab events until counts are are greater then
    % K_buf_backfill. Initialize a cell array to hold the backfilled history 
    % for each lowEventCount pixel
    bf_sorted_t   = sorted_t(:); 
    bf_unique_idx = unique_idx(:);
    bf_group_ends = group_ends(:);
    bf_pos        = pos(:);         % <--- NEW: Track the start indices
    
    for p = 1:length(lowEventCount)
        target_pixel = lowEventCount(p);
        pixel_history = [];
        
        % Search backward through history (skip index 1 since it's the current frame)
        for f = 2:length(event_storage)
            frame_data = event_storage{f};
            if isempty(frame_data); continue; end
            
            uidx = find(frame_data.unique_idx == target_pixel, 1);
            if ~isempty(uidx)
                if uidx == 1
                    start_bound = 1;
                else
                    start_bound = frame_data.group_ends(uidx - 1) + 1;
                end
                end_bound   = frame_data.group_ends(uidx);
                
                frame_pixel_events = frame_data.sorted_t(start_bound:end_bound);
                pixel_history = [frame_pixel_events(:).', pixel_history]; %#ok<AGROW>
            end
            
            % Check if historical + current events meet the required K_buf
            % (First, find how many events the current frame already has for this pixel)
            curr_uidx = find(bf_unique_idx == target_pixel, 1);
            
            % SAFE INDEXING: Use an explicit if/else block to prevent evaluating index 0
            if curr_uidx == 1
                curr_count = bf_group_ends(1);
            else
                curr_count = bf_group_ends(curr_uidx) - bf_group_ends(curr_uidx-1);
            end
            
            if (length(pixel_history) + curr_count) >= K_buf_backfill
                % Trim history so total combined events equals exactly K_buf
                needed = K_buf_backfill - curr_count;
                pixel_history = pixel_history(end - needed + 1:end);
                break;
            end
        end
        
        % If we found historical events, inject them into the vectors
        if ~isempty(pixel_history)
            % Recalculate c_start/c_end
            if curr_uidx == 1
                c_start = 1;
            else
                c_start = bf_group_ends(curr_uidx - 1) + 1;
            end
            c_end = bf_group_ends(curr_uidx);
            
            updated_slice = [pixel_history(:); bf_sorted_t(c_start:c_end)];
            num_added = length(pixel_history);
            
            % Vertically stack the pieces
            bf_sorted_t = [bf_sorted_t(1:c_start-1); updated_slice; bf_sorted_t(c_end+1:end)];
            
            % Shift all subsequent group ends forward
            bf_group_ends(curr_uidx:end) = bf_group_ends(curr_uidx:end) + num_added;
            
            % ---> NEW: Shift all subsequent group starts (pos) forward <---
            if curr_uidx < length(bf_pos)
                bf_pos(curr_uidx+1:end) = bf_pos(curr_uidx+1:end) + num_added;
            end
            
            % ---> NEW: Update the global counts array for downstream <---
            counts(target_pixel) = counts(target_pixel) + num_added;
        end
    end

    % Calculate the window statistics using the updated backfilled arrays
    [bf_t_mean, bf_t_std, bf_t_max, bf_t_min, bf_t_mean_diff, bf_t_std_diff] = ...
        stats.computeNeighborhoodStats(bf_sorted_t, bf_unique_idx, bf_pos, ...
        bf_group_ends, imgSz);
end



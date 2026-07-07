function [bf_t_mean, bf_t_std, bf_t_mean_diff, bf_t_std_diff, counts] = ...
    backfillIEIStatisticsOptimized(sorted_t, unique_idx, pos, group_ends, imgSz, ...
        counts, K_buf_backfill, event_storage, kernel_size, ...
        t_mean, t_std, t_mean_diff, t_std_diff)
    % BACKFILLIEISTATISTICSOPTIMIZED Spatial-temporal backfill of event history
    % for low-count pixels, with incremental statistics update.
    % 
    % [bf_t_mean, bf_t_std, bf_t_mean_diff, bf_t_std_diff, counts] = ...
    %     backfillIEIStatisticsOptimized(sorted_t, unique_idx, pos, group_ends, imgSz, ...
    %     counts, K_buf_backfill, event_storage, kernel_size, ...
    %     t_mean, t_std, t_mean_diff, t_std_diff)
    % 
    % Inputs:
    %   sorted_t       - [M x 1] Current-frame timestamps (ABSOLUTE time),
    %                    grouped by ascending linear pixel index, ascending in
    %                    time within each group
    %   unique_idx     - [K x 1] Unique linear pixel index per group, ascending
    %   pos            - [K x 1] Group start indices into sorted_t
    %   group_ends     - [K x 1] Group end indices into sorted_t
    %   imgSz          - [1 x 2] Image dimensions [nRows, nCols]
    %   counts         - [imgSz] Per-pixel counts (may be pre-masked; masked
    %                    pixels are never backfilled and never resurrected)
    %   K_buf_backfill - Target minimum events per active pixel
    %   event_storage  - {1 x F} cell of structs with fields sorted_t
    %                    (ABSOLUTE time), unique_idx, group_ends, and
    %                    optionally evt_idx ([events x 1] per-event linear
    %                    pixel index, cached by the producer). Element 1 =
    %                    current frame (skipped); 2..F progressively older
    %   kernel_size    - Odd spatial search window (1 = same pixel only)
    %   t_mean, t_std, t_mean_diff, t_std_diff
    %                  - [imgSz] current-frame statistics maps from
    %                    stats.computeNeighborhoodStats (main.m line ~201).
    %                    Returned unchanged except at pixels that receive
    %                    backfilled events

    arguments
        sorted_t       (:,1) double
        unique_idx     (:,1) double {mustBeInteger, mustBePositive}
        pos            (:,1) double {mustBeInteger, mustBePositive}
        group_ends     (:,1) double {mustBeInteger, mustBePositive}
        imgSz          (1,2) double {mustBeInteger, mustBePositive}
        counts         (:,:) double
        K_buf_backfill (1,1) double {mustBeInteger, mustBePositive}
        event_storage  (1,:) cell
        kernel_size    (1,1) double {mustBeInteger, mustBePositive}
        t_mean         (:,:) double
        t_std          (:,:) double
        t_mean_diff    (:,:) double
        t_std_diff     (:,:) double
    end

    assert(mod(kernel_size, 2) == 1, 'kernel_size must be odd');
    assert(isequal(size(counts), [imgSz(1), imgSz(2)]), ...
        'counts must match imgSz');

    % Pass-through defaults: untouched pixels keep the caller's statistics
    bf_t_mean      = t_mean;
    bf_t_std       = t_std;
    bf_t_mean_diff = t_mean_diff;
    bf_t_std_diff  = t_std_diff;

    % Backfill only pixels that are alive (count >= 1; upstream masking
    % zeroed rejected pixels) but below the target count
    target_pixels = find(counts(:) >= 1 & counts(:) < K_buf_backfill);
    if isempty(target_pixels)
        return;   % nothing to do -- and nothing to recompute
    end

    rem_needed = zeros(prod(imgSz), 1);
    rem_needed(target_pixels) = K_buf_backfill - counts(target_pixels);

    num_frames    = numel(event_storage);
    kept_hist_idx = cell(num_frames, 1);
    kept_hist_t   = cell(num_frames, 1);

    if kernel_size > 1
        half_k   = floor(kernel_size / 2);
        [dr, dc] = ndgrid(-half_k:half_k, -half_k:half_k);
        dr = dr(:)'; dc = dc(:)';
    end

    % Walk history newest-to-oldest so the most recent events fill first
    for f = 2:num_frames
        if ~any(rem_needed(target_pixels) > 0), break; end

        frame = event_storage{f};
        if isempty(frame), continue; end

        % Per-event pixel indices: use the producer's cached copy if present
        if isfield(frame, 'evt_idx')
            f_idx = frame.evt_idx;
        else
            f_idx = repelem(frame.unique_idx, diff([0; frame.group_ends]));
        end
        f_t = frame.sorted_t;   % ABSOLUTE time (see contract)

        % PRUNE FIRST: one O(events) lookup against the dilated needy-pixel
        % mask replaces expanding every event by kernel_size^2. conv2 of a
        % binary mask with a ones kernel is standard binary dilation with a
        % square structuring element (Gonzalez & Woods, Digital Image
        % Processing); recomputed per frame so it tightens as needs are met
        need_map = reshape(rem_needed > 0, imgSz);
        if kernel_size > 1
            reach = conv2(double(need_map), ones(kernel_size), 'same') > 0;
        else
            reach = need_map;
        end
        keep_near_target = reach(f_idx);
        f_idx = f_idx(keep_near_target);
        f_t   = f_t(keep_near_target);
        if isempty(f_idx), continue; end

        % Kernel expansion -- survivors only
        if kernel_size > 1
            [r, c] = ind2sub(imgSz, f_idx);
            r_exp  = r(:) + dr;
            c_exp  = c(:) + dc;
            t_exp  = repmat(f_t(:), 1, kernel_size^2);

            valid = r_exp >= 1 & r_exp <= imgSz(1) & ...
                    c_exp >= 1 & c_exp <= imgSz(2);
            f_idx = sub2ind(imgSz, r_exp(valid), c_exp(valid));
            f_t   = t_exp(valid);
        end

        % Keep only copies that land on pixels still in need
        valid_mask = rem_needed(f_idx) > 0;
        f_idx = f_idx(valid_mask);
        f_t   = f_t(valid_mask);
        if isempty(f_idx), continue; end

        if kernel_size > 1
            % Deduplicate exact (pixel, timestamp) pairs from the expansion,
            % which would otherwise inject spurious zero IEIs
            dedup = unique([f_idx, f_t], 'rows');
            f_idx = dedup(:, 1);
            f_t   = dedup(:, 2);
        end

        % Per pixel, keep the most recent rem_needed events
        % (sort by pixel, then time descending; small arrays after pruning)
        [~, sort_order] = sortrows([f_idx, -f_t]);
        f_idx = f_idx(sort_order);
        f_t   = f_t(sort_order);

        [u_idx, start_pos] = unique(f_idx, 'first');
        [~,     end_pos]   = unique(f_idx, 'last');

        keep_counts = min(end_pos - start_pos + 1, rem_needed(u_idx));
        keep_mask   = false(size(f_idx));
        for k = 1:max(keep_counts)
            valid_groups = keep_counts >= k;
            keep_mask(start_pos(valid_groups) + k - 1) = true;
        end

        kept_hist_idx{f} = f_idx(keep_mask);
        kept_hist_t{f}   = f_t(keep_mask);

        rem_needed(u_idx) = rem_needed(u_idx) - keep_counts;
    end

    kept_hist_idx = vertcat(kept_hist_idx{:});
    kept_hist_t   = vertcat(kept_hist_t{:});
    if isempty(kept_hist_idx)
        return;   % history had nothing to give; maps pass through
    end

    % ---- INCREMENTAL UPDATE: recompute only pixels that received events --
    recv_pixels = unique(kept_hist_idx);              % sorted receiving pixels
    [recv_found, recv_group_k] = ismember(recv_pixels, unique_idx);
    assert(all(recv_found), ...
        'backfill produced events at pixels outside unique_idx');

    % Gather the receiving pixels' current-frame events with a vectorized
    % run expansion (concatenated pos(k):group_ends(k) index runs)
    recv_starts = pos(recv_group_k);
    recv_ends   = group_ends(recv_group_k);
    recv_counts = recv_ends - recv_starts + 1;

    gather_idx  = ones(sum(recv_counts), 1);
    run_heads   = cumsum([1; recv_counts(1:end-1)]);
    gather_idx(run_heads) = [recv_starts(1); ...
        recv_starts(2:end) - recv_ends(1:end-1)];
    gather_idx  = cumsum(gather_idx);

    current_recv_idx = repelem(recv_pixels, recv_counts);
    current_recv_t   = sorted_t(gather_idx);

    % Merge current + backfilled events for the receiving pixels only
    merged_idx = [current_recv_idx; kept_hist_idx];
    merged_t   = [current_recv_t;   kept_hist_t];
    [~, merge_order] = sortrows([merged_idx, merged_t]);
    merged_idx = merged_idx(merge_order);
    merged_t   = merged_t(merge_order);

    [~, merged_pos]        = unique(merged_idx, 'first');
    [~, merged_group_ends] = unique(merged_idx, 'last');

    % Same frame-local origin as the caller's maps (see header note)
    merged_t = merged_t - min(sorted_t);

    % Per-pixel statistics on the merged sub-stream; recv_pixels /
    % merged_pos / merged_group_ends form a valid grouping triple, so the
    % canonical stats function applies directly (its maps are zero away
    % from recv_pixels, which we never read)
    [recv_t_mean, recv_t_std, ~, ~, recv_t_mean_diff, recv_t_std_diff] = ...
        stats.computeNeighborhoodStats(merged_t, recv_pixels, ...
        merged_pos, merged_group_ends, imgSz);

    bf_t_mean(recv_pixels)      = recv_t_mean(recv_pixels);
    bf_t_std(recv_pixels)       = recv_t_std(recv_pixels);
    bf_t_mean_diff(recv_pixels) = recv_t_mean_diff(recv_pixels);
    bf_t_std_diff(recv_pixels)  = recv_t_std_diff(recv_pixels);
    counts(recv_pixels)         = merged_group_ends - merged_pos + 1;
end
function [bf_t_mean,bf_t_mean_diff, bf_t_std_diff, counts] = backfillIEIStatisticsOptimized(sorted_t, unique_idx, pos, ...
    group_ends, imgSz, counts, K_buf_backfill, event_storage)
% BACKFILLIEISTATISTICS Top up low-event-count pixels with historical events, then compute windowed IEI statistics.
%   Drop-in, output-identical replacement for the reference implementation.
%   History lookups use O(1) scatter LUTs (no per-pixel linear scans) and the
%   backfilled timestamp vector is assembled in a single pre-allocated pass
%   instead of being reallocated for every backfilled pixel.
%
%   Inputs:
%     sorted_t        - (N x 1)  Timestamps, grouped per pixel; ascending within each group.
%     unique_idx      - (U x 1)  Linear pixel indices present in the current frame.
%     pos             - (U x 1)  Group start indices. Retained for signature compatibility; unused.
%     group_ends      - (U x 1)  Cumulative group end indices into sorted_t.
%     imgSz           - (1 x 2)  Sensor resolution [rows cols].
%     counts          - (imgSz)  Per-pixel current event counts (any shape; read column-major).
%     K_buf_backfill  - (scalar) Target minimum event count per pixel.
%     event_storage   - (1 x F cell) Frame ring buffer. Index 1 is the current frame;
%                        2..F are progressively older frames. Each cell has fields
%                        .sorted_t, .unique_idx, .group_ends with the same layout as above.
%
%   Outputs:
%     bf_t_mean_diff  - (U x 1)  Mean inter-event interval over the backfilled window.
%     bf_t_std_diff   - (U x 1)  Std of the inter-event interval over the backfilled window.
%
%   Notes:
%     - Bit-for-bit identical to the previous implementation. For each low-count
%       pixel the backfill set is the newest `needed` historical events (newest
%       frame first, most-recent tail within the frame that crosses the
%       threshold), merged with the current-frame slice and re-sorted ascending.
%     - Pixels not found in history, or already at/above K_buf_backfill, are
%       left unchanged.
%     - Validated against the reference via a randomized differential test
%       (image size, K, frame count, sparsity, empty/absent-history frames,
%       adversarial timestamp overlap); arrays passed to
%       stats.computeNeighborhoodStats match exactly.
%
%   See also stats.computeNeighborhoodStats

    % ----------------------------------------------------------------------
    % 1. Setup
    % ----------------------------------------------------------------------
    U      = numel(unique_idx);
    nPix   = prod(imgSz);
    starts = [1; group_ends(1:end-1) + 1];          % current-frame group starts

    % ----------------------------------------------------------------------
    % 2. Identify low-count pixels that exist in the current frame
    % ----------------------------------------------------------------------
    cur_lut             = zeros(nPix, 1);
    cur_lut(unique_idx) = 1:U;                       % pixel index -> position in unique_idx

    lowPix  = find((counts(:) >= 1) & (counts(:) < K_buf_backfill));
    lowUidx = cur_lut(lowPix);
    keep    = lowUidx > 0;                           % drop any not present this frame
    lowPix  = lowPix(keep);
    lowUidx = lowUidx(keep);

    cur_counts = group_ends(lowUidx) - starts(lowUidx) + 1;
    needed     = K_buf_backfill - cur_counts;        % events still required per low pixel
    numLow     = numel(lowPix);

    collected  = cell(numLow, 1);                    % historical events gathered per low pixel
    rem_needed = needed;                             % remaining requirement (mutated below)

    % ----------------------------------------------------------------------
    % 3. Gather newest historical events (newest frames first)
    % ----------------------------------------------------------------------
    frame_lut = zeros(nPix, 1);                      % reusable scatter LUT (touched entries reset each frame)
    for f = 2:numel(event_storage)
        if ~any(rem_needed > 0); break; end
        fd = event_storage{f};
        if isempty(fd); continue; end

        frame_lut(fd.unique_idx) = 1:numel(fd.unique_idx);
        f_starts = [1; fd.group_ends(1:end-1) + 1];

        active = find(rem_needed > 0);
        fu     = frame_lut(lowPix(active));          % position in this frame (0 if absent)
        hit    = fu > 0;
        active = active(hit);
        fu     = fu(hit);

        for k = 1:numel(active)
            p     = active(k);
            u     = fu(k);
            slice = fd.sorted_t(f_starts(u):fd.group_ends(u));   % ascending in time
            take  = rem_needed(p);
            if numel(slice) > take
                slice = slice(end-take+1:end);       % keep the most recent `take`
            end
            collected{p}  = [collected{p}; slice(:)];
            rem_needed(p) = rem_needed(p) - numel(slice);
        end

        frame_lut(fd.unique_idx) = 0;                % reset touched entries for reuse
    end

    % ----------------------------------------------------------------------
    % 4. Assemble the backfilled timestamp vector in one pre-allocated pass
    % ----------------------------------------------------------------------
    nAdded_per_low = cellfun(@numel, collected);            % events added per low pixel
    counts(lowPix) = counts(lowPix) + nAdded_per_low;       % <-- your update

    add           = zeros(U, 1);
    add(lowUidx)  = nAdded_per_low;
    bf_group_ends = group_ends + cumsum(add);
    bf_pos        = [1; bf_group_ends(1:end-1) + 1];
    N_new         = bf_group_ends(end);

    bf_sorted_t = zeros(N_new, 1, 'like', sorted_t);

    % Scatter the original events to their shifted positions (correct for every
    % unchanged group; backfilled groups are fully overwritten just below).
    offset_before = [0; cumsum(add(1:end-1))];               % additions preceding each group
    old_counts    = group_ends - starts + 1;
    event_offset  = repelem(offset_before, old_counts);      % per-event shift
    bf_sorted_t((1:numel(sorted_t)).' + event_offset) = sorted_t;

    % Overwrite the slices that received history with the merged, sorted result.
    for k = 1:numLow
        h = collected{k};
        if isempty(h); continue; end
        u      = lowUidx(k);
        merged = sort([h; sorted_t(starts(u):group_ends(u))]);
        bf_sorted_t(bf_pos(u):bf_group_ends(u)) = merged;
    end

    % ----------------------------------------------------------------------
    % 5. Windowed statistics on the backfilled arrays
    % ----------------------------------------------------------------------
    [bf_t_mean, ~, ~, ~, bf_t_mean_diff, bf_t_std_diff] = ...
        stats.computeNeighborhoodStats(bf_sorted_t, unique_idx, bf_pos, ...
        bf_group_ends, imgSz);
end
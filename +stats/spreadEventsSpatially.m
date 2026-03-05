function [sorted_x, sorted_y, sorted_t, unique_idx, ...
    pos, group_ends] = ...
    spreadEventsSpatially(x, y, t, imgSz, radius)
% SPREADEVENTSSPATIALLY  Spatial event propagation via kernel replication.
%
%   [SORTED_X, SORTED_Y, SORTED_T, UNIQUE_IDX, POS, GROUP_ENDS] =
%   SPREADEVENTSSPATIALLY(X, Y, T, IMGSZ, RADIUS) replicates each
%   event to its spatial neighbours within a (2R+1) x (2R+1) kernel,
%   producing an expanded event set suitable for neighbourhood
%   statistics computation.
%
%   This effectively performs a spatial convolution on the event
%   stream: replicated events at each pixel location allow
%   computeNeighborhoodStats to compute local spatial averages and
%   local IEI coherence measures.
%
%   Inputs:
%     x, y   - [N x 1] Coordinate vectors (row, col). Converted
%              to int32 internally if needed.
%     t      - [N x 1] Timestamp vector [s].
%     imgSz  - [1 x 2] Image dimensions [nRows, nCols].
%     radius - Scalar integer expansion radius. A radius of 2
%              produces a 5x5 kernel (-2:+2 offsets).
%
%   Outputs:
%     sorted_x, sorted_y, sorted_t - Expanded and sorted event
%              vectors. Sorted by linear pixel index, then by
%              timestamp within each pixel group.
%     unique_idx - [K x 1] Linear pixel indices for each active
%              pixel in the expanded set.
%     pos        - [K x 1] Start index of each group in sorted_*.
%     group_ends - [K x 1] End index of each group in sorted_*.
%
%   Algorithm:
%     1. Build a (2R+1)^2 offset kernel.
%     2. Replicate each event K times (one per kernel offset).
%     3. Filter out-of-bounds events.
%     4. Sort by (pixel index, timestamp) and compute grouping
%        indices for downstream stats computation.
%
%   Notes:
%     - Output format is designed to feed directly into
%       stats.computeNeighborhoodStats.
%     - Memory usage scales as O(N * (2R+1)^2). For large N and
%       R, this can be significant.
%     - Coordinates: x = row, y = col, sub2ind(imgSz, x, y).
%
%   See also: stats.computeNeighborhoodStats

    % ----------------------------------------------------------------
    % 1. Define kernel offsets
    % ----------------------------------------------------------------
    range = -radius:radius;
    [dY, dX] = meshgrid(range, range);
    off_x = int32(dX(:));
    off_y = int32(dY(:));
    num_replicas = length(off_x);

    % ----------------------------------------------------------------
    % 2. Replicate events (vectorized)
    % ----------------------------------------------------------------
    n_events = length(x);

    if ~isa(x, 'int32'), x = int32(x); end
    if ~isa(y, 'int32'), y = int32(y); end

    % Index map: replicate each event num_replicas times
    idx_map = floor((0:n_events*num_replicas - 1)' ...
        / num_replicas) + 1;

    x_big = x(idx_map);
    y_big = y(idx_map);
    t_big = t(idx_map);

    % Apply spatial offsets
    off_x_big = repmat(off_x, n_events, 1);
    off_y_big = repmat(off_y, n_events, 1);

    x_new = x_big + off_x_big;
    y_new = y_big + off_y_big;

    % ----------------------------------------------------------------
    % 3. Filter out-of-bounds events
    % ----------------------------------------------------------------
    valid_mask = x_new >= 1 & x_new <= imgSz(1) & ...
                 y_new >= 1 & y_new <= imgSz(2);

    x_valid = x_new(valid_mask);
    y_valid = y_new(valid_mask);
    t_valid = t_big(valid_mask);

    % ----------------------------------------------------------------
    % 4. Sort by pixel index and timestamp, compute grouping
    % ----------------------------------------------------------------
    linear_idx = sub2ind(imgSz, x_valid, y_valid);

    % Multi-column sort: pixel ID first, then timestamp
    [~, sort_order] = sortrows(...
        [double(linear_idx), double(t_valid)]);

    sorted_x = x_valid(sort_order);
    sorted_y = y_valid(sort_order);
    sorted_t = t_valid(sort_order);
    sorted_idx = linear_idx(sort_order);

    % Grouping indices for stats computation
    [unique_idx, pos, ~] = unique(sorted_idx);
    group_ends = [pos(2:end) - 1; length(sorted_idx)];

end
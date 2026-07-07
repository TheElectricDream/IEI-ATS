function [t_mean, t_std, t_max, t_min, t_mean_diff, t_std_diff, fcn_time] = ...
    computeNeighborhoodStats(t_sorted, unique_idx, pos, ...
    group_ends, img_size)
    % COMPUTENEIGHBORHOODSTATS Per-pixel timestamp and inter-event interval (IEI)
    % statistics from grouped event streams.
    % 
    % [t_mean, t_std, t_max, t_min, t_mean_diff, t_std_diff] = ...
    %     computeNeighborhoodStats(t_sorted, unique_idx, pos, ...
    %     group_ends, img_size)
    % 
    % Inputs:
    %   t_sorted   - [M x 1] Event timestamps grouped by linear pixel
    %                index and sorted ascending within each group
    %   linear_idx - [M x 1] Linear pixel index of each group
    %   unique_idx - [K x 1] Unique pixel index of each group
    %   pos        - [K x 1] Start index of each group in t_sorted
    %                Must be ascending; groups must contiguously
    %                partition t_sorted
    %   group_ends - [K x 1] End index of each group in t_sorted
    %   img_size   - [1 x 2] Image dimensions [nRows, nCols]
    % 
    % Outputs (all [img_size] double maps; zero at pixels with no events):
    %   t_mean      - Mean timestamp
    %   t_std       - Sample (Bessel-corrected) standard deviation of
    %                 timestamps. Zero for single-event pixels
    %   t_max       - Maximum timestamp
    %   t_min       - Minimum timestamp
    %   t_mean_diff - Mean IEI. 
    %   t_std_diff  - Sample standard deviation of IEIs. 

    arguments
        t_sorted   (:,1) double
        unique_idx (:,1) double {mustBeInteger, mustBePositive}
        pos        (:,1) double {mustBeInteger, mustBePositive}
        group_ends (:,1) double {mustBeInteger, mustBePositive}
        img_size   (1,2) double {mustBeInteger, mustBePositive}
    end

    % Assertions
    assert(numel(pos) == numel(unique_idx) && ...
           numel(pos) == numel(group_ends), ...
        'pos, group_ends, and unique_idx must have equal length');
    assert(isempty(pos) || issorted(pos, 'strictascend'), ...
        'pos must be strictly ascending');
    assert(all(group_ends >= pos), ...
        'each group_ends(k) must be >= pos(k)');

    % Start internal timer
    tic;

    % Initialize output maps
    t_mean      = zeros(img_size);
    t_std       = zeros(img_size);
    t_max       = zeros(img_size);
    t_min       = zeros(img_size);
    t_mean_diff = zeros(img_size);
    t_std_diff  = zeros(img_size);

    K = numel(unique_idx);
    if K == 0
        return;
    end

    % Get the minimum & maximum time for each grouping (x,y)
    minimum_sorted_t = t_sorted(pos);
    maximum_sorted_t = t_sorted(group_ends);

    % Number of events and inter-event intervals in each grouping
    event_counts = group_ends - pos + 1;
    iei_counts   = event_counts - 1;
 
    % Remove each grouping's first timestamp so the shifted values are
    % non-negative and bounded 
    t_sorted_shifted    = t_sorted - repelem(minimum_sorted_t, event_counts);
    cumsum_t_shifted    = [0; cumsum(t_sorted_shifted)];
    cumsum_t_shifted_sq = [0; cumsum(t_sorted_shifted.^2)];
 
    event_group_t_sums  = cumsum_t_shifted(group_ends + 1) ...
                        - cumsum_t_shifted(pos);
    event_group_t_sumsq = cumsum_t_shifted_sq(group_ends + 1) ...
                        - cumsum_t_shifted_sq(pos);
 
    % Calculate the mean time for each grouping (x,y)
    event_group_means = event_group_t_sums ./ event_counts + minimum_sorted_t;
 
    % Bessel-corrected variance of the timestamps, and zero for single-event
    % groupings
    event_group_t_vars = (event_group_t_sumsq ...
        - event_group_t_sums.^2 ./ event_counts) ./ max(event_counts - 1, 1);
    event_group_t_vars(event_counts == 1) = 0;
    event_group_stds = sqrt(max(event_group_t_vars, 0));
 
    % Within-group diffs of grouping k
    t_sorted_diffs = diff(t_sorted);
 
    % Calculate the telescoping sum of within-group IEIs
    event_group_iei_sums = maximum_sorted_t - minimum_sorted_t;
 
    cumsum_iei_sq = [0; cumsum(t_sorted_diffs.^2)];
    event_group_iei_sumsq = cumsum_iei_sq(group_ends) - cumsum_iei_sq(pos);
 
    % Calculate the mean IEI
    event_group_mean_diffs = zeros(K, 1);
    has_intervals = iei_counts >= 1;
    event_group_mean_diffs(has_intervals) = ...
        event_group_iei_sums(has_intervals) ./ iei_counts(has_intervals);
 
    % Bessel-corrected standard deviation of the IEIs, and zero for 
    % groupings with fewer then two points
    event_group_iei_vars = zeros(K, 1);
    has_multiple_intervals = iei_counts >= 2;
    event_group_iei_vars(has_multiple_intervals) = ...
        (event_group_iei_sumsq(has_multiple_intervals) ...
        - event_group_iei_sums(has_multiple_intervals).^2 ...
        ./ iei_counts(has_multiple_intervals)) ...
        ./ (iei_counts(has_multiple_intervals) - 1);
    event_group_std_diffs = sqrt(max(event_group_iei_vars, 0));

    % Scatter into image maps
    t_mean(unique_idx)      = event_group_means;
    t_std(unique_idx)       = event_group_stds;
    t_max(unique_idx)       = maximum_sorted_t;
    t_min(unique_idx)       = minimum_sorted_t;
    t_mean_diff(unique_idx) = event_group_mean_diffs;
    t_std_diff(unique_idx)  = event_group_std_diffs;

    % Get total function time
    fcn_time = toc; 

end
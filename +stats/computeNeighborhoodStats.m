function [t_mean, t_std, t_max, t_min, t_mean_diff, t_std_diff] = ...
    computeNeighborhoodStats(sorted_t, unique_idx, pos, ...
    group_ends, imgSz)
% COMPUTENEIGHBORHOODSTATS Per-pixel timestamp and inter-event interval (IEI) statistics from grouped event streams.
%
%     [t_mean, t_std, t_max, t_min, t_mean_diff, t_std_diff] = ...
%         computeNeighborhoodStats(sorted_t, unique_idx, pos, ...
%         group_ends, imgSz)
%
%     Inputs:
%       sorted_t   - [M x 1] Event timestamps grouped by linear pixel
%                    index and sorted ascending within each group.
%       unique_idx - [K x 1] Linear pixel index of each group.
%       pos        - [K x 1] Start index of each group in sorted_t.
%                    Must be ascending; groups must contiguously
%                    partition sorted_t.
%       group_ends - [K x 1] End index of each group in sorted_t.
%       imgSz      - [1 x 2] Image dimensions [nRows, nCols].
%
%     Outputs (all [imgSz] double maps; zero at pixels with no events):
%       t_mean      - Mean timestamp.
%       t_std       - Sample (Bessel-corrected) standard deviation of
%                     timestamps. Zero for single-event pixels.
%       t_max       - Maximum timestamp.
%       t_min       - Minimum timestamp.
%       t_mean_diff - Mean IEI. Zero for single-event pixels.
%       t_std_diff  - Sample standard deviation of IEIs. Zero for
%                     pixels with fewer than three events (fewer than
%                     two IEIs).
%
%     Notes:
%       Fully vectorized; no per-pixel loops. Min/max are read from
%       group boundaries; means use prefix-sum (cumsum) subtraction;
%       standard deviations use the one-pass identity
%       var(x) = E[x^2] - E[x]^2 with Bessel correction; the IEI mean
%       uses the telescoping identity
%       mean(diff(x)) = (x(end) - x(1)) / (n - 1); IEI squared sums are
%       scattered per group with accumarray after masking diffs that
%       cross group boundaries.
%
%       The one-pass variance formula can suffer catastrophic
%       cancellation when magnitudes are large relative to spread
%       (Higham 2002, "Accuracy and Stability of Numerical Algorithms,"
%       2nd ed., Ch. 1). Timestamps here are window-relative, keeping
%       magnitudes moderate; a max(., 0) clamp guards residual
%       floating-point noise.
%
%       t_mean_diff feeds the coherence IEI regularity rule and, via
%       EMA smoothing, the IEI-ATS per-pixel tau mapping.
%
%     See also COHERENCE.COMPUTECOHERENCEMASK.

    % ----------------------------------------------------------------
    % 0. Initialize output maps
    % ----------------------------------------------------------------
    t_mean      = zeros(imgSz);
    t_std       = zeros(imgSz);
    t_max       = zeros(imgSz);
    t_min       = zeros(imgSz);
    t_mean_diff = zeros(imgSz);
    t_std_diff  = zeros(imgSz);

    K = numel(unique_idx);
    if K == 0
        return;
    end

    % ----------------------------------------------------------------
    % 1. Group counts
    % ----------------------------------------------------------------
    counts = group_ends - pos + 1;         % [K x 1] events per pixel

    % ----------------------------------------------------------------
    % 2. Min / max — direct indexing (sorted within group)
    % ----------------------------------------------------------------
    mn = sorted_t(pos);
    mx = sorted_t(group_ends);

    % ----------------------------------------------------------------
    % 3. Mean — prefix-sum subtraction
    % ----------------------------------------------------------------
    % Prepend zero so that cs(end+1) - cs(start) gives the group sum
    % without special-casing the first group.
    cs      = [0; cumsum(sorted_t)];
    grp_sum = cs(group_ends + 1) - cs(pos);
    mu      = grp_sum ./ counts;

    % ----------------------------------------------------------------
    % 4. Std — prefix-sum of squares + algebraic variance
    % ----------------------------------------------------------------
    %   var(x) = [ sum(x^2) - n * mu^2 ] / (n - 1)
    cs2      = [0; cumsum(sorted_t .^ 2)];
    grp_sum2 = cs2(group_ends + 1) - cs2(pos);

    denom    = max(counts - 1, 1);         % protect single-event pixels
    variance = (grp_sum2 - counts .* mu.^2) ./ denom;
    sd       = sqrt(max(variance, 0));     % clamp float noise
    sd(counts == 1) = 0;

    % ----------------------------------------------------------------
    % 5. IEI mean — telescoping sum identity
    % ----------------------------------------------------------------
    %   mean(diff(x)) = (x(end) - x(1)) / (numel(x) - 1)
    diff_counts = counts - 1;
    mu_d = (mx - mn) ./ max(diff_counts, 1);
    mu_d(counts <= 1) = 0;

    % ----------------------------------------------------------------
    % 6. IEI std — vectorized diff with cross-group masking
    % ----------------------------------------------------------------
    M     = numel(sorted_t);
    d_all = diff(sorted_t);               % [M-1 x 1]

    % 6a. Build per-element group labels via cumsum delta trick:
    %     place a 1 at every group start, cumsum maps each element
    %     to its group index.
    delta     = zeros(M, 1);
    delta(pos) = 1;
    g         = cumsum(delta);            % g(i) = group index of sorted_t(i)

    % 6b. Identify within-group diffs (exclude boundary diffs where
    %     consecutive elements belong to different groups).
    g_diff = g(1:end-1);                  % group label for each diff
    within = true(M - 1, 1);
    if K > 1
        within(group_ends(1:end-1)) = false;
    end

    % 6c. Scatter sum(d^2) per group via accumarray
    d_w    = d_all(within);
    gw     = g_diff(within);
    sum_d2 = accumarray(gw, d_w .^ 2, [K, 1]);

    %   var(d) = [ sum(d^2)/n - mu_d^2 ] * n/(n-1)    (Bessel)
    n_iei = max(diff_counts, 1);
    var_d = (sum_d2 ./ n_iei - mu_d .^ 2) .* n_iei ./ max(n_iei - 1, 1);
    sd_d  = sqrt(max(var_d, 0));
    sd_d(counts <= 2) = 0;               % need >= 2 IEIs for std

    % ----------------------------------------------------------------
    % 7. Scatter into image maps
    % ----------------------------------------------------------------
    t_mean(unique_idx)      = mu;
    t_std(unique_idx)       = sd;
    t_max(unique_idx)       = mx;
    t_min(unique_idx)       = mn;
    t_mean_diff(unique_idx) = mu_d;
    t_std_diff(unique_idx)  = sd_d;

end
function [bf_t_mean, bf_t_std, bf_t_mean_diff, bf_t_std_diff, counts] = ...
    backfillIEIStatisticsPooled(imgSz, counts, K_buf_backfill, map_storage, ...
        kernel_size, t_mean, t_std, t_mean_diff, t_std_diff)
% BACKFILLIEISTATISTICSPOOLED Pooled-neighborhood per-pixel mean IEI via
% constant-time-per-pixel map filters. Alternative code path to
% backfillIEIStatisticsOptimized for A/B comparison.
%
% ESTIMATOR: merging event streams from m pixels sums their rates
% (superposition of point processes; Cox & Smith, "On the superposition
% of renewal processes", Biometrika 41(1-2):91-99, 1954), so the raw
% merged mean IEI (t_max - t_min)/(n - 1) measures the reciprocal of the
% NEIGHBORHOOD AGGREGATE rate -- roughly (per-pixel IEI)/m -- not the
% per-pixel quantity IEI-ATS consumes. This function therefore
% renormalizes by the number of distinct contributing active pixels:
%
%     per-pixel-equivalent mean IEI = m * (t_max - t_min) / (n - 1)
%
% i.e., the reciprocal of the neighborhood-average per-pixel rate. All
% four ingredients (min of first timestamps, max of last timestamps, sum
% of counts, sum of activity indicators) are associative, so each is one
% separable box/min/max filter pass: van Herk (Pattern Recognition
% Letters 13(7):517-521, 1992) and Gil & Werman (IEEE TPAMI 15(5):
% 504-507, 1993) for min/max; summed-area tables (Crow, SIGGRAPH 1984)
% for the sums. O(N_pixels) per frame, independent of kernel size and
% event count.
%
% INTERPRETATION CAVEATS (state these in any writeup):
%   a. The estimate is the reciprocal of the neighborhood MEAN rate -- a
%      harmonic-type average of neighbor IEIs, <= their arithmetic
%      average when rates are heterogeneous. It is a rate regularizer
%      under a local-homogeneity assumption, not a measurement of the
%      pixel's own IEI.
%   b. IEI spread is NOT recoverable by pooling: superpositions of
%      renewal processes are not renewal except in the Poisson case, and
%      sparse superpositions tend to Poisson, so pooled interval
%      distributions wash toward exponential regardless of per-pixel
%      regularity. t_std_diff therefore passes through unchanged, as do
%      t_mean and t_std.
%
% SEMANTIC DIFFERENCES vs backfillIEIStatisticsOptimized:
%   1. Pools ALL history events in the window -- no top-K recency cap
%      (order statistics are not associative-decomposable). counts at
%      backfilled pixels may exceed K_buf_backfill, which here is purely
%      a sparsity threshold selecting WHERE to intervene.
%   2. Spatial borrowing applies to HISTORY frames only; the current
%      frame contributes each pixel's own events.
%   3. Note the exact version with kernel_size > 1 shares the
%      superposition deflation (bounded by contributing pixels within its
%      kernel) and does NOT renormalize; for a like-for-like A/B run it
%      with kernel_size = 1 (pure temporal extension, unbiased) or apply
%      the same m-renormalization there.
%
% Inputs:
%   imgSz          - [1 x 2] Image dimensions [nRows, nCols]
%   counts         - [imgSz] Per-pixel counts (may be pre-masked; masked
%                    pixels are never targeted)
%   K_buf_backfill - Sparsity threshold: pixels with
%                    1 <= counts < K_buf_backfill are backfilled
%   map_storage    - {1 x F} cell of structs with fields t_min, t_max
%                    (ABSOLUTE-time first/last event timestamp maps;
%                    values at zero-count pixels are ignored) and counts
%                    (raw, unfiltered). Element 1 = current frame;
%                    2..F progressively older
%   kernel_size    - Odd spatial pooling window (1 = same pixel only)
%   t_mean, t_std, t_mean_diff, t_std_diff
%                  - [imgSz] current-frame statistics maps. t_mean_diff
%                    is overwritten at backfilled pixels; the rest pass
%                    through (caveat b)
%
% Outputs:
%   bf_t_mean, bf_t_std, bf_t_mean_diff, bf_t_std_diff - maps (see above)
%   counts - pooled event counts at backfilled pixels; unchanged elsewhere

    arguments
        imgSz          (1,2) double {mustBeInteger, mustBePositive}
        counts         (:,:) double
        K_buf_backfill (1,1) double {mustBeInteger, mustBePositive}
        map_storage    (1,:) cell
        kernel_size    (1,1) double {mustBeInteger, mustBePositive}
        t_mean         (:,:) double
        t_std          (:,:) double
        t_mean_diff    (:,:) double
        t_std_diff     (:,:) double
    end

    assert(mod(kernel_size, 2) == 1, 'kernel_size must be odd');
    assert(isequal(size(counts), [imgSz(1), imgSz(2)]), ...
        'counts must match imgSz');
    assert(~isempty(map_storage) && ~isempty(map_storage{1}), ...
        'map_storage{1} (current frame) must be populated');

    % Pass-through defaults: untouched pixels keep the caller's statistics
    bf_t_mean      = t_mean;
    bf_t_std       = t_std;
    bf_t_mean_diff = t_mean_diff;
    bf_t_std_diff  = t_std_diff;

    % Same targeting rule as the exact version: alive but below threshold
    target_pixels = find(counts(:) >= 1 & counts(:) < K_buf_backfill);
    if isempty(target_pixels)
        return;
    end

    % ------------- Combine history frames elementwise (2..F) -------------
    % Empty pixels are neutralized with +/-Inf so they never win the
    % min/max; counts identify emptiness, so stored values there are moot
    hist_t_min  = inf(imgSz);
    hist_t_max  = -inf(imgSz);
    hist_counts = zeros(imgSz);

    for f = 2:numel(map_storage)
        frame = map_storage{f};
        if isempty(frame), continue; end

        frame_active = frame.counts > 0;
        frame_t_min  = frame.t_min;
        frame_t_max  = frame.t_max;
        frame_t_min(~frame_active) = Inf;
        frame_t_max(~frame_active) = -Inf;

        hist_t_min  = min(hist_t_min, frame_t_min);
        hist_t_max  = max(hist_t_max, frame_t_max);
        hist_counts = hist_counts + frame.counts;
    end

    % A history pixel contributes if it fired in ANY buffered frame; each
    % distinct pixel counts once toward the rate renormalization
    hist_active = hist_counts > 0;

    % --------- Neighborhood pooling of history (separable filters) -------
    if kernel_size > 1
        nb_t_min  = movmin(movmin(hist_t_min, kernel_size, 1), ...
                           kernel_size, 2);
        nb_t_max  = movmax(movmax(hist_t_max, kernel_size, 1), ...
                           kernel_size, 2);
        nb_counts = conv2(ones(kernel_size, 1), ones(1, kernel_size), ...
                          hist_counts, 'same');
        nb_active = conv2(ones(kernel_size, 1), ones(1, kernel_size), ...
                          double(hist_active), 'same');
    else
        nb_t_min  = hist_t_min;
        nb_t_max  = hist_t_max;
        nb_counts = hist_counts;
        nb_active = double(hist_active);
    end

    % ------- Merge with the current frame's own-pixel events only --------
    current = map_storage{1};
    current_active = current.counts > 0;
    current_t_min  = current.t_min;
    current_t_max  = current.t_max;
    current_t_min(~current_active) = Inf;
    current_t_max(~current_active) = -Inf;

    combined_t_min  = min(nb_t_min, current_t_min);
    combined_t_max  = max(nb_t_max, current_t_max);
    combined_counts = nb_counts + current.counts;

    % Distinct contributing pixels: all history-active pixels in the
    % window, plus the center pixel itself if it is only current-active
    % (a center pixel active in history is already counted in nb_active)
    contributing_pixels = nb_active + double(current_active & ~hist_active);

    % ------------- Renormalized telescoping estimator --------------------
    % Aggregate rate of the merged set = (n-1)/span; dividing by the
    % number of contributing pixels gives the neighborhood-average
    % per-pixel rate, whose reciprocal is the per-pixel-equivalent IEI
    pooled_iei_per_pixel = contributing_pixels ...
        .* (combined_t_max - combined_t_min) ...
        ./ max(combined_counts - 1, 1);

    % Write only at target pixels with at least one pooled interval;
    % targets with no history support keep their pass-through values
    write_mask = false(imgSz);
    write_mask(target_pixels) = combined_counts(target_pixels) >= 2;

    bf_t_mean_diff(write_mask) = pooled_iei_per_pixel(write_mask);
    counts(write_mask)         = combined_counts(write_mask);
end
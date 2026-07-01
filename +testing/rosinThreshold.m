function [auto_threshold, diagnostics] = rosinThreshold(score_map, varargin)
% ROSINTHRESHOLD  Unimodal histogram threshold via Rosin's corner method.
%
%   AUTO_THRESHOLD = ROSINTHRESHOLD(SCORE_MAP) computes a threshold
%   for a score distribution consisting of a single dominant peak
%   (noise) at the low end plus a long tail (signal), following
%   Rosin's unimodal thresholding method [1]. A straight line is
%   drawn from the histogram peak to the first empty bin following
%   the last non-empty bin; the threshold is placed at the bin of
%   maximum perpendicular distance from this line, after
%   normalizing both axes to [0, 1] for scale invariance.
%
%   This is the histogram-domain analogue of the geometric elbow
%   used in FINDELBOWTHRESHOLD, and the appropriate choice when the
%   score histogram is unimodal-plus-tail rather than bimodal
%   (where Otsu or Kittler-Illingworth would apply).
%
%   [AUTO_THRESHOLD, DIAGNOSTICS] = ROSINTHRESHOLD(SCORE_MAP) also
%   returns a struct for inspection and plotting.
%
%   [...] = ROSINTHRESHOLD(SCORE_MAP, 'Name', Value) accepts
%   optional name-value arguments:
%     'NumBins'      - (256) Number of histogram bins.
%     'SmoothWidth'  - (5) Moving-average window (in bins) applied
%                      to the histogram before corner detection.
%                      Rosin recommends smoothing because sparse
%                      tail bins are noisy [1]. Set to 1 to disable.
%     'TailQuantile' - (1.0) Anchor of the far end of the chord.
%                      1.0 reproduces Rosin's original definition
%                      (first empty bin after the last non-empty
%                      bin). Values < 1 anchor the chord at that
%                      quantile of the nonzero scores instead,
%                      which makes the chord robust to a handful
%                      of isolated extreme-score pixels (e.g. hot
%                      pixels surviving upstream filtering) that
%                      would otherwise stretch the chord and pull
%                      the detected corner toward the peak.
%
%   Inputs:
%     score_map - [H x W] 2D map of per-pixel scores (single or
%                 double). Zero and NaN entries are inactive and
%                 are excluded from the histogram.
%
%   Outputs:
%     auto_threshold - Scalar threshold at the histogram corner
%                      (bin-center value). NaN if degenerate.
%     diagnostics    - Struct with fields:
%       .bin_centers   - [1 x B] Histogram bin centers.
%       .counts        - [1 x B] Raw histogram counts.
%       .counts_smooth - [1 x B] Smoothed counts used for detection.
%       .peak_idx      - Index of the histogram peak bin.
%       .end_idx       - Index of the chord's far anchor bin.
%       .perp_dist     - [1 x B] Perpendicular distance from the
%                        chord (zero outside [peak_idx, end_idx]).
%       .thresh_idx    - Index of the detected corner bin.
%       .chord_x       - [1 x 2] Chord endpoints, x (score values).
%       .chord_y       - [1 x 2] Chord endpoints, y (counts).
%
%   Notes:
%     - The perpendicular distance is computed on axes normalized
%       to [0, 1] over the chord span, so the result is invariant
%       to rescaling of either scores or counts.
%     - If the dominant peak does not sit at the low end (i.e. the
%       distribution is not unimodal-plus-right-tail), the method's
%       assumption is violated; a warning is raised and NaN is
%       returned rather than an arbitrary threshold.
%     - Degenerate cases (fewer than 10 active pixels, or fewer
%       than 3 bins between peak and tail anchor) return NaN with
%       a warning.
%
%   References:
%     [1] P. L. Rosin, "Unimodal thresholding," Pattern
%         Recognition, vol. 34, no. 11, pp. 2083-2096, 2001.
%         DOI: 10.1016/S0031-3203(00)00136-9
%
%   See also: findElbowThreshold, histcounts

    % ----------------------------------------------------------------
    % 0. Parse arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addRequired(p, 'score_map', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'NumBins',      256, @(x) isscalar(x) && x >= 8);
    addParameter(p, 'SmoothWidth',  1,   @(x) isscalar(x) && x >= 1);
    addParameter(p, 'TailQuantile', 1.0, @(x) isscalar(x) && x > 0 && x <= 1);
    parse(p, score_map, varargin{:});

    n_bins   = p.Results.NumBins;
    smooth_w = p.Results.SmoothWidth;
    tail_q   = p.Results.TailQuantile;

    auto_threshold = NaN;
    diagnostics    = struct();

    % ----------------------------------------------------------------
    % 1. Histogram of active scores
    % ----------------------------------------------------------------
    vals = score_map(score_map > 0 & ~isnan(score_map));
    vals = vals(:);

    if numel(vals) < 10
        warning('rosinThreshold:tooFewActive', ...
            'Score map has fewer than 10 nonzero pixels. Returning NaN.');
        return;
    end

    [counts, edges] = histcounts(vals, n_bins);
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

    % ----------------------------------------------------------------
    % 2. Smooth (sparse tail bins are noisy; see Rosin [1], Sec. 4)
    % ----------------------------------------------------------------
    if smooth_w > 1
        counts_s = movmean(counts, smooth_w);
    else
        counts_s = counts;
    end

    % ----------------------------------------------------------------
    % 3. Chord anchors: peak and tail end
    % ----------------------------------------------------------------
    [peak_count, peak_idx] = max(counts_s);

    if peak_idx > round(n_bins / 2)
        warning('rosinThreshold:peakNotAtLowEnd', ...
            ['Histogram peak is in the upper half of the score range; ', ...
             'distribution is not unimodal-plus-right-tail. Returning NaN.']);
        return;
    end

    if tail_q < 1
        % Robust anchor: bin containing the requested quantile.
        % Sort-based to avoid a Statistics Toolbox dependency.
        sorted_v = sort(vals);
        q_val    = sorted_v(max(1, round(tail_q * numel(sorted_v))));
        end_idx = find(bin_centers >= q_val, 1, 'first');
        if isempty(end_idx), end_idx = n_bins; end
    else
        % Rosin's original anchor: first empty bin after the last
        % non-empty bin (y = 0 there by construction)
        last_filled = find(counts_s > 0, 1, 'last');
        end_idx     = min(last_filled + 1, n_bins);
    end

    if end_idx - peak_idx < 3
        warning('rosinThreshold:degenerateSpan', ...
            'Fewer than 3 bins between peak and tail anchor. Returning NaN.');
        return;
    end

    % ----------------------------------------------------------------
    % 4. Normalize both axes to [0, 1] over the chord span
    % ----------------------------------------------------------------
    idx_span = peak_idx:end_idx;
    x_span   = bin_centers(idx_span);
    y_span   = counts_s(idx_span);

    x_range = x_span(end) - x_span(1);
    y_range = peak_count - y_span(end);

    if x_range == 0 || y_range <= 0
        warning('rosinThreshold:degenerateRange', ...
            'Zero span on one axis between peak and tail anchor. Returning NaN.');
        return;
    end

    x_n = (x_span - x_span(1)) / x_range;
    y_n = (y_span - y_span(end)) / y_range;

    % ----------------------------------------------------------------
    % 5. Corner: max perpendicular distance from the chord
    % ----------------------------------------------------------------
    v_chord = [x_n(end) - x_n(1), y_n(end) - y_n(1)];
    v_hat   = v_chord / norm(v_chord);

    w_x = x_n - x_n(1);
    w_y = y_n - y_n(1);
    perp_span = abs(w_x * v_hat(2) - w_y * v_hat(1));

    [~, rel_idx] = max(perp_span);
    thresh_idx   = idx_span(rel_idx);

    auto_threshold = bin_centers(thresh_idx);

    % ----------------------------------------------------------------
    % 6. Pack diagnostics
    % ----------------------------------------------------------------
    if nargout > 1
        perp_full = zeros(1, n_bins);
        perp_full(idx_span) = perp_span;

        diagnostics.bin_centers   = bin_centers;
        diagnostics.counts        = counts;
        diagnostics.counts_smooth = counts_s;
        diagnostics.peak_idx      = peak_idx;
        diagnostics.end_idx       = end_idx;
        diagnostics.perp_dist     = perp_full;
        diagnostics.thresh_idx    = thresh_idx;
        diagnostics.chord_x       = [bin_centers(peak_idx), bin_centers(end_idx)];
        diagnostics.chord_y       = [counts_s(peak_idx),   counts_s(end_idx)];
    end
end

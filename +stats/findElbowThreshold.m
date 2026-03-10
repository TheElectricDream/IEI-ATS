function [auto_threshold, diagnostics] = findElbowThreshold(score_map, varargin)
% FINDELBOWTHRESHOLD  Data-driven threshold via geometric elbow detection.
%
%   AUTO_THRESHOLD = FINDELBOWTHRESHOLD(SCORE_MAP) sweeps a range of
%   candidate thresholds over SCORE_MAP, computes the mean nearest-
%   neighbour (NN) distance among surviving pixels at each threshold,
%   and returns the threshold at the geometric elbow of the resulting
%   dispersion curve. The elbow is defined as the point of maximum
%   perpendicular distance from the chord connecting the first and
%   last points of the normalized curve, following the Kneedle
%   convention [1].
%
%   [AUTO_THRESHOLD, DIAGNOSTICS] = FINDELBOWTHRESHOLD(SCORE_MAP)
%   also returns a struct of intermediate data for inspection and
%   plotting (see Outputs below).
%
%   [...] = FINDELBOWTHRESHOLD(SCORE_MAP, 'Name', Value, ...) accepts
%   optional name-value arguments described below.
%
%   Inputs:
%     score_map - [H x W] 2D map of per-pixel scores (single or
%                 double). Any map that produces a scalar per pixel
%                 can be used: coherence map, similarity map, trace
%                 map, etc. Zero and NaN entries are treated as
%                 inactive pixels and excluded from the analysis.
%
%   Name-Value Arguments:
%     'NumThresholds'  - (50)   Number of uniformly spaced candidate
%                        thresholds to evaluate between the lower and
%                        upper percentile bounds.
%     'PrctileLo'      - (2)    Lower percentile of nonzero values
%                        for the sweep range. Avoids outlier-driven
%                        range inflation at the low end.
%     'PrctileHi'      - (98)   Upper percentile for the sweep range.
%     'MinPoints'       - (5)    Minimum number of surviving pixels
%                        required to compute a valid NN distance.
%                        Thresholds that yield fewer points are
%                        marked NaN and excluded from elbow detection.
%     'MaxSampleSize'   - (8000) Maximum number of pixels to use for
%                        the KNN computation at each threshold level.
%                        If more pixels survive, a uniform random
%                        subsample of this size is drawn. This bounds
%                        the cost of the KD-tree build to O(N log N)
%                        with N <= MaxSampleSize.
%     'KNeighbours'     - (1)    Number of nearest neighbours whose
%                        mean distance forms the dispersion metric.
%                        K=1 is standard mean NN distance. K>1
%                        averages over multiple neighbours, which is
%                        more robust to isolated outliers but smooths
%                        the curve.
%
%   Outputs:
%     auto_threshold - Scalar threshold at the geometric elbow.
%     diagnostics    - Struct with fields:
%       .th_vec          - [1 x M] Evaluated thresholds (valid only).
%       .mean_nn         - [1 x M] Mean NN distance at each threshold.
%       .n_points        - [1 x M] Number of surviving pixels.
%       .th_norm         - [1 x M] Thresholds normalized to [0, 1].
%       .nn_norm         - [1 x M] NN distances normalized to [0, 1].
%       .perp_dist       - [1 x M] Perpendicular distance from chord.
%       .elbow_idx       - Index into the valid vectors.
%       .th_vec_full     - [1 x N_th] Full sweep vector (incl. NaN).
%       .mean_nn_full    - [1 x N_th] Full NN vector (incl. NaN).
%
%   Algorithm:
%     1. Extract nonzero values from the score map and build a
%        percentile-bounded sweep vector of candidate thresholds.
%     2. For each candidate threshold, find all pixels >= threshold,
%        extract their (row, col) coordinates, subsample if needed,
%        and compute the mean distance to the K nearest neighbours
%        via knnsearch (MATLAB's Statistics Toolbox KD-tree).
%     3. Prune thresholds where fewer than MinPoints survived.
%     4. Normalize both axes to [0, 1] so that the perpendicular
%        distance is scale-invariant (critical — see [1] §3).
%     5. Compute the perpendicular distance of each point from the
%        chord connecting the first and last normalized points.
%     6. The elbow is the point of maximum perpendicular distance.
%
%   Notes:
%     - The expected curve shape is: low NN distance at low thresholds
%       (dense noise coverage), rising steeply once signal-only
%       clusters remain. The elbow marks the transition from
%       "removing noise" to "removing signal."
%     - If the curve is non-monotonic (U-shaped), the elbow detector
%       still returns the point of maximum chord deviation. Consider
%       also inspecting the minimum of the raw curve as an
%       alternative operating point.
%     - Computational cost: O(N_th * N_sub * log(N_sub)) where
%       N_sub = min(n_surviving, MaxSampleSize).
%
%   References:
%     [1] V. Satopaa, J. Albrecht, D. Irwin, and B. Raghavan,
%         "Finding a 'Kneedle' in a Haystack: Detecting Knee Points
%         in System Behavior," Proc. 31st Int. Conf. Distributed
%         Computing Systems Workshops, pp. 166-171, 2011.
%         DOI: 10.1109/ICDCSW.2011.20
%     [2] P. C. Hansen, "Analysis of Discrete Ill-Posed Problems by
%         Means of the L-Curve," SIAM Review, vol. 34, no. 4,
%         pp. 561-580, 1992.
%
%   See also: knnsearch, prctile, coherence.computeCoherenceMask

    % ----------------------------------------------------------------
    % 0. Parse optional arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'NumThresholds',  50,   @(x) isscalar(x) && x >= 3);
    addParameter(p, 'PrctileLo',      2,    @(x) isscalar(x) && x >= 0);
    addParameter(p, 'PrctileHi',      98,   @(x) isscalar(x) && x <= 100);
    addParameter(p, 'MinPoints',      5,    @(x) isscalar(x) && x >= 2);
    addParameter(p, 'MaxSampleSize',  8000, @(x) isscalar(x) && x >= 10);
    addParameter(p, 'KNeighbours',    1,    @(x) isscalar(x) && x >= 1);
    parse(p, varargin{:});
    opts = p.Results;

    N_th       = opts.NumThresholds;
    pct_lo     = opts.PrctileLo;
    pct_hi     = opts.PrctileHi;
    min_pts    = opts.MinPoints;
    max_sample = opts.MaxSampleSize;
    K_nn       = opts.KNeighbours;

    % ----------------------------------------------------------------
    % 1. Build sweep vector from nonzero map values
    % ----------------------------------------------------------------
    vals = score_map(score_map > 0 & ~isnan(score_map));

    if numel(vals) < min_pts
        warning('findElbowThreshold:tooFewActive', ...
            'Score map has fewer than %d nonzero pixels. Returning NaN.', ...
            min_pts);
        auto_threshold = NaN;
        diagnostics = struct();
        return;
    end

    lo = prctile(vals, pct_lo);
    hi = prctile(vals, pct_hi);

    if lo >= hi
        % Degenerate case: nearly all values identical
        auto_threshold = lo;
        diagnostics = struct();
        return;
    end

    th_vec_full  = linspace(lo, hi, N_th);
    mean_nn_full = nan(1, N_th);
    n_pts_full   = nan(1, N_th);

    % ----------------------------------------------------------------
    % 2. Sweep: compute mean NN distance at each threshold
    % ----------------------------------------------------------------
    % K+1 because knnsearch returns self as the first neighbour
    K_query = K_nn + 1;

    for k = 1:N_th
        mask = score_map >= th_vec_full(k);
        [r, c] = find(mask);
        n = numel(r);
        n_pts_full(k) = n;

        if n < min_pts
            continue;
        end

        % Subsample if needed to bound KD-tree cost
        if n > max_sample
            idx = randperm(n, max_sample);
            r = r(idx);
            c = c(idx);
        end

        [~, D] = knnsearch([r c], [r c], 'K', K_query);
        % Mean over neighbours 2:end (skip self at column 1)
        mean_nn_full(k) = mean(mean(D(:, 2:end), 2));
    end

    % ----------------------------------------------------------------
    % 3. Prune invalid entries
    % ----------------------------------------------------------------
    good = ~isnan(mean_nn_full);

    if sum(good) < 3
        warning('findElbowThreshold:tooFewValid', ...
            'Fewer than 3 valid threshold evaluations. Returning median.');
        auto_threshold = median(th_vec_full);
        diagnostics = struct();
        return;
    end

    th_g = th_vec_full(good);
    nn_g = mean_nn_full(good);
    np_g = n_pts_full(good);

    % ----------------------------------------------------------------
    % 4. Normalize both axes to [0, 1]
    % ----------------------------------------------------------------
    th_n = (th_g - th_g(1)) / (th_g(end) - th_g(1));

    nn_range = max(nn_g) - min(nn_g);
    if nn_range == 0
        % Flat curve — no meaningful elbow
        auto_threshold = th_g(round(end/2));
        diagnostics = struct();
        return;
    end
    nn_n = (nn_g - min(nn_g)) / nn_range;

    % ----------------------------------------------------------------
    % 5. Geometric elbow: max perpendicular distance from chord
    % ----------------------------------------------------------------
    p1 = [th_n(1),   nn_n(1)];
    p2 = [th_n(end), nn_n(end)];
    v  = p2 - p1;
    v_hat = v / norm(v);

    perp_dist = zeros(size(th_n));
    for k = 1:numel(th_n)
        w = [th_n(k), nn_n(k)] - p1;
        % Signed perpendicular distance (magnitude only)
        perp_dist(k) = abs(w(1)*v_hat(2) - w(2)*v_hat(1));
    end

    [~, elbow_idx] = max(perp_dist);
    auto_threshold = th_g(elbow_idx);

    % ----------------------------------------------------------------
    % 6. Pack diagnostics
    % ----------------------------------------------------------------
    if nargout > 1
        diagnostics.th_vec       = th_g;
        diagnostics.mean_nn      = nn_g;
        diagnostics.n_points     = np_g;
        diagnostics.th_norm      = th_n;
        diagnostics.nn_norm      = nn_n;
        diagnostics.perp_dist    = perp_dist;
        diagnostics.elbow_idx    = elbow_idx;
        diagnostics.th_vec_full  = th_vec_full;
        diagnostics.mean_nn_full = mean_nn_full;
    end

end
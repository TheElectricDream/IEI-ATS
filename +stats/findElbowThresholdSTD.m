function [auto_threshold, diagnostics] = findElbowThresholdSTD(score_map, varargin)
% FINDELBOWTHRESHOLD  Data-driven threshold via geometric elbow detection.
%
%   AUTO_THRESHOLD = FINDELBOWTHRESHOLD(SCORE_MAP) sweeps a range of
%   candidate thresholds over SCORE_MAP, computes the spatial spread
%   (standard deviation) of surviving pixels at each threshold, and
%   returns the threshold at the geometric elbow of the resulting
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
%     'NumThresholds'  - (50)   Number of uniformly spaced sample
%                        points along the sorted pixel array at which
%                        the dispersion curve is evaluated.
%     'PrctileLo'      - (2)    Lower percentile of nonzero values.
%                        Pixels below this are excluded from the
%                        sweep to avoid outlier-driven range
%                        inflation at the low end.
%     'PrctileHi'      - (98)   Upper percentile for the sweep range.
%     'MinPoints'       - (5)    Minimum number of surviving pixels
%                        required for a valid spread estimate.
%
%   Outputs:
%     auto_threshold - Scalar threshold at the geometric elbow.
%     diagnostics    - Struct with fields:
%       .th_vec          - [1 x M] Threshold values at sample points.
%       .spread          - [1 x M] Spatial spread (px) at each level.
%       .n_points        - [1 x M] Number of surviving pixels.
%       .th_norm         - [1 x M] Thresholds normalized to [0, 1].
%       .spread_norm     - [1 x M] Spread normalized to [0, 1].
%       .perp_dist       - [1 x M] Perpendicular distance from chord.
%       .elbow_idx       - Index into the sample vectors.
%
%   Algorithm:
%     1. Extract nonzero pixel coordinates and scores.
%     2. Sort descending by score.
%     3. Compute cumulative sums of (r, c, r^2, c^2) along the
%        sorted order. At position k the "surviving cloud" is the
%        set of pixels with score >= sorted_score(k), i.e., the
%        first k entries. The spatial variance at k is computed in
%        O(1) from the cumulative sums:
%          var_r(k) = cs_r2(k)/k - (cs_r(k)/k)^2
%        and the total spread is sqrt(var_r + var_c).
%     4. Subsample the curve at NumThresholds evenly spaced indices.
%     5. Normalize both axes to [0, 1] and find the geometric elbow
%        (max perpendicular distance from the chord).
%
%   Computational cost:
%     O(N log N) for the sort, O(N) for cumulative sums and
%     variance computation. No KD-tree construction. Typical
%     execution on a 640 x 480 coherence map: < 5 ms.
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
%   See also: coherence.computeCoherenceMask, plot.plotElbowDiagnostics

    % ----------------------------------------------------------------
    % 0. Parse optional arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'NumThresholds', 50,  @(x) isscalar(x) && x >= 3);
    addParameter(p, 'PrctileLo',     2,   @(x) isscalar(x) && x >= 0);
    addParameter(p, 'PrctileHi',     98,  @(x) isscalar(x) && x <= 100);
    addParameter(p, 'MinPoints',     5,   @(x) isscalar(x) && x >= 2);
    parse(p, varargin{:});
    opts = p.Results;

    N_th    = opts.NumThresholds;
    pct_lo  = opts.PrctileLo;
    pct_hi  = opts.PrctileHi;
    min_pts = opts.MinPoints;

    % ----------------------------------------------------------------
    % 1. Extract active pixels
    % ----------------------------------------------------------------
    valid = score_map > 0 & ~isnan(score_map);
    [r, c] = find(valid);
    v = score_map(valid);

    if numel(v) < min_pts
        warning('findElbowThreshold:tooFewActive', ...
            'Score map has fewer than %d nonzero pixels. Returning NaN.', ...
            min_pts);
        auto_threshold = NaN;
        diagnostics    = struct();
        return;
    end

    % ----------------------------------------------------------------
    % 2. Percentile clipping and sort descending by score
    % ----------------------------------------------------------------
    lo = prctile(v, pct_lo);
    hi = prctile(v, pct_hi);

    keep = v >= lo & v <= hi;
    r = r(keep);
    c = c(keep);
    v = v(keep);

    [v, order] = sort(v, 'descend');
    r = double(r(order));
    c = double(c(order));
    N = numel(v);

    if N < min_pts
        auto_threshold = lo;
        diagnostics    = struct();
        return;
    end

    % ----------------------------------------------------------------
    % 3. Cumulative sums for incremental variance
    % ----------------------------------------------------------------
    %   As we sweep from high -> low threshold (index 1 -> N), the
    %   cloud grows. At index k the cloud contains the k highest-
    %   scored pixels. Variance in O(1):
    %     var_r(k) = cs_r2(k)/k - (cs_r(k)/k)^2
    %
    cs_r  = cumsum(r);
    cs_c  = cumsum(c);
    cs_r2 = cumsum(r .* r);
    cs_c2 = cumsum(c .* c);
    n_vec = (1:N)';

    var_r = cs_r2 ./ n_vec - (cs_r ./ n_vec).^2;
    var_c = cs_c2 ./ n_vec - (cs_c ./ n_vec).^2;

    % Guard against floating-point negative variance at very small k
    var_r = max(var_r, 0);
    var_c = max(var_c, 0);

    spread_full = sqrt(var_r + var_c);

    % ----------------------------------------------------------------
    % 4. Subsample at N_th evenly spaced indices
    % ----------------------------------------------------------------
    sample_idx = unique(round(linspace(max(min_pts, 1), N, N_th)));
    M = numel(sample_idx);

    th_vec   = v(sample_idx)';
    spread_s = spread_full(sample_idx)';
    n_pts    = n_vec(sample_idx)';

    % Remove degenerate samples (zero spread or too few points)
    good = spread_s > 0 & n_pts >= min_pts;
    if sum(good) < 3
        auto_threshold = median(th_vec);
        diagnostics    = struct();
        return;
    end

    th_vec   = th_vec(good);
    spread_s = spread_s(good);
    n_pts    = n_pts(good);
    M        = numel(th_vec);

    % Flip to ascending threshold order for the elbow geometry
    % (low threshold on left  -> many points -> high spread,
    %  high threshold on right -> few points  -> low spread)
    th_vec   = fliplr(th_vec);
    spread_s = fliplr(spread_s);
    n_pts    = fliplr(n_pts);

    % ----------------------------------------------------------------
    % 5. Normalize both axes to [0, 1]
    % ----------------------------------------------------------------
    th_range = th_vec(end) - th_vec(1);
    sp_range = max(spread_s) - min(spread_s);

    if th_range == 0 || sp_range == 0
        auto_threshold = th_vec(round(M/2));
        diagnostics    = struct();
        return;
    end

    th_n = (th_vec - th_vec(1))       / th_range;
    sp_n = (spread_s - min(spread_s)) / sp_range;

    % ----------------------------------------------------------------
    % 6. Geometric elbow: max perpendicular distance from chord
    % ----------------------------------------------------------------
    v_chord = [th_n(end) - th_n(1), sp_n(end) - sp_n(1)];
    v_hat   = v_chord / norm(v_chord);

    % Vectorized cross-product magnitude (2D)
    w_x = th_n - th_n(1);
    w_y = sp_n - sp_n(1);
    perp_dist = abs(w_x * v_hat(2) - w_y * v_hat(1));

    [~, elbow_idx] = max(perp_dist);
    auto_threshold = th_vec(elbow_idx);

    % ----------------------------------------------------------------
    % 7. Pack diagnostics
    % ----------------------------------------------------------------
    if nargout > 1
        diagnostics.th_vec      = th_vec;
        diagnostics.spread      = spread_s;
        diagnostics.n_points    = n_pts;
        diagnostics.th_norm     = th_n;
        diagnostics.spread_norm = sp_n;
        diagnostics.perp_dist   = perp_dist;
        diagnostics.elbow_idx   = elbow_idx;
    end

end
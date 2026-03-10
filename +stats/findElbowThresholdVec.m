function [auto_threshold, diagnostics] = findElbowThresholdVec(score_map, varargin)
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
%   The implementation builds a single KD-tree over all active pixels
%   and precomputes K nearest global neighbours. The threshold sweep
%   then uses vectorized logical indexing against the sorted rank
%   array to determine which precomputed neighbours are active at
%   each threshold level, avoiding per-threshold tree construction.
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
%     'NumThresholds'   - (50)  Number of uniformly spaced sample
%                         points along the sorted pixel array at
%                         which the dispersion curve is evaluated.
%     'PrctileLo'       - (2)   Lower percentile of nonzero values.
%                         Pixels below this are excluded from the
%                         sweep to avoid outlier-driven range
%                         inflation.
%     'PrctileHi'       - (98)  Upper percentile for the sweep range.
%     'MinPoints'        - (5)   Minimum number of surviving pixels
%                         required for a valid NN estimate.
%     'KPrecompute'      - (10)  Number of nearest neighbours to
%                         precompute per pixel. Higher values handle
%                         very sparse active sets at high thresholds
%                         but increase memory. 10 is sufficient for
%                         most event-camera score maps.
%
%   Outputs:
%     auto_threshold - Scalar threshold at the geometric elbow.
%     diagnostics    - Struct with fields:
%       .th_vec       - [1 x M] Threshold values at sample points
%                       (ascending order).
%       .mean_nn      - [1 x M] Mean NN distance [px] at each level.
%       .n_points     - [1 x M] Number of pixels with a valid NN.
%       .th_norm      - [1 x M] Thresholds normalized to [0, 1].
%       .nn_norm      - [1 x M] NN distances normalized to [0, 1].
%       .perp_dist    - [1 x M] Perpendicular distance from chord.
%       .elbow_idx    - Index into the sample vectors.
%
%   Algorithm:
%     1. Extract nonzero pixel coordinates and scores.
%     2. Sort descending by score and record each pixel's rank.
%     3. Build one KD-tree on all N active pixels. Query each
%        pixel's K nearest neighbours (excluding self). Store
%        neighbour distances D (N x K, sorted ascending by
%        distance) and neighbour ranks R (N x K).
%     4. Subsample at NumThresholds evenly spaced indices along
%        the sorted order. At sample index k the active set is
%        pixels with rank <= k (the k highest-scored pixels).
%     5. For each sample k, the NN distance of pixel i (rank <= k)
%        to its nearest ACTIVE neighbour is:
%          min over j where R(i,j) <= k of D(i,j)
%        Since D is sorted ascending, masking inactive neighbours
%        to Inf and taking the row-wise minimum is O(k * K).
%     6. Normalize both axes to [0, 1] and find the geometric
%        elbow (max perpendicular distance from chord).
%
%   Computational cost:
%     O(N log N) for the sort + KD-tree build + K-NN query (once).
%     O(sum_k k * K) for the sweep, which is approximately
%     O(N * K * N_th / 2). With defaults (N=30k, K=10, N_th=50):
%     ~7.5 M vectorized operations, typically < 20 ms in MATLAB.
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
%   See also: knnsearch, coherence.computeCoherenceMask,
%             plot.plotElbowDiagnostics

    % ----------------------------------------------------------------
    % 0. Parse optional arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'NumThresholds',  50,  @(x) isscalar(x) && x >= 3);
    addParameter(p, 'PrctileLo',      2,   @(x) isscalar(x) && x >= 0);
    addParameter(p, 'PrctileHi',      98,  @(x) isscalar(x) && x <= 100);
    addParameter(p, 'MinPoints',      5,   @(x) isscalar(x) && x >= 2);
    addParameter(p, 'KPrecompute',    10,  @(x) isscalar(x) && x >= 1);
    parse(p, varargin{:});
    opts = p.Results;

    N_th    = opts.NumThresholds;
    pct_lo  = opts.PrctileLo;
    pct_hi  = opts.PrctileHi;
    min_pts = opts.MinPoints;
    K_pre   = opts.KPrecompute;

    % ----------------------------------------------------------------
    % 1. Extract active pixels
    % ----------------------------------------------------------------
    valid_mask = score_map > 0 & ~isnan(score_map);
    [r, c] = find(valid_mask);
    v = score_map(valid_mask);

    if numel(v) < min_pts
        warning('findElbowThreshold:tooFewActive', ...
            'Score map has fewer than %d nonzero pixels. Returning NaN.', ...
            min_pts);
        auto_threshold = NaN;
        diagnostics    = struct();
        return;
    end

    % ----------------------------------------------------------------
    % 2. Percentile clipping + sort descending by score
    % ----------------------------------------------------------------
    lo = prctile(v, pct_lo);
    hi = prctile(v, pct_hi);

    keep = v >= lo & v <= hi;
    r = double(r(keep));
    c = double(c(keep));
    v = double(v(keep));

    [v, sort_order] = sort(v, 'descend');
    r = r(sort_order);
    c = c(sort_order);
    N = numel(v);

    if N < min_pts
        auto_threshold = lo;
        diagnostics    = struct();
        return;
    end

    % Clamp K_pre to N-1 (can't have more neighbours than points)
    K_pre = min(K_pre, N - 1);

    % ----------------------------------------------------------------
    % 3. Build KD-tree ONCE; precompute K nearest neighbours
    % ----------------------------------------------------------------
    %   knnsearch returns K+1 results with self as the first match.
    %   We strip self to get the K true nearest neighbours.
    %
    %   D: N x K_pre — distances, sorted ascending per row.
    %   IDX: N x K_pre — indices into the sorted (r,c,v) arrays.
    %
    coords = [r, c];
    [IDX, D] = knnsearch(coords, coords, 'K', K_pre + 1);
    IDX = IDX(:, 2:end);   % strip self-match
    D   = D(:, 2:end);

    % ----------------------------------------------------------------
    % 4. Compute neighbour ranks in the sorted order
    % ----------------------------------------------------------------
    %   Since r, c, v are already in sort_order, the "rank" of pixel
    %   at position i is simply i. The rank of its j-th precomputed
    %   neighbour is IDX(i,j) — because IDX indexes into the sorted
    %   arrays directly.
    %
    %   At threshold sample k, a neighbour at rank IDX(i,j) is
    %   "active" iff IDX(i,j) <= k.
    %
    nbr_ranks = IDX;   % N x K_pre, already in sorted-array indices

    % ----------------------------------------------------------------
    % 5. Sweep: vectorized NN distance at each threshold sample
    % ----------------------------------------------------------------
    sample_idx = unique(round(linspace(max(min_pts, 1), N, N_th)));
    M = numel(sample_idx);

    mean_nn_vec = nan(1, M);
    n_valid_vec = nan(1, M);
    th_vec      = v(sample_idx)';

    for s = 1:M
        k = sample_idx(s);

        % Logical mask: which precomputed neighbours of pixels 1:k
        % are themselves in the active set (rank <= k)?
        active_nbr = nbr_ranks(1:k, :) <= k;   % k x K_pre logical

        % Copy distances and mask out inactive neighbours
        D_k = D(1:k, :);                        % k x K_pre
        D_k(~active_nbr) = Inf;

        % Row-wise minimum gives NN distance to nearest active neighbour
        nn_dists = min(D_k, [], 2);              % k x 1

        % Pixels with no active precomputed neighbour get Inf; exclude
        has_nn = isfinite(nn_dists);
        n_valid_vec(s) = sum(has_nn);

        if n_valid_vec(s) >= min_pts
            mean_nn_vec(s) = mean(nn_dists(has_nn));
        end
    end

    % ----------------------------------------------------------------
    % 6. Prune invalid samples
    % ----------------------------------------------------------------
    good = ~isnan(mean_nn_vec) & n_valid_vec >= min_pts;

    if sum(good) < 3
        warning('findElbowThreshold:tooFewValid', ...
            'Fewer than 3 valid threshold evaluations. Returning median.');
        auto_threshold = median(th_vec);
        diagnostics    = struct();
        return;
    end

    th_g  = th_vec(good);
    nn_g  = mean_nn_vec(good);
    np_g  = n_valid_vec(good);

    % Flip to ascending threshold order for elbow geometry
    % (low threshold → many points → low NN on the left,
    %  high threshold → few points → high NN on the right)
    th_g = fliplr(th_g);
    nn_g = fliplr(nn_g);
    np_g = fliplr(np_g);

    % ----------------------------------------------------------------
    % 7. Normalize both axes to [0, 1]
    % ----------------------------------------------------------------
    th_range = th_g(end) - th_g(1);
    nn_range = max(nn_g) - min(nn_g);

    if th_range == 0 || nn_range == 0
        auto_threshold = th_g(round(numel(th_g) / 2));
        diagnostics    = struct();
        return;
    end

    th_n = (th_g - th_g(1))       / th_range;
    nn_n = (nn_g - min(nn_g))     / nn_range;

    % ----------------------------------------------------------------
    % 8. Geometric elbow: max perpendicular distance from chord
    % ----------------------------------------------------------------
    v_chord = [th_n(end) - th_n(1), nn_n(end) - nn_n(1)];
    v_hat   = v_chord / norm(v_chord);

    % Vectorized 2D cross-product magnitude
    w_x = th_n - th_n(1);
    w_y = nn_n - nn_n(1);
    perp_dist = abs(w_x * v_hat(2) - w_y * v_hat(1));

    [~, elbow_idx] = max(perp_dist);
    auto_threshold = th_g(elbow_idx);

    % ----------------------------------------------------------------
    % 9. Pack diagnostics
    % ----------------------------------------------------------------
    if nargout > 1
        diagnostics.th_vec    = th_g;
        diagnostics.mean_nn   = nn_g;
        diagnostics.n_points  = np_g;
        diagnostics.th_norm   = th_n;
        diagnostics.nn_norm   = nn_n;
        diagnostics.perp_dist = perp_dist;
        diagnostics.elbow_idx = elbow_idx;
    end

end
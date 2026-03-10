function [auto_threshold, diagnostics] = findElbowThresholdCheap(score_map, varargin)
% FINDELBOWTHRESHOLD  Data-driven threshold via geometric elbow detection.
%
%   AUTO_THRESHOLD = FINDELBOWTHRESHOLD(SCORE_MAP) sweeps candidate
%   thresholds over SCORE_MAP, computes the mean local binary
%   variance of the surviving pixel mask at each threshold, and
%   returns the threshold at the geometric elbow of the resulting
%   curve. The elbow is defined as the point of maximum perpendicular
%   distance from the chord connecting the curve endpoints after
%   normalizing both axes to [0, 1], following the Kneedle
%   convention [1].
%
%   The local binary variance at each pixel is p(1-p), where p is
%   the fraction of active pixels within a local window. This is
%   computed via a single conv2 call per threshold — no KD-trees or
%   nearest-neighbour searches. It captures local spatial clustering:
%   mixed noise-and-signal regions produce high p(1-p) while
%   homogeneous regions (pure cluster or pure background) produce
%   values near zero. The mean across the image summarizes the
%   overall degree of spatial mixing at each threshold level.
%
%   [AUTO_THRESHOLD, DIAGNOSTICS] = FINDELBOWTHRESHOLD(SCORE_MAP)
%   also returns a struct for inspection and plotting.
%
%   [...] = FINDELBOWTHRESHOLD(SCORE_MAP, N_TH) specifies the number
%   of candidate thresholds (default: 50).
%
%   [...] = FINDELBOWTHRESHOLD(SCORE_MAP, N_TH, WIN_SZ) specifies
%   the local window half-size in pixels (default: 5, giving an
%   11 x 11 window). Larger windows measure clustering at coarser
%   spatial scales.
%
%   Inputs:
%     score_map - [H x W] 2D map of per-pixel scores (single or
%                 double). Zero and NaN entries are inactive.
%     N_TH      - (Optional, default 50) Number of thresholds.
%     WIN_SZ    - (Optional, default 5) Window half-size [pixels].
%
%   Outputs:
%     auto_threshold - Scalar threshold at the geometric elbow.
%     diagnostics    - Struct with fields:
%       .th_vec       - [1 x M] Candidate thresholds (ascending).
%       .mean_lbv     - [1 x M] Mean local binary variance.
%       .th_norm      - [1 x M] Thresholds normalized to [0, 1].
%       .lbv_norm     - [1 x M] LBV normalized to [0, 1].
%       .perp_dist    - [1 x M] Perpendicular distance from chord.
%       .elbow_idx    - Index into the vectors above.
%
%   Algorithm:
%     1. Extract unique nonzero score values and build a linearly
%        spaced sweep vector from min to max.
%     2. At each candidate threshold, binarize the map (>= th).
%     3. Compute local density p via box-filter convolution:
%          p = conv2(mask, kernel, 'same') / numel(kernel)
%     4. Mean local binary variance = mean(p .* (1 - p)) over
%        the full image. This is high when the mask is spatially
%        mixed (noise + signal) and low when it is homogeneous.
%     5. Normalize both axes to [0, 1] and find the elbow
%        (max perpendicular distance from chord).
%
%   Computational cost:
%     O(N_TH x H x W) via conv2. For 640 x 480 with N_TH = 50:
%     ~50 convolutions, typically < 10 ms total in MATLAB.
%
%   References:
%     [1] V. Satopaa, J. Albrecht, D. Irwin, and B. Raghavan,
%         "Finding a 'Kneedle' in a Haystack: Detecting Knee Points
%         in System Behavior," Proc. 31st Int. Conf. Distributed
%         Computing Systems Workshops, pp. 166-171, 2011.
%         DOI: 10.1109/ICDCSW.2011.20
%
%   See also: conv2, coherence.computeCoherenceMask,
%             plot.plotElbowDiagnostics

    % ----------------------------------------------------------------
    % 0. Parse arguments (positional for simplicity)
    % ----------------------------------------------------------------
    N_th   = 50;
    win_hz = 5;

    if nargin >= 2 && ~isempty(varargin{1})
        N_th = varargin{1};
    end
    if nargin >= 3 && ~isempty(varargin{2})
        win_hz = varargin{2};
    end

    % ----------------------------------------------------------------
    % 1. Build threshold sweep from actual data range
    % ----------------------------------------------------------------
    vals = score_map(score_map > 0 & ~isnan(score_map));

    if numel(vals) < 10
        warning('findElbowThreshold:tooFewActive', ...
            'Score map has fewer than 10 nonzero pixels. Returning NaN.');
        auto_threshold = NaN;
        diagnostics    = struct();
        return;
    end

    th_vec = linspace(min(vals), max(vals), N_th);

    % ----------------------------------------------------------------
    % 2. Box-filter kernel (uniform weights)
    % ----------------------------------------------------------------
    win_side = 2 * win_hz + 1;
    kernel   = ones(win_side) / win_side^2;

    % ----------------------------------------------------------------
    % 3. Sweep: mean local binary variance at each threshold
    % ----------------------------------------------------------------
    mean_lbv = zeros(1, N_th);

    for k = 1:N_th
        mask = double(score_map >= th_vec(k));
        p    = conv2(mask, kernel, 'same');
        lbv  = p .* (1 - p);
        mean_lbv(k) = mean(lbv(:));
    end

    % ----------------------------------------------------------------
    % 4. Prune flat leading/trailing regions
    % ----------------------------------------------------------------
    good = mean_lbv > 0;
    if sum(good) < 3
        auto_threshold = th_vec(round(N_th / 2));
        diagnostics    = struct();
        return;
    end

    th_g  = th_vec(good);
    lbv_g = mean_lbv(good);

    % ----------------------------------------------------------------
    % 5. Normalize both axes to [0, 1]
    % ----------------------------------------------------------------
    th_range  = th_g(end) - th_g(1);
    lbv_range = max(lbv_g) - min(lbv_g);

    if th_range == 0 || lbv_range == 0
        auto_threshold = th_g(round(numel(th_g) / 2));
        diagnostics    = struct();
        return;
    end

    th_n  = (th_g  - th_g(1))      / th_range;
    lbv_n = (lbv_g - min(lbv_g))   / lbv_range;

    % ----------------------------------------------------------------
    % 6. Geometric elbow: max perpendicular distance from chord
    % ----------------------------------------------------------------
    v_chord = [th_n(end) - th_n(1), lbv_n(end) - lbv_n(1)];
    v_hat   = v_chord / norm(v_chord);

    w_x = th_n - th_n(1);
    w_y = lbv_n - lbv_n(1);
    perp_dist = abs(w_x * v_hat(2) - w_y * v_hat(1));

    [~, elbow_idx] = max(perp_dist);
    auto_threshold = th_g(elbow_idx);

    % ----------------------------------------------------------------
    % 7. Pack diagnostics
    % ----------------------------------------------------------------
    if nargout > 1
        diagnostics.th_vec    = th_g;
        diagnostics.mean_lbv  = lbv_g;
        diagnostics.th_norm   = th_n;
        diagnostics.lbv_norm  = lbv_n;
        diagnostics.perp_dist = perp_dist;
        diagnostics.elbow_idx = elbow_idx;
    end

end
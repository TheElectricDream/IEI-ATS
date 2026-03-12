function m = computeFrameMetrics(filter_mask, sorted_x, sorted_y, ...
    sorted_t, p_signed, imgSz, counts, ...
    normalized_output_frame, prev_output_frame, ...
    n_passed, n_total)
% COMPUTEFRAMEMETRICS  Per-frame quantitative evaluation of filter + accumulator.
%
%   M = COMPUTEFRAMEMETRICS(FILTER_MASK, SORTED_X, SORTED_Y,
%   SORTED_T, P_SIGNED, IMGSZ, COUNTS, NORMALIZED_OUTPUT_FRAME,
%   PREV_OUTPUT_FRAME, N_PASSED, N_TOTAL) computes a comprehensive
%   set of quality metrics for a single processed frame, evaluating
%   both the upstream filter's selectivity and the downstream
%   accumulator's output quality.
%
%   This function is designed to be called once per frame inside
%   main.m, after both the filtering and accumulation stages have
%   completed. It is agnostic to the specific filter or accumulator
%   used — all methods produce the required inputs through the
%   existing pipeline interface.
%
%   Inputs:
%     filter_mask             - [imgSz] Binary or continuous filter
%                               mask from the selected filter.
%                               Values > 0 indicate retained pixels.
%     sorted_x, sorted_y     - [N x 1] Pixel coordinates (row, col),
%                               sorted by linear pixel index.
%     sorted_t               - [N x 1] Timestamps [s], sorted to
%                               match sorted_x/y.
%     p_signed               - [N x 1] Signed polarity (+1 or -1)
%                               for each event.
%     imgSz                  - [1 x 2] Image dimensions [nRows, nCols].
%     counts                 - [imgSz] Per-pixel event count for the
%                               current frame.
%     normalized_output_frame - [imgSz] Accumulated surface output
%                               from the selected accumulator,
%                               mapped to [0, 1].
%     prev_output_frame      - [imgSz] Previous frame's accumulated
%                               surface. Pass zeros(imgSz) or [] on
%                               the first frame.
%     n_passed               - Scalar. Number of events that passed
%                               the filter. For 'COHERENCE', compute
%                               from the mask and event locations.
%     n_total                - Scalar. Total events in this frame.
%
%   Outputs:
%     m - Scalar struct with fields:
%
%       --- Filter selectivity metrics ---
%       .srr                 - Signal Retention Rate: n_passed / n_total.
%                              Fraction of events surviving the filter.
%       .pixel_retention     - Fraction of active pixels retained:
%                              (# pixels with mask > 0 AND events) /
%                              (# pixels with events).
%       .polarity_balance    - Polarity balance of retained events:
%                              |mean(p_retained)|. Values near 0
%                              indicate balanced +/- (specular or
%                              oscillating regions). Values near 1
%                              indicate consistent polarity (clean
%                              leading/trailing edge). Key metric
%                              for IEI-ATS vs. STCC comparison.
%       .n_passed            - Absolute count of passed events.
%       .n_total             - Absolute count of total events.
%       .n_active_pixels     - Number of pixels with events.
%       .n_retained_pixels   - Number of pixels retained by filter.
%
%       --- Surface quality metrics ---
%       .edge_sharpness      - Mean gradient magnitude along detected
%                              edges in the accumulated surface.
%                              Higher values indicate sharper edges.
%                              Uses Sobel gradient via imgradient().
%       .mean_gradient       - Mean gradient magnitude over the full
%                              image. Captures overall spatial detail.
%       .contrast_iqr        - Interquartile range of nonzero pixel
%                              values in the surface. Measures
%                              dynamic range utilization. Stable
%                              values across frames indicate good
%                              normalization.
%       .noise_floor         - Mean |surface value - 0.5| in the
%                              background region (pixels with zero
%                              event count). For unsigned surfaces
%                              (HOTS, AGD, Zhu, METS), background
%                              should decay to 0. For signed
%                              surfaces (IEI-ATS), background
%                              should be at 0.5 (mid-gray).
%                              Lower values indicate better noise
%                              suppression.
%       .surface_fill        - Fraction of pixels with |surface -
%                              midpoint| > 0.01 (non-trivial value).
%                              Captures how much of the frame is
%                              occupied by signal.
%
%       --- Temporal coherence metrics ---
%       .temporal_ssim       - Structural similarity index between
%                              the current and previous frame.
%                              Values near 1 indicate smooth temporal
%                              transitions. NaN on the first frame.
%       .temporal_mae        - Mean absolute error between current
%                              and previous frame. Captures raw
%                              temporal difference without the
%                              structural weighting of SSIM.
%
%       --- Spatial coherence metrics ---
%       .n_clusters          - Number of connected components in
%                              the binarized filter mask. Fewer,
%                              larger clusters indicate better
%                              spatial coherence. Many small clusters
%                              suggest noise leakage.
%       .mean_cluster_size   - Mean size (in pixels) of connected
%                              components. Larger clusters indicate
%                              spatially coherent signal retention.
%       .largest_cluster_frac - Fraction of retained pixels in the
%                              largest connected component. High
%                              values indicate that the filter is
%                              selecting a single dominant object.
%
%   Algorithm:
%     1. Filter selectivity: event-level and pixel-level retention
%        statistics, polarity analysis of retained events.
%     2. Surface quality: Sobel gradient for edge sharpness,
%        intensity statistics for contrast and noise floor.
%     3. Temporal coherence: SSIM and MAE between consecutive frames.
%     4. Spatial coherence: connected component analysis of the
%        binary filter mask.
%
%   Notes:
%     - All metrics are computed from data already available in
%       main.m after the filter and accumulator stages. No
%       additional passes over the event stream are required.
%     - The noise_floor metric uses a heuristic for background
%       detection (zero event count). For more rigorous evaluation,
%       use ground truth object masks from the EVOS position data.
%     - Temporal SSIM uses MATLAB's built-in ssim() function from
%       the Image Processing Toolbox, which is already a dependency.
%     - For the IEI-ATS signed output, the neutral midpoint is 0.5
%       (zero surface maps to mid-gray via the symmetric sigmoid).
%       For all unsigned accumulators, the neutral midpoint is 0.
%       The noise_floor and surface_fill metrics account for this
%       by detecting the output type from the surface statistics.
%     - Connected component analysis uses bwconncomp on a binary
%       version of the filter mask (mask > 0).
%
%   Example:
%     % Inside main.m, after filtering and accumulation:
%     frame_metrics(frameIndex) = metrics.computeFrameMetrics(...
%         filter_mask, sorted_x, sorted_y, sorted_t, p_signed, ...
%         imgSz, counts, normalized_output_frame, ...
%         prev_output_frame_for_metrics, n_passed, n_total);
%
%   See also: metrics.summarizeMetrics, ssim, imgradient, bwconncomp

    % ================================================================
    % 0. Handle edge cases
    % ================================================================
    if isempty(sorted_x) || n_total == 0
        m = emptyMetrics();
        return;
    end

    % Detect surface type: signed (IEI-ATS, midpoint 0.5) vs.
    % unsigned (all others, midpoint 0). Heuristic: if > 30% of
    % nonzero pixels are in the range [0.45, 0.55], it is signed.
    nz_vals = normalized_output_frame(normalized_output_frame ~= 0);
    if ~isempty(nz_vals) && ...
            mean(nz_vals > 0.45 & nz_vals < 0.55) > 0.30
        midpoint = 0.5;
    else
        midpoint = 0.0;
    end

    % ================================================================
    % 1. Filter selectivity metrics
    % ================================================================

    % --- Signal Retention Rate (event-level) ---
    m.srr = n_passed / max(n_total, 1);

    % --- Pixel-level retention ---
    active_pixels = counts > 0;
    retained_pixels = (filter_mask > 0) & active_pixels;

    m.n_active_pixels  = sum(active_pixels(:));
    m.n_retained_pixels = sum(retained_pixels(:));
    m.pixel_retention  = m.n_retained_pixels / ...
        max(m.n_active_pixels, 1);

    % --- Polarity balance of retained events ---
    % Look up the filter mask value at each event location and
    % select events at retained pixels.
    linear_idx = sub2ind(imgSz, sorted_x(:), sorted_y(:));
    event_retained = filter_mask(linear_idx) > 0;

    p_retained = p_signed(event_retained);
    if ~isempty(p_retained)
        m.polarity_balance = abs(mean(p_retained));
    else
        m.polarity_balance = NaN;
    end

    m.n_passed = n_passed;
    m.n_total  = n_total;

    % ================================================================
    % 2. Surface quality metrics
    % ================================================================

    % --- Edge sharpness and mean gradient ---
    S = double(normalized_output_frame);

    [Gmag, ~] = imgradient(S, 'sobel');

    m.mean_gradient = mean(Gmag(:));

    % Edge sharpness: mean gradient at edge pixels.
    % Detect edges as pixels above the 90th percentile of gradient.
    if any(Gmag(:) > 0)
        edge_th = prctile(Gmag(Gmag > 0), 90);
        edge_mask = Gmag >= edge_th;
        m.edge_sharpness = mean(Gmag(edge_mask));
    else
        m.edge_sharpness = 0;
    end

    % --- Contrast (IQR of nonzero surface values) ---
    % Deviation from midpoint captures usable dynamic range.
    deviation = abs(S(S ~= 0) - midpoint);
    if numel(deviation) >= 4
        q = prctile(deviation, [25 75]);
        m.contrast_iqr = q(2) - q(1);
    else
        m.contrast_iqr = 0;
    end

    % --- Noise floor ---
    % Background: pixels with zero event count in this frame.
    bg_mask = counts == 0;
    if any(bg_mask(:))
        bg_vals = S(bg_mask);
        m.noise_floor = mean(abs(bg_vals - midpoint));
    else
        m.noise_floor = NaN;
    end

    % --- Surface fill ---
    % Fraction of pixels with non-trivial value
    m.surface_fill = mean(abs(S(:) - midpoint) > 0.01);

    % ================================================================
    % 3. Temporal coherence metrics
    % ================================================================
    if ~isempty(prev_output_frame) && any(prev_output_frame(:) ~= 0)
        S_prev = double(prev_output_frame);

        % SSIM (requires Image Processing Toolbox)
        % DynamicRange set to 1.0 since surfaces are in [0, 1].
        m.temporal_ssim = ssim(S, S_prev, 'DynamicRange', 1.0);

        % Mean absolute error
        m.temporal_mae = mean(abs(S(:) - S_prev(:)));
    else
        m.temporal_ssim = NaN;
        m.temporal_mae  = NaN;
    end

    % ================================================================
    % 4. Spatial coherence metrics
    % ================================================================
    binary_mask = filter_mask > 0;

    if any(binary_mask(:))
        cc = bwconncomp(binary_mask, 8);
        cluster_sizes = cellfun(@numel, cc.PixelIdxList);

        m.n_clusters          = cc.NumObjects;
        m.mean_cluster_size   = mean(cluster_sizes);
        m.largest_cluster_frac = max(cluster_sizes) / ...
            sum(cluster_sizes);
    else
        m.n_clusters          = 0;
        m.mean_cluster_size   = 0;
        m.largest_cluster_frac = 0;
    end

end

% ================================================================
% Local function: empty metrics struct (for skip frames)
% ================================================================
function m = emptyMetrics()
% EMPTYMETRICS  Return a metrics struct with NaN/zero fields.
%   Used when a frame has no events and metrics cannot be computed.

    m.srr                  = NaN;
    m.pixel_retention      = NaN;
    m.polarity_balance     = NaN;
    m.n_passed             = 0;
    m.n_total              = 0;
    m.n_active_pixels      = 0;
    m.n_retained_pixels    = 0;
    m.edge_sharpness       = NaN;
    m.mean_gradient        = NaN;
    m.contrast_iqr         = NaN;
    m.noise_floor          = NaN;
    m.surface_fill         = NaN;
    m.temporal_ssim        = NaN;
    m.temporal_mae         = NaN;
    m.n_clusters           = 0;
    m.mean_cluster_size    = 0;
    m.largest_cluster_frac = 0;
end
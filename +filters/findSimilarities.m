function [similarity_score, cv_map, regularity_map] = findSimilarities(sorted_x, sorted_y, std_map, mean_map, imgSz)
% findSimilarities Per-pixel inter-event-interval regularity scores.
%
%   [similarity_score, cv_map, regularity_map] = findSimilarities( ...
%       sorted_x, sorted_y, std_map, mean_map, imgSz)
%
%   Scores how temporally regular each pixel's firing is, then gates that
%   score by IEI magnitude so that "regular AND reasonably fast" pixels
%   (signal) are separated from "regular but glacially slow" pixels
%   (the long-gap backfill noise that inflates the mean IEI).
%
%   Regularity functional:
%     Uses the coefficient of local variation
%       LV = (1/(n-1)) * sum_i 3 (T_i - T_{i+1})^2 / (T_i + T_{i+1})^2
%     (Shinomoto, Shima & Tanji, "Differences in spiking patterns among
%      cortical neurons," Neural Computation 15(12):2823-2842, 2003,
%      doi:10.1162/089976603322518759). LV = 0 perfectly regular,
%      LV ~ 1 Poisson, LV > 1 bursty/irregular.
%
%     IMPORTANT SCOPE: std_map/mean_map are aggregate per-pixel IEI
%     statistics, not the raw interval sequence, so true multi-interval
%     LV cannot be reconstructed here for dense pixels. For the two-IEI
%     backfill regime (counts == 1 or 2 in updateWindowedIEI) there is
%     exactly one interval pair, and LV collapses to a closed form in the
%     CV alone:  LV = 1.5 * CV^2  (verified algebraically and numerically;
%     e1 = m+d, e2 = m-d => CV^2 = 2 d^2/m^2, LV = 3 d^2/m^2). For pixels
%     with >= 3 intervals this expression is a MONOTONIC remap of CV, not
%     their true LV; it is exact precisely where the saturation pile-up
%     occurred (the 2-IEI pixels).
%
%   Why not 1/(1+CV): for two intervals CV saturates at sqrt(2), so the
%   old score 1/(1+CV) was floored at 1/(1+sqrt2) ~ 0.414 and crushed all
%   irregular 2-IEI pixels into [0.41, 0.55], piling them on the decision
%   boundary. LV (range [0, 3) at n=2) does not saturate and spreads them.
%
%   Magnitude gate (Rule 3, steps 4-5 of the original spec): a soft
%   threshold at the median observed IEI penalises very large mean IEI,
%   so a noise pixel that fires twice across a long silence with two
%   similar gaps (CV ~ 0, LV ~ 0, "perfectly regular") is still rejected
%   on magnitude. Without this gate, scale-invariant LV/CV cannot catch
%   balanced-but-huge-gap noise, and the inflating mean is never penalised.
%
%   Inputs:
%     sorted_x, sorted_y - event coords (sorted by linear pixel index).
%     std_map, mean_map  - per-pixel IEI std and mean maps (imgSz);
%                          0 at pixels with no observed IEI.
%     imgSz              - [nrows, ncols].
%
%   Outputs:
%     similarity_score - N-by-1 per-event regularity score in [0,1],
%                        high = regular & reasonably fast.
%     cv_map           - per-pixel CV (raw), for visualisation/tuning.
%     regularity_map   - per-pixel combined score (LV-based * magnitude).
%

    % --- tuning -----------------------------------------------------
    mag_sharpness = 2;          % soft-threshold exponent (Hill slope)

    % --- 1. observation mask ----------------------------------------
    obs_mask = mean_map > 0;

    % --- 2. per-pixel CV --------------------------------------------
    cv_map = std_map ./ max(mean_map, eps);
    cv_map(~obs_mask) = Inf;            % unobserved -> score 0 downstream

    % --- 3. LV-based regularity -------------------------------------
    %   LV = 1.5 * CV^2  (exact for the 2-IEI backfill regime; monotonic
    %   CV remap for >= 3 IEIs). reg in (0,1], high = regular.
    lv_map  = 1.5 .* cv_map.^2;         % Inf where unobserved
    reg_map = 1 ./ (1 + lv_map);        % -> 0 where unobserved

    % --- 4. IEI magnitude gate (soft threshold at median IEI) -------
    if any(obs_mask(:))
        med_iei = median(mean_map(obs_mask));
    else
        med_iei = eps;
    end
    %   1 for fast pixels, 0.5 at the median, -> 0 for very slow pixels.
    mag_map = 1 ./ (1 + (mean_map ./ max(med_iei, eps)).^mag_sharpness);
    mag_map(~obs_mask) = 0;

    % --- 5. combine -------------------------------------------------
    regularity_map = reg_map .* mag_map;

    % --- 6. per-event lookup ----------------------------------------
    linear_idx = sub2ind(imgSz, sorted_x(:), sorted_y(:));
    similarity_score = regularity_map(linear_idx);

end
% 
% function [similarity_score, cv_map, regularity_map] = findSimilarities(sorted_x, sorted_y, std_map, mean_map, imgSz)
% % findSimilarities Compute per-pixel inter-event-interval regularity scores.
% %
% %   [similarity_score, cv_map, regularity_map] = findSimilarities( ...
% %       sorted_x, sorted_y, std_map, mean_map, imgSz)
% %
% %   Measures how temporally regular each pixel's firing pattern is by
% %   computing the Coefficient of Variation (CV = sigma/mu) directly from
% %   the EMA-tracked per-pixel IEI statistics. Real signal pixels fire at
% %   consistent intervals, yielding low CV. Noise pixels fire erratically,
% %   yielding high CV. A pure-noise region where all pixels have similarly
% %   large IEIs may still produce low CV, so the score also incorporates
% %   IEI magnitude to disambiguate: only pixels with both low CV AND
% %   reasonable mean IEI receive high scores.
% %
% %   Inputs:
% %     sorted_x  - Vector of x-coordinates (row indices) for events in the
% %                 current frame, sorted by linear pixel index.
% %     sorted_y  - Vector of y-coordinates (column indices) for events.
% %     std_map   - 2D matrix (imgSz) of EMA-tracked per-pixel IEI standard
% %                 deviation. Pixels with no observed events should be 0.
% %     mean_map  - 2D matrix (imgSz) of EMA-tracked per-pixel IEI mean.
% %                 Pixels with no observed events should be 0.
% %     imgSz     - Two-element vector [nrows, ncols].
% %
% %   Outputs:
% %     similarity_score - N-by-1 vector of scores in [0, 1] for each input
% %                        event. High values indicate temporally regular
% %                        pixels (likely real). Low values indicate
% %                        irregular or inactive pixels (likely noise).
% %     cv_map           - 2D matrix (imgSz) of per-pixel coefficient of
% %                        variation values. Useful for visualization,
% %                        parameter tuning, and downstream feature
% %                        detection (e.g. Harris on CV field).
% %     regularity_map   - 2D matrix (imgSz) of the final combined score
% %                        before event lookup. Useful for visualization.
% %
% %   Algorithm:
% %     1. Build an observation mask from pixels with nonzero mean IEI.
% %     2. Compute per-pixel CV = std_map / mean_map.
% %     3. Convert CV to a regularity score in [0,1] via 1/(1 + CV).
% %     4. Compute an IEI magnitude score that penalizes very large
% %        intervals using a soft threshold at the median observed IEI.
% %     5. Combine: score = regularity * magnitude_score.
% %     6. Look up the combined score at each event location.
% %
% %   Notes:
% %     - This function uses the pre-computed EMA-tracked per-pixel IEI
% %       statistics (std_map, mean_map) rather than computing local
% %       spatial statistics via convolution. Per-pixel temporal CV is the
% %       principled metric for Rule 3: it directly measures whether a
% %       pixel's own firing pattern is regular over time. Spatial context
% %       is handled independently by Rule 1 (spatial density).
% %     - The observation mask uses mean_map > 0 to identify pixels that
% %       have received at least one event update.
% %
% 
%     % ----------------------------------------------------------------
%     % 1. Observation mask
%     % ----------------------------------------------------------------
%     obs_mask = mean_map > 0;
% 
%     % ----------------------------------------------------------------
%     % 2. Per-pixel Coefficient of Variation
%     % ----------------------------------------------------------------
%     cv_map = std_map ./ max(mean_map, eps);
% 
%     % Unobserved pixels: set CV to Inf so they map to score = 0
%     cv_map(~obs_mask) = Inf;
% 
%     % ----------------------------------------------------------------
%     % 3. Regularity score: low CV -> high score
%     %    Monotonic mapping from [0, inf) to (0, 1]:
%     %      CV = 0   -> score = 1.0   (perfectly uniform firing)
%     %      CV = 1   -> score = 0.5   (std equals mean)
%     %      CV -> inf -> score -> 0   (highly irregular)
%     % ----------------------------------------------------------------
%     regularity_map = 1 ./ (1 + cv_map);
% 
%     % ----------------------------------------------------------------
%     % 6. Look up per-event scores
%     % ----------------------------------------------------------------
%     linear_idx = sub2ind(imgSz, sorted_x(:), sorted_y(:));
%     similarity_score = regularity_map(linear_idx);
% 
% end
% 

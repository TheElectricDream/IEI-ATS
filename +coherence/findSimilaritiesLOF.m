function [similarity_score] = findSimilaritiesLOF(...
    x_valid, y_valid, cv_valid, imgSz)
% FINDSIMILARITIESLOF  Local Outlier Factor-style similarity scoring (experimental).
%
%   SIMILARITY_SCORE = FINDSIMILARITIESLOF(X_VALID, Y_VALID,
%   CV_VALID, IMGSZ) computes a simplified LOF-like score for each
%   event by comparing its local neighbourhood density to that of
%   its neighbours in a 3D normalized (x, y, feature) space.
%
%   Inputs:
%     x_valid  - [N x 1] Row coordinates (pixels).
%     y_valid  - [N x 1] Column coordinates (pixels).
%     cv_valid - [N x 1] Third feature dimension (e.g., CV or IEI).
%     imgSz    - [1 x 2] Image dimensions [nRows, nCols].
%
%   Outputs:
%     similarity_score - [N x 1] Scores > 0. Values near 1 indicate
%                        consistent local density (inlier). Values > 1
%                        indicate relatively sparser points (outlier).
%
%   Algorithm:
%     1. Normalize x, y to [0, 1] by imgSz. cv_valid used as-is.
%     2. Build KD-tree over the 3D points.
%     3. K=5 nearest neighbour search.
%     4. Compute mean distance for each point (inverse density).
%     5. Compute ratio: mean of neighbours' densities / own density.
%
%   Notes:
%     - This is an experimental alternative to findSimilarities.
%       It is not used in the main IEI-ATS pipeline.
%     - K=5 neighbours is hardcoded. For denser data, increase K.
%     - Coordinates: x = row, y = col.
%
%   See also: coherence.findSimilarities,
%             coherence.computeCoherenceMask

    % ----------------------------------------------------------------
    % 0. Normalize and build KD-tree
    % ----------------------------------------------------------------
    x_norm = single(x_valid(:) ./ imgSz(1));
    y_norm = single(y_valid(:) ./ imgSz(2));
    cv_norm = single(cv_valid(:));

    points3D = [x_norm, y_norm, cv_norm];

    tree = createns(points3D, 'NsMethod', 'kdtree', ...
        'Distance', 'euclidean', 'BucketSize', 50);

    % ----------------------------------------------------------------
    % 1. K=5 nearest neighbour search
    % ----------------------------------------------------------------
    [idx, D] = knnsearch(tree, points3D, 'K', 5);

    % ----------------------------------------------------------------
    % 2. Compute simplified LOF ratio
    % ----------------------------------------------------------------
    % Mean distance per point (inverse density proxy)
    mean_dist = mean(D, 2);
    mean_dist(mean_dist == 0) = eps;

    % Mean distance of each point's neighbours
    neighbor_mean_dists = mean_dist(idx);
    local_context_density = mean(neighbor_mean_dists, 2);

    % Ratio: if > 1, the point is sparser than its neighbours
    similarity_score = local_context_density ./ mean_dist;

end
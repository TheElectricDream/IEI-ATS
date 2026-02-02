function [similarity_score] = findSimilarities(x_valid, y_valid, cv_valid, imgSz)
    % This function finds the spatial neighbors of given points within a specified radius.
    % Inputs:
    %   x_valid - x-coordinates of events
    %   y_valid - y-coordinates of events
    %   t_valid - time stamp of events
    %   r_s - search radius for neighbors
    % Outputs:
    %   neighbor_db - a cell array where each cell contains the indices of neighbors for each point

    % Normalize to [0,1]
    x_norm = single(x_valid ./ imgSz(1)); 
    y_norm = single(y_valid ./ imgSz(2)); 
    cv_norm = single(cv_valid);

    % Stack the data
    points3D = [x_norm(:), y_norm(:), cv_norm(:)];

    % Create the tree object for the rangesearch function
    tree = createns(points3D, 'NsMethod', 'kdtree', 'Distance','euclidean', 'BucketSize', 50);

    % Use MATLABs rangesearch to find the nearest neighbours
    [idx, D] = knnsearch(tree, points3D, 'K', 200);

    % 2. Compute Mean Distance for every point (Inverse Density)
    % D is N x K. We take the mean across the rows (neighbors).
    mean_dist = mean(D, 2); 

    % Handle divide-by-zero if duplicates exist
    mean_dist(mean_dist == 0) = eps; 

    % 3. Look up the Mean Distance of the Neighbors
    % idx is N x K. We want mean_dist(idx).
    neighbor_mean_dists = mean_dist(idx); % Now N x K matrix

    % 4. Compute the Ratio (Simplified LOF)
    % "Am I as close to my neighbors as they are to theirs?"
    % If Ratio > 1, I am essentially "further out" than expected (Outlier)
    % If Ratio ~ 1, I fit in well (Similar)
    local_context_density = mean(neighbor_mean_dists, 2);
    similarity_score = local_context_density ./ mean_dist;

end
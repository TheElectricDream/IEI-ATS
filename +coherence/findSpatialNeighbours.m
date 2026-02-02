function [neighbor_db, D] = findSpatialNeighbours(x_valid, y_valid, t_valid, r_s, imgSz, t_interval)
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
    t_norm = single((t_valid-min(t_valid)) ./ t_interval);

    % Stack the data
    points3D = [x_norm(:), y_norm(:), t_norm(:)];

    % Create the tree object for the rangesearch function
    tree = createns(points3D, 'NsMethod', 'kdtree', 'Distance','euclidean', 'BucketSize', 50);

    % Use MATLABs rangesearch to find the nearest neighbours
    [idx, D] = rangesearch(tree, points3D, r_s);

    % Return the nearest neighbours
    neighbor_db = idx';
end
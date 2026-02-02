function [similarity_score] = findSimilaritiesTemporal(current_xy, current_vals, prev_xy, prev_vals)
% FINDSIMILARITIES Temporal Consistency Check
% Checks if current event intervals are similar to the intervals of 
% spatial neighbors from the PREVIOUS frame.
%
% Inputs:
%   current_xy   - [N x 2] matrix of (x,y) coordinates for current frame
%   current_vals - [N x 1] vector of values (e.g., mean time interval)
%   prev_xy      - [M x 2] matrix of (x,y) coordinates from previous frame
%   prev_vals    - [M x 1] vector of values from previous frame

    % Handle initialization case (Frame 1)
    if isempty(prev_xy)
        similarity_score = ones(size(current_vals)); 
        return;
    end

    % 1. Create Search Tree on PREVIOUS Frame Data
    % We search in the *previous* cloud to find neighbors for *current* points
    tree = createns(prev_xy, 'NsMethod', 'kdtree');

    % 2. Find 10 Nearest Neighbors in Previous Frame
    % idx contains indices into prev_vals
    K = 200;

    % Ensure we don't ask for more neighbors than exist
    K = min(K, size(prev_xy, 1)); 
    [idx, ~] = knnsearch(tree, current_xy, 'K', K);

    % 3. Retrieve Neighbor Values
    % idx is [N x K]. We get the corresponding values from history.
    neighbor_vals = prev_vals(idx); 

    % 4. Compute Similarity
    % Calculate the mean interval of the historical neighbors
    mean_neighbor_val = mean(neighbor_vals, 2, 'omitnan');

    % Calculate percent difference (or absolute difference)
    % A score of 1 means perfect match. 0 means very different.
    % Using a Lorentzian-like kernel for similarity: 1 / (1 + diff)
    
    % Option A: Absolute Difference
    % diff_val = abs(current_vals - mean_neighbor_val);
    
    % Option B: Relative Difference (Scale Invariant)
    % Good for time intervals which can vary largely in magnitude
    diff_val = abs(current_vals - mean_neighbor_val) ./ (mean_neighbor_val + eps);
    
    % Convert difference to similarity score [0, 1]
    % If diff is 0 (perfect match), score is 1.
    % If diff is 1 (100% variance), score is 0.5.
    similarity_score = 1 ./ (1 + diff_val);

end

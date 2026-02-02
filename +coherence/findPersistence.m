function [closestRows, closestCols, sortedIdx,...
            t_row, t_col, t_val, ...
            b_rows, b_cols, b_vals, minDist] = findPersistence(mapA, mapB, imgSz, targetRow, targetCol)

    % Define Image Dimensions 
    height = imgSz(1);
    width = imgSz(2);
    
    % Get all valid points from mapB
    % We need Rows, Cols, and the Values at those points
    [b_rows, b_cols] = find(~isnan(mapB));
    b_vals = mapB(sub2ind([height, width], b_rows, b_cols));
    
    % Get the point from mapA
    t_row = targetRow;
    t_col = targetCol;
    t_val = mapA(t_row, t_col);
    
    % Normalize Candidates (mapB points) to range [0, 1]
    % Spatial Normalization
    n_b_rows = (b_rows - 1) / (height - 1);
    n_b_cols = (b_cols - 1) / (width - 1);
    
    % Time Normalization 
    min_v = min(b_vals);
    max_v = max(b_vals);
    range_v = max_v - min_v;
    
    % Handle case where all values are identical (div by zero)
    if range_v == 0, range_v = 1; end 
    
    n_b_vals = (b_vals - min_v) / range_v;
    
    % Normalize Target (mapA point) using mapB's scaling factors
    n_t_row = (t_row - 1) / (height - 1);
    n_t_col = (t_col - 1) / (width - 1);
    n_t_val = (t_val - min_v) / range_v;
    
    % Calculate 3D Euclidean Distance (Squared)
    % D^2 = (dx)^2 + (dy)^2 + (d_val)^2
    sq_dist = (n_b_rows - n_t_row).^2 + ...
              (n_b_cols - n_t_col).^2 + ...
              (n_b_vals - n_t_val).^2;
    
    % Find 10 Closest
    [minDist, sortedIdx] = mink(sq_dist, 1);
    
    % Return Original Coordinates 
    closestRows = b_rows(sortedIdx);
    closestCols = b_cols(sortedIdx);
    % closestVals = b_vals(sortedIdx);
    
end
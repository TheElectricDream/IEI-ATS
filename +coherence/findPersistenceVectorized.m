function minDist = findPersistenceVectorized(mapA, mapB, imgSz, rowsExists, colsExists)

    height = imgSz(1);
    width  = imgSz(2);
    
    % --- Prepare Candidates (Reference Set) from mapB ---
    [brows, bcols] = find(~isnan(mapB));
    if isempty(brows)
        minDist = inf(numel(rowsExists), 1);
        return;
    end
    bvals = mapB(sub2ind([height, width], brows, bcols));
    
    % Normalize candidates
    minv   = min(bvals);
    maxv   = max(bvals);
    rangev = maxv - minv;
    if rangev == 0, rangev = 1; end
    
    nbrows = (brows - 1) / (height - 1);
    nbcols = (bcols - 1) / (width  - 1);
    nbvals = (bvals - minv) / rangev;
    
    % Create Reference matrix [N_candidates x 3]
    X = [nbrows, nbcols, nbvals];
    
    % --- Prepare Targets (Query Set) from mapA ---
    tvals = mapA(sub2ind([height, width], rowsExists, colsExists));
    
    ntrows = (rowsExists - 1) / (height - 1);
    ntcols = (colsExists - 1) / (width  - 1);
    ntvals = (tvals - minv) / rangev;
    
    % Create Query matrix [N_targets x 3]
    Y = [ntrows(:), ntcols(:), ntvals(:)];
    
    % --- Find Nearest Neighbors ---
    % knnsearch finds the nearest neighbor in X for each point in Y
    [~, minDist] = knnsearch(X, Y); 
end

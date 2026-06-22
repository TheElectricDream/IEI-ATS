function refine_lines(map)
    % refine_lines - Compares 3 line-tightening algorithms on a binary map
    % Input: map - A 2D binary matrix (e.g., 480x640)

    % Ensure the map is a logical binary format
    if ~islogical(map)
        map = map > 0;
    end

    [height, width] = size(map);

    % Extract (x, y) coordinates from the map matrix
    % MATLAB's find returns (row, col) which maps to (y, x)
    [y_coords, x_coords] = find(map);
    points = [x_coords, y_coords];
    num_points = size(points, 1);

    if num_points < 2
        disp('Not enough points in the map to form lines.');
        return;
    end

    % Thresholds
    snap_threshold = 1.0;
    
    %% ==========================================
    % APPROACH 1: HOUGH-SNAP
    % ==========================================
    % Requires Image Processing Toolbox
    [H, theta, rho] = hough(map);
    P = houghpeaks(H, 1, 'threshold', ceil(0.3 * max(H(:))));
    lines = houghlines(map, theta, rho, P, 'FillGap', 20, 'MinLength', 50);

    hough_snapped = [];
    if ~isempty(lines)
        % Extract line segment endpoints
        x1 = lines(1).point1(1); y1 = lines(1).point1(2);
        x2 = lines(1).point2(1); y2 = lines(1).point2(2);

        % Calculate Line Equation: Ax + By + C = 0
        A = y2 - y1;
        B = x1 - x2;
        C = x2*y1 - x1*y2;
        norm_sq = A^2 + B^2;

        for i = 1:num_points
            pt = points(i, :);
            dist = abs(A*pt(1) + B*pt(2) + C) / sqrt(norm_sq);
            
            if dist <= snap_threshold
                % Orthogonal projection onto the line
                x_new = pt(1) - A * (A*pt(1) + B*pt(2) + C) / norm_sq;
                y_new = pt(2) - B * (A*pt(1) + B*pt(2) + C) / norm_sq;
                hough_snapped = [hough_snapped; x_new, y_new]; %#ok<AGROW>
            end
        end
    end

    %% ==========================================
    % APPROACH 2: RANSAC (Custom Implementation)
    % ==========================================
    max_iters = 1000;
    best_inliers = false(num_points, 1);
    best_line = [0, 0, 0]; % A, B, C

    for i = 1:max_iters
        % Randomly sample 2 distinct points
        idx = randperm(num_points, 2);
        p1 = points(idx(1), :);
        p2 = points(idx(2), :);

        A = p2(2) - p1(2);
        B = p1(1) - p2(1);
        C = p2(1)*p1(2) - p1(1)*p2(2);
        norm_val = sqrt(A^2 + B^2);

        if norm_val < 1e-5; continue; end

        % Measure orthogonal distances of all points to the candidate line
        dists = abs(A*points(:,1) + B*points(:,2) + C) / norm_val;
        inliers = dists < snap_threshold;

        % Keep the line with the most inliers
        if sum(inliers) > sum(best_inliers)
            best_inliers = inliers;
            best_line = [A, B, C];
        end
    end

    ransac_snapped = [];
    if sum(best_inliers) > 0
        A = best_line(1); B = best_line(2); C = best_line(3);
        norm_sq = A^2 + B^2;
        
        inlier_pts = points(best_inliers, :);
        for i = 1:size(inlier_pts, 1)
            pt = inlier_pts(i, :);
            x_new = pt(1) - A * (A*pt(1) + B*pt(2) + C) / norm_sq;
            y_new = pt(2) - B * (A*pt(1) + B*pt(2) + C) / norm_sq;
            ransac_snapped = [ransac_snapped; x_new, y_new]; %#ok<AGROW>
        end
    end

    %% ==========================================
    % APPROACH 3: MOVING LEAST SQUARES (MLS)
    % ==========================================
    radius = 10.0;
    min_neighbors = 20;
    mls_snapped = [];

    for i = 1:num_points
        % Find local neighbors using basic distance formula
        dists = sqrt((points(:,1) - points(i,1)).^2 + (points(:,2) - points(i,2)).^2);
        neighbors = points(dists <= radius, :);

        if size(neighbors, 1) >= min_neighbors
            % Center the local neighborhood data
            mu = mean(neighbors, 1);
            centered = neighbors - mu;
            
            % Perform Singular Value Decomposition (SVD) for local PCA
            [~, ~, V] = svd(centered, 'econ');
            dir_vec = V(:, 1)'; % Direction of maximum variance (local line)

            % Project the current point onto the local trendline
            pt_centered = points(i, :) - mu;
            proj_len = dot(pt_centered, dir_vec);
            snapped_pt = mu + proj_len * dir_vec;

            mls_snapped = [mls_snapped; snapped_pt]; %#ok<AGROW>
        end
    end

    %% ==========================================
    % 4. PLOTTING THE COMPARISON
    % ==========================================
    figure('Name', 'Line Tightener Algorithms', 'Position', [100, 100, 1000, 700]);
    
    % Original Map
    subplot(2, 2, 1);
    scatter(points(:,1), points(:,2), 10, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
    title('Original Map (Signal + Noise)');
    axis([0 width 0 height]);
    set(gca, 'YDir', 'reverse'); % Flip Y axis to match image coordinates

    % Hough-Snap
    subplot(2, 2, 2);
    if ~isempty(hough_snapped)
        scatter(hough_snapped(:,1), hough_snapped(:,2), 10, 'b', 'filled');
    end
    title('Approach 1: Hough-Snap');
    axis([0 width 0 height]);
    set(gca, 'YDir', 'reverse');

    % RANSAC
    subplot(2, 2, 3);
    if ~isempty(ransac_snapped)
        scatter(ransac_snapped(:,1), ransac_snapped(:,2), 10, 'g', 'filled');
    end
    title('Approach 2: RANSAC');
    axis([0 width 0 height]);
    set(gca, 'YDir', 'reverse');

    % MLS
    subplot(2, 2, 4);
    if ~isempty(mls_snapped)
        scatter(mls_snapped(:,1), mls_snapped(:,2), 10, 'r', 'filled');
    end
    title('Approach 3: Local MLS');
    axis([0 width 0 height]);
    set(gca, 'YDir', 'reverse');
end
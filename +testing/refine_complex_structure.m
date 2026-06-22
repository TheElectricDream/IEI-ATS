function refine_complex_structure(map)
    % refine_complex_structure - Uses Sequential RANSAC to tighten multiple lines
    % Input: map - A 2D binary matrix
    
    if nargin < 1
        % Generate a synthetic complex test map (a triangle with noise) if no map is provided
        map = generate_test_map();
    end

    if ~islogical(map)
        map = map > 0;
    end

    [height, width] = size(map);
    [y_coords, x_coords] = find(map);
    remaining_points = [x_coords, y_coords];
    
    % --- ALGORITHM PARAMETERS ---
    snap_threshold = 15.0;       % Max distance to snap points to a line
    min_points_for_line = 1;    % Stop if we can't find a line with at least this many points
    max_lines = 2000;              % Safety break to prevent infinite loops
    ransac_iters = 500;
    
    final_snapped_cloud = [];
    line_count = 0;

    % --- SEQUENTIAL RANSAC LOOP ---
    while size(remaining_points, 1) > min_points_for_line && line_count < max_lines
        
        num_points = size(remaining_points, 1);
        best_inliers = false(num_points, 1);
        best_line = [0, 0, 0]; % A, B, C

        % 1. Run RANSAC to find the BEST current line
        for i = 1:ransac_iters
            idx = randperm(num_points, 2);
            p1 = remaining_points(idx(1), :);
            p2 = remaining_points(idx(2), :);

            A = p2(2) - p1(2);
            B = p1(1) - p2(1);
            C = p2(1)*p1(2) - p1(1)*p2(2);
            norm_val = sqrt(A^2 + B^2);

            if norm_val < 1e-5; continue; end

            dists = abs(A*remaining_points(:,1) + B*remaining_points(:,2) + C) / norm_val;
            inliers = dists < snap_threshold;

            if sum(inliers) > sum(best_inliers)
                best_inliers = inliers;
                best_line = [A, B, C];
            end
        end

        % 2. Check if the best found line has enough points to be valid
        if sum(best_inliers) < min_points_for_line
            break; % No more dense lines exist, just noise left.
        end

        line_count = line_count + 1;

        % 3. Snap the inliers to the detected line
        A = best_line(1); B = best_line(2); C = best_line(3);
        norm_sq = A^2 + B^2;
        
        inlier_pts = remaining_points(best_inliers, :);
        for i = 1:size(inlier_pts, 1)
            pt = inlier_pts(i, :);
            x_new = pt(1) - A * (A*pt(1) + B*pt(2) + C) / norm_sq;
            y_new = pt(2) - B * (A*pt(1) + B*pt(2) + C) / norm_sq;
            final_snapped_cloud = [final_snapped_cloud; x_new, y_new]; %#ok<AGROW>
        end

        % 4. REMOVE the inliers from the pool so the next iteration finds the *next* line
        remaining_points(best_inliers, :) = [];
    end

    % --- PLOT RESULTS ---
    [orig_y, orig_x] = find(map);
    
    figure('Name', 'Sequential Line Tightener', 'Position', [100, 100, 900, 450]);
    
    subplot(1, 2, 1);
    scatter(orig_x, orig_y, 10, [0.6 0.6 0.6], 'filled');
    title('Original Complex Map');
    axis([0 width 0 height]); set(gca, 'YDir', 'reverse');
    
    subplot(1, 2, 2);
    if ~isempty(final_snapped_cloud)
        scatter(final_snapped_cloud(:,1), final_snapped_cloud(:,2), 10, 'b', 'filled');
    end
    title(sprintf('Tightened (Found %d Lines)', line_count));
    axis([0 width 0 height]); set(gca, 'YDir', 'reverse');
end

function map = generate_test_map()
    % Helper function to create a noisy triangle to test multi-line logic
    map = false(480, 640);
    
    % Line 1 (Bottom edge)
    x1 = round(linspace(100, 500, 200)); y1 = repmat(380, 1, 200) + round(5*randn(1,200));
    % Line 2 (Left edge)
    x2 = round(linspace(100, 300, 200)) + round(5*randn(1,200)); y2 = round(linspace(380, 100, 200));
    % Line 3 (Right edge)
    x3 = round(linspace(300, 500, 200)) + round(5*randn(1,200)); y3 = round(linspace(100, 380, 200));
    
    all_x = [x1, x2, x3, randi([1, 640], 1, 200)]; % + random background noise
    all_y = [y1, y2, y3, randi([1, 480], 1, 200)];
    
    valid = all_x >= 1 & all_x <= 640 & all_y >= 1 & all_y <= 480;
    idx = sub2ind([480, 640], all_y(valid), all_x(valid));
    map(idx) = true;
end
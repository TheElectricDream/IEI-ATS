function [] = visualizePersistance(t_row, t_col, t_val, ...
    b_rows, b_cols, b_vals, sortedIdx)
    
    % 1. Create a Figure
    figure('Color', 'white');
    clf; hold on; grid on; rotate3d on;
    
    % 2. Separate "Closest" from "Rest"
    % Create a logical mask for the top 10 indices
    isClosest = false(size(b_rows));
    isClosest(sortedIdx) = true;
    
    % Split the data
    rest_rows = b_rows(~isClosest);
    rest_cols = b_cols(~isClosest);
    rest_vals = b_vals(~isClosest);
    
    closest_rows = b_rows(isClosest);
    closest_cols = b_cols(isClosest);
    closest_vals = b_vals(isClosest);
    
    % 3. Plot "Rest of mapB" (Context)
    % Style: Small blue dots, high transparency to reduce clutter
    h1 = scatter3(rest_cols(rest_vals~=0), rest_rows(rest_vals~=0), rest_vals(rest_vals~=0), 10, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.2); 
    
    % 4. Plot "10 Closest Points"
    % Style: Large green circles, solid
    h2 = scatter3(closest_cols, closest_rows, closest_vals, 100, 'g', 'filled', ...
        'MarkerEdgeColor', 'k'); 
    
    % 5. Plot "Target Point from mapA"
    % Style: Large red pentagram (star)
    h3 = scatter3(t_col, t_row, t_val, 200, 'r', 'p', 'filled', ...
        'MarkerEdgeColor', 'k');
    
    % 6. Formatting
    xlabel('Column Index (X)');
    ylabel('Row Index (Y)');
    zlabel('Time Value');
    title('3D Proximity Search: Space + Time');
    view(3); 
    
    % Optional: Draw lines connecting target to the 10 closest
    % This helps visually measure the "distance"
    for i = 1:length(closest_rows)
        plot3([t_col, closest_cols(i)], ...
              [t_row, closest_rows(i)], ...
              [t_val, closest_vals(i)], 'k-', 'LineWidth', 0.5);
    end
    set(gca, 'FontSize', 16);
    hold off;

end

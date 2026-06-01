function [] = visualizePersistance(t_row, t_col, t_val, ...
    b_rows, b_cols, b_vals, sortedIdx)
% VISUALIZEPERSISTANCE  3D plot of persistence neighbourhood search.
%
%   VISUALIZEPERSISTANCE(T_ROW, T_COL, T_VAL, B_ROWS, B_COLS,
%   B_VALS, SORTEDIDX) visualizes a target point from mapA and its
%   nearest neighbours from mapB, with lines connecting them.
%   Includes a zoomed-in inset for clarity.

    % --- 1. Professional Style Definitions ---
    % Use muted, publication-ready RGB triplets instead of bright standard colors
    colorContext = [0.65, 0.65, 0.65]; % Subtle grey
    colorMatch   = [0.00, 0.45, 0.74]; % Deep blue
    colorTarget  = [0.85, 0.33, 0.10]; % Burnt orange
    
    % Initialize figure with a specific size to accommodate the inset
    figure('Color', 'white', 'Position', [100, 100, 900, 650]);
    set(gcf, 'DefaultTextFontName', 'Times New Roman', ...
             'DefaultAxesFontName', 'Times New Roman');

    % --- 2. Data Preparation ---
    isClosest = false(size(b_rows));
    isClosest(sortedIdx) = true;
    rest_rows = b_rows(~isClosest);
    rest_cols = b_cols(~isClosest);
    rest_vals = b_vals(~isClosest);
    closest_rows = b_rows(isClosest);
    closest_cols = b_cols(isClosest);
    closest_vals = b_vals(isClosest);
    active = rest_vals ~= 0;

    % --- 3. Main Axes ---
    axMain = axes('Position', [0.1, 0.1, 0.85, 0.85]);
    hold(axMain, 'on'); grid(axMain, 'on'); rotate3d(axMain, 'on');
    
    % Context points (Transparent grey)
    scatter3(axMain, rest_cols(active), rest_rows(active), rest_vals(active), ...
        15, colorContext, 'filled', 'MarkerFaceAlpha', 0.15);
        
    % Closest points (Medium blue circles with clean white edges)
    scatter3(axMain, closest_cols, closest_rows, closest_vals, ...
        60, colorMatch, 'o', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
        
    % Target point (Burnt orange diamond with black edge)
    scatter3(axMain, t_col, t_row, t_val, ...
        80, colorTarget, 'd', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        
    % Connecting lines (Dark grey, semi-transparent)
    for i = 1:length(closest_rows)
        plot3(axMain, [t_col, closest_cols(i)], ...
                      [t_row, closest_rows(i)], ...
                      [t_val, closest_vals(i)], ...
            '-', 'Color', [0.3 0.3 0.3 0.6], 'LineWidth', 1.2);
    end
    
    xlabel(axMain, 'X [px]');
    ylabel(axMain, 'Y [px]');
    zlabel(axMain, '\rho_{norm} [-]');
    
    view(axMain, 3);
    axis(axMain, 'tight');
    box(axMain, 'on');
    set(axMain, 'FontSize', 14);

    % --- 4. Zoomed-In Inset Axes ---
    % Positioned in the upper right corner
    axInset = axes('Position', [0.65, 0.60, 0.25, 0.25]);
    hold(axInset, 'on'); grid(axInset, 'on'); box(axInset, 'on');
    
    % Re-plot the elements in the inset (slightly higher alpha for context points)
    scatter3(axInset, rest_cols(active), rest_rows(active), rest_vals(active), ...
        15, colorContext, 'filled', 'MarkerFaceAlpha', 0.3);
    scatter3(axInset, closest_cols, closest_rows, closest_vals, ...
        60, colorMatch, 'o', 'filled', 'MarkerEdgeColor', 'w');
    scatter3(axInset, t_col, t_row, t_val, ...
        80, colorTarget, 'd', 'filled', 'MarkerEdgeColor', 'k');
        
    for i = 1:length(closest_rows)
        plot3(axInset, [t_col, closest_cols(i)], ...
                       [t_row, closest_rows(i)], ...
                       [t_val, closest_vals(i)], ...
            '-', 'Color', [0.3 0.3 0.3 0.6], 'LineWidth', 1.2);
    end
    
    % Set spatial limits to tightly frame the target point
    % (You may need to tweak these margin variables depending on your array dimensions)
    xyMargin = 15; 
    zMargin = max(closest_vals) - min(closest_vals) + 0.1; 
    
    xlim(axInset, [t_col - xyMargin, t_col + xyMargin]);
    ylim(axInset, [t_row - xyMargin, t_row + xyMargin]);
    zlim(axInset, [t_val - zMargin, t_val + zMargin]);
    
    view(axInset, 3);
    title(axInset, 'Local Neighborhood', 'FontSize', 16, 'FontWeight', 'normal');
    
    % Remove tick labels on the inset to keep it looking clean
    set(axInset, 'FontSize', 10, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
    exportgraphics(gcf,'Temporal-Persistance-Sample-Nom-Rot.pdf')
end
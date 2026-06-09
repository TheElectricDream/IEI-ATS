function [] = showPersistance(t_row, t_col, t_val, b_rows, b_cols, b_vals, sortedIdx)
% VISUALIZEPERSISTANCE  3D plot of persistence neighbourhood search.
% Styled to perfectly match showScatterPlotOfRuleMaps3D aesthetics with a zoomed inset.

    % 1. Normalize coordinates to [0 1] 
    x_target = t_col / 480;
    y_target = t_row / 640;
    z_target = t_val;

    x_match = b_cols(sortedIdx) ./ 480;
    y_match = b_rows(sortedIdx) ./ 640;
    z_match = b_vals(sortedIdx);

    % Filter out the match point from the background to prevent color overlap
    bg_mask = true(size(b_cols));
    bg_mask(sortedIdx) = false;
    x_bg = b_cols(bg_mask) ./ 480;
    y_bg = b_rows(bg_mask) ./ 640;
    z_bg = b_vals(bg_mask);

    % 2. Setup Figure
    fig = figure();

    % =========================================================================
    % 3. MAIN AXES SETUP & PLOTTING
    % =========================================================================
    ax = axes('Parent', fig);
    hold(ax, 'on');

    % Plot Data
    scatter3(ax, x_bg, y_bg, z_bg, 100, [0.7 0.7 0.7], '.');      % Background: Grey
    scatter3(ax, x_match, y_match, z_match, 400, 'r', '.');       % Match: Red
    scatter3(ax, x_target, y_target, z_target, 800, 'g', '.');    % Root: Green
    
    for i = 1:length(x_match)
        plot3(ax, [x_target, x_match(i)], ...
                  [y_target, y_match(i)], ...
                  [z_target, z_match(i)], ...
              '-k', 'LineWidth', 1.5);                            % Connecting Line
    end

    % Main Aesthetics
    xlabel(ax, 'X_{norm} [-]');
    ylabel(ax, 'Y_{norm} [-]');
    zlabel(ax, 't_{norm} [-]');
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    zlim(ax, [0 1]);
    
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    view(ax, 3);
    set(ax, 'FontSize', 16, 'FontName', 'Times New Roman');

    % =========================================================================
    % 4. INSET AXES SETUP (ZOOMED WINDOW)
    % =========================================================================
    % Create axes in the top-right corner [left bottom width height]
    axInset = axes('Parent', fig, 'Position', [0.65 0.65 0.25 0.25]);
    hold(axInset, 'on');

    % Re-plot data on the inset
    scatter3(axInset, x_bg, y_bg, z_bg, 100, [0.7 0.7 0.7], '.');
    scatter3(axInset, x_match, y_match, z_match, 400, 'r', '.');   
    scatter3(axInset, x_target, y_target, z_target, 800, 'g', '.'); 
    for i = 1:length(x_match)
        plot3(axInset, [x_target, x_match(i)], ...
                       [y_target, y_match(i)], ...
                       [z_target, z_match(i)], ...
              '-k', 'LineWidth', 1.5);
    end

    % Calculate bounds based strictly on the target and match coordinates
    margin = 0.01; % Adjust this for a tighter or looser zoom
    xlim(axInset, [min([x_target, min(x_match)]) - margin, max([x_target, max(x_match)]) + margin]);
    ylim(axInset, [min([y_target, min(y_match)]) - margin, max([y_target, max(y_match)]) + margin]);
    zlim(axInset, [min([z_target, min(z_match)]) - margin, max([z_target, max(z_match)]) + margin]);

    % Inset Aesthetics: Clean appearance without tick labels
    grid(axInset, 'on');
    box(axInset, 'on');
    axis(axInset,'square');
    view(axInset, 3);
    title(axInset, 'Zoomed View', 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
    set(axInset, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
    
    % =========================================================================
    % 5. FIGURE FINALIZATION & EXPORT
    % =========================================================================
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    
    journal.exportTight3DScatterPlots(gcf, '/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/Temporal-Persistance-Sample-Nom-Rot-NOM.pdf')
end

%  TO CALL FUNCTION
% target_idx = 1000;
% linearIndexA = validIndicesA(target_idx);
% % Extract the coordinates and value of the chosen target point
% [t_row, t_col] = ind2sub(imgSz, linearIndexA);
% t_val = mapA(linearIndexA);
% % --- 5. Find the mapping index (sortedIdx) ---
% % Extract the matching coordinates found by the vectorized function
% match_row = closestRows(target_idx);
% match_col = closestCols(target_idx);
% % Find where this specific match lives inside the b_rows/b_cols arrays
% sortedIdx = find(b_rows == match_row & b_cols == match_col, 1);
% % --- 6. Call the Visualization Function ---
% journal.showPersistance(t_row, t_col, t_val, ...
%     b_rows, b_cols, b_vals, sortedIdx);
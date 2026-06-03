function [] = showScatterPlotOfHotPixelAccumulatorMap(map)
% MAPTOSCATTERPLOT  3D scatter plot from a 2D value map with threshold bounds.
%
%   MAPTOSCATTERPLOT(MAP) converts a 2D map to a point
%   cloud and displays it as a 3D scatter plot. Points within the hot 
%   pixel threshold are colored red, and nominal points are colored green, 
%   matching their respective bounding volumes.
%
%   See also: process.generateMeshFromFrame, plot.mapToSurfPlot

    % --- Threshold Configuration ---
    % Adjust these values to match the specific accumulator thresholds 
    % defining your hot-pixel region.
    z_hot_lower = 0.2; % Lower bound of the hot pixel region
    z_hot_upper = 0.8; % Upper bound of the hot pixel region
    
    [pointCloud] = process.generateMeshFromFrame(map');
    x_trimmed = pointCloud(~isnan(pointCloud(:,3)), 1)./640;
    y_trimmed = pointCloud(~isnan(pointCloud(:,3)), 2)./480;
    z_trimmed = pointCloud(~isnan(pointCloud(:,3)), 3);
    
    % Dynamically cap the top and bottom of the green regions to the data range
    z_min_data = min(z_trimmed);
    z_max_data = max(z_trimmed);

    fig = figure();
    ax  = axes('Parent', fig);
    
    % Ensure OpenGL renderer is used for proper transparency export to PDF
    set(fig, 'Renderer', 'opengl'); 
    
    hold(ax, 'on'); % Retain the axes to overlay patches and scatter points

    % --- Draw Bounding Volumes ---
    x_bounds = [0, 1]; % Normalized X limits
    y_bounds = [0, 1]; % Normalized Y limits
    
    % 1. Bottom Nominal Region (Faint Green)
    if z_min_data < z_hot_lower
        drawFaintCube(ax, x_bounds, y_bounds, [z_min_data, z_hot_lower], 'g', 0.03);
    end
    
    % 2. Hot Pixel Region (Faint Red)
    drawFaintCube(ax, x_bounds, y_bounds, [z_hot_lower, z_hot_upper], 'r', 0.05);
    
    % 3. Top Nominal Region (Faint Green) - Optional, depends on your filter logic
    if z_max_data > z_hot_upper
        drawFaintCube(ax, x_bounds, y_bounds, [z_hot_upper, max(z_max_data, z_hot_upper+0.1)], 'g', 0.03);
    end

    % --- Point Cloud Coloring Logic ---
    % Initialize an N x 3 matrix for RGB colors, defaulting to a nominal green
    pointColors = repmat([0.15 0.7 0.15], length(z_trimmed), 1); 
    
    % Create a logical mask for points that fall strictly within the hot pixel bounds
    isHotPixel = (z_trimmed >= z_hot_lower) & (z_trimmed <= z_hot_upper);
    
    % Overwrite the color for hot pixels to a distinct red
    pointColors(isHotPixel, :) = repmat([0.85 0.15 0.15], sum(isHotPixel), 1);

    % --- Plot Scatter Data ---
    % Plotting this after the patches ensures the points are rendered on top/inside.
    % We pass the Nx3 pointColors matrix as the color argument.
    scatter3(ax, x_trimmed, y_trimmed, z_trimmed, 100, pointColors, '.');
    
    % --- Formatting ---
    view(ax, 3); % Forces the 3D perspective to prevent locking to 2D
    
    xlabel(ax, 'X_{norm} [-]')
    ylabel(ax, 'Y_{norm} [-]')
    zlabel(ax, '\rho_{norm} [-]')
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    
    % Adjust axis limits tightly to the normalized window
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    
    % Note: colormap(ax, jet) and colorbar(ax) have been removed 
    % since we are explicitly defining discrete point colors.
    
    set(ax, 'FontSize', 16, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    
    hold(ax, 'off');
    
    % Export to PDF
    exportgraphics(fig, '/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/Hot-Pixel-Accumulator-Map-Nominal-Rot.pdf');
end

% --- Local Helper Function ---
function drawFaintCube(ax, x_bounds, y_bounds, z_bounds, colorStr, alphaVal)
    % Generates an 8-vertex 3D box and draws it using patch for transparency support.
    
    % Define the 8 vertices of a rectangular prism
    x = [x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1) x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1)];
    y = [y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2) y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2)];
    z = [z_bounds(1) z_bounds(1) z_bounds(1) z_bounds(1) z_bounds(2) z_bounds(2) z_bounds(2) z_bounds(2)];
    
    % Define the 6 faces connecting the vertices
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    
    % Draw the patch onto the specified axes
    patch(ax, 'Vertices', [x' y' z'], 'Faces', faces, ...
        'FaceColor', colorStr, 'FaceAlpha', alphaVal, ...
        'EdgeColor', 'k', 'EdgeAlpha', 1.0, ...  % Solid black edges
        'LineWidth', 1.5, ...                    % Increased line weight for visibility
        'LineStyle', '-');
end
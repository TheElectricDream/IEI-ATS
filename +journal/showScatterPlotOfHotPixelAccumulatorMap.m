function [] = showScatterPlotOfHotPixelAccumulatorMap(map, name, show)
% MAPTOSCATTERPLOT  3D scatter plot from a 2D value map with threshold bounds.
%
%   MAPTOSCATTERPLOT(MAP) converts a 2D map to a point
%   cloud and displays it as a 3D scatter plot. Points above the hot 
%   pixel upper threshold are colored green, and nominal points are colored red, 
%   matching their respective bounding volumes.
%
%   See also: process.generateMeshFromFrame, plot.mapToSurfPlot

    % --- Threshold Configuration ---
    % Only an upper bound is used to define the hot-pixel region.
    z_hot_upper = 0.8; % Upper bound of the hot pixel region
    
    [pointCloud] = process.generateMeshFromFrame(map');
    x_trimmed = pointCloud(~isnan(pointCloud(:,3)), 1)./640;
    y_trimmed = pointCloud(~isnan(pointCloud(:,3)), 2)./480;
    z_trimmed = pointCloud(~isnan(pointCloud(:,3)), 3);
    
    % Dynamically cap the top of the green region to the data range
    z_min_data = min(z_trimmed);
    z_max_data = max(z_trimmed);

    if show
        fig = figure();
    else
        fig = figure('Visible', 'off'); % Create a figure without displaying it
    end
    ax  = axes('Parent', fig);
    
    % Ensure OpenGL renderer is used for proper transparency export to PDF
    set(fig, 'Renderer', 'opengl'); 
    
    hold(ax, 'on'); % Retain the axes to overlay patches and scatter points

    % --- Draw Bounding Volumes ---
    x_bounds = [0, 1]; % Normalized X limits
    y_bounds = [0, 1]; % Normalized Y limits
    
    % Prepare handles for legend entries
    hBelow = []; 
    hAbove = [];
    
    % 1. Nominal Region (Faint Red) - everything below the upper bound
    if z_min_data < z_hot_upper
        hBelow = drawFaintCube(ax, x_bounds, y_bounds, [z_min_data, min(z_hot_upper, z_max_data)], 'r', 0.05);
    end
    
    % 2. Hot Pixel Region (Faint Green) - above the upper bound up to data max
    if z_max_data > z_hot_upper
        hAbove = drawFaintCube(ax, x_bounds, y_bounds, [z_hot_upper, z_max_data], 'g', 0.03);
    end

    % --- Point Cloud Coloring Logic ---
    % Initialize an N x 3 matrix for RGB colors, defaulting to a nominal red
    pointColors = repmat([0.85 0.15 0.15], length(z_trimmed), 1); 
    
    % Create a logical mask for points that are above the hot pixel upper bound
    isHotPixel = (z_trimmed > z_hot_upper);
    
    % Overwrite the color for hot pixels to a distinct green
    pointColors(isHotPixel, :) = repmat([0.15 0.7 0.15], sum(isHotPixel), 1);

    % --- Plot Scatter Data ---
    scatter3(ax, x_trimmed, y_trimmed, z_trimmed, 100, pointColors, '.');
    
    % --- Formatting ---
    view(ax, 3); % Forces the 3D perspective to prevent locking to 2D
    
    xlabel(ax, 'X_{norm} [-]')
    ylabel(ax, 'Y_{norm} [-]')
    zlabel(ax, 'S(x,y) [-]')
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    
    % Adjust axis limits tightly to the normalized window
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    
    set(ax, 'FontSize', 16, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    
    % Add legend for the transparent regions. Use the patch handles if they exist.
    legendEntries = {};
    legendHandles = [];
    if ~isempty(hAbove)
        legendHandles(end+1) = hAbove; %#ok<AGROW>
        legendEntries{end+1} = 'Above Threshold'; %#ok<AGROW>
    end
    if ~isempty(hBelow)
        legendHandles(end+1) = hBelow; %#ok<AGROW>
        legendEntries{end+1} = 'Below Threshold'; %#ok<AGROW>
    end
    if ~isempty(legendHandles)
        % Create legend without altering axis children order
        lh = legend(ax, legendHandles, legendEntries, 'Location', 'northeast');
        set(lh, 'Box', 'on', 'FontSize', 16);
    end

    hold(ax, 'off');
    
    % Export to PDF
    journal.exportTight3DScatterPlots(gcf, ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name]);
end

% --- Local Helper Function ---
function h = drawFaintCube(ax, x_bounds, y_bounds, z_bounds, colorStr, alphaVal)
    % Generates an 8-vertex 3D box and draws it using patch for transparency support.
    
    % Define the 8 vertices of a rectangular prism
    x = [x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1) x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1)];
    y = [y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2) y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2)];
    z = [z_bounds(1) z_bounds(1) z_bounds(1) z_bounds(1) z_bounds(2) z_bounds(2) z_bounds(2) z_bounds(2)];
    
    % Define the 6 faces connecting the vertices
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    
    % Draw the patch onto the specified axes
    h = patch(ax, 'Vertices', [x' y' z'], 'Faces', faces, ...
        'FaceColor', colorStr, 'FaceAlpha', alphaVal, ...
        'EdgeColor', 'k', 'EdgeAlpha', 1.0, ...  % Solid black edges
        'LineWidth', 1.5, ...                    % Increased line weight for visibility
        'LineStyle', '-');
end
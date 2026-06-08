function [] = showScatterPlotOfLeakyBucketMap(map)
% SHOWSCATTERPLOTOFLEAKYBUCKETMAP  3D scatter plot of the leaky-bucket accumulator.
%
%   SHOWSCATTERPLOTOFLEAKYBUCKETMAP(MAP) converts a 2D accumulator map to a point
%   cloud and displays it as a 3D scatter plot. Points exceeding the persistent 
%   hot-pixel threshold (> 2.0) are colored red, and nominal points are colored green.
%
%   See also: process.generateMeshFromFrame, plot.mapToSurfPlot

    % --- Threshold Configuration ---
    % Defective pixels are flagged when the accumulator exceeds 2.0
    z_threshold = 2.0; 
    
    [pointCloud] = process.generateMeshFromFrame(map');
    x_trimmed = pointCloud(~isnan(pointCloud(:,3)), 1)./640;
    y_trimmed = pointCloud(~isnan(pointCloud(:,3)), 2)./480;
    z_trimmed = pointCloud(~isnan(pointCloud(:,3)), 3);
    
    % Dynamically cap the top of the red region to the data range, 
    % ensuring it renders cleanly even if the max value is very high.
    z_max_data = max(max(z_trimmed), z_threshold + 1.0); 

    fig = figure();
    ax  = axes('Parent', fig);
    
    % Ensure OpenGL renderer is used for proper transparency export to PDF
    set(fig, 'Renderer', 'opengl'); 
    
    hold(ax, 'on'); % Retain the axes to overlay patches and scatter points

    % --- Draw Bounding Volumes ---
    x_bounds = [0, 1]; % Normalized X limits
    y_bounds = [0, 1]; % Normalized Y limits
    
    % 1. Bottom Nominal Region (Faint Green, 0 to 2.0)
    drawFaintCube(ax, x_bounds, y_bounds, [0, z_threshold], 'g', 0.03);
    
    % 2. Hot Pixel Region (Faint Red, > 2.0)
    drawFaintCube(ax, x_bounds, y_bounds, [z_threshold, z_max_data], 'r', 0.05);

    % --- Point Cloud Coloring Logic ---
    % Initialize an N x 3 matrix for RGB colors, defaulting to a nominal green
    pointColors = repmat([0.15 0.7 0.15], length(z_trimmed), 1); 
    
    % Create a logical mask for points that exceed the threshold
    isHotPixel = (z_trimmed > z_threshold);
    
    % Overwrite the color for hot pixels to a distinct red
    pointColors(isHotPixel, :) = repmat([0.85 0.15 0.15], sum(isHotPixel), 1);

    % --- Plot Scatter Data ---
    scatter3(ax, x_trimmed, y_trimmed, z_trimmed, 100, pointColors, '.');
    
    % --- Formatting ---
    view(ax, 3); % Forces the 3D perspective to prevent locking to 2D
    
    xlabel(ax, 'X_{norm} [-]')
    ylabel(ax, 'Y_{norm} [-]')
    zlabel(ax, 'A_k [-]') % Updated to denote the accumulator value
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    
    % Adjust axis limits tightly to the normalized window
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    zlim(ax, [0 z_max_data]); % Clamp Z to our defined maximum
    
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    
    set(ax, 'FontSize', 16, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    
    hold(ax, 'off');
    
    % Export to PDF
    exportgraphics(fig, '/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/Leaky-Bucket-Accumulator-Map-Nominal-Rot.pdf');
end

% --- Local Helper Function ---
function drawFaintCube(ax, x_bounds, y_bounds, z_bounds, colorStr, alphaVal)
    % Generates an 8-vertex 3D box and draws it using patch for transparency support.
    
    x = [x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1) x_bounds(1) x_bounds(2) x_bounds(2) x_bounds(1)];
    y = [y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2) y_bounds(1) y_bounds(1) y_bounds(2) y_bounds(2)];
    z = [z_bounds(1) z_bounds(1) z_bounds(1) z_bounds(1) z_bounds(2) z_bounds(2) z_bounds(2) z_bounds(2)];
    
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    
    patch(ax, 'Vertices', [x' y' z'], 'Faces', faces, ...
        'FaceColor', colorStr, 'FaceAlpha', alphaVal, ...
        'EdgeColor', 'k', 'EdgeAlpha', 1.0, ...  
        'LineWidth', 1.5, ...                    
        'LineStyle', '-');
end
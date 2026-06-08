function [] = show2DScatterPlotOfHotPixelAccumulatorMap(map, name, show)
% MAPTO2DSCATTERPLOT  2D scatter plot from a 2D value map with threshold colors.
%
%   Converts a 2D map to a point cloud and displays it as a 2D scatter plot. 
%   Points above the hot pixel upper threshold are colored green, and 
%   nominal points are colored red.
%
%   See also: process.generateMeshFromFrame

    % --- Threshold Configuration ---
    z_hot_upper = 0.9; % Upper bound of the hot pixel region
    
    [pointCloud] = process.generateMeshFromFrame(map');
    
    % Extract and normalize coordinates
    x_trimmed = pointCloud(~isnan(pointCloud(:,3)), 1)./640;
    z_trimmed = pointCloud(~isnan(pointCloud(:,3)), 3);
    
    % Mirror Y-axis so the image orientation is correct
    y_trimmed = 1 - (pointCloud(~isnan(pointCloud(:,3)), 2)./480); 
    
    if show
        fig = figure();
    else
        fig = figure('Visible', 'off'); 
    end
    ax = axes('Parent', fig);
    
    hold(ax, 'on'); 
    
    % --- Point Cloud Grouping & Plotting ---
    % Create logical masks for point categorization
    isHotPixel = (z_trimmed > z_hot_upper);
    isNominal = ~isHotPixel;
    
    hNominal = []; 
    hHot = [];
    
    % Plot nominal points (Red)
    if any(isNominal)
        hNominal = scatter(ax, x_trimmed(isNominal), y_trimmed(isNominal), 500, [0.85 0.15 0.15], '.');
    end
    
    % Plot hot pixels (Green)
    if any(isHotPixel)
        hHot = scatter(ax, x_trimmed(isHotPixel), y_trimmed(isHotPixel), 500, [0.15 0.7 0.15], '.');
    end
    
    % --- Formatting ---
    xlabel(ax, 'X_{norm} [-]')
    ylabel(ax, 'Y_{norm} [-]')
    grid(ax, 'on');
    box(ax, 'on');
    
    % Adjust axis limits tightly to the normalized window
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    
    % Enforce the 640x480 (4:3) aspect ratio on the 2D plot box
    pbaspect(ax, [4 3 1]);
    
    set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    
    % --- Add Legend ---
    legendEntries = {};
    legendHandles = [];
    
    % Only add to legend if the points actually exist in the data
    if ~isempty(hHot)
        legendHandles(end+1) = hHot;
        legendEntries{end+1} = 'Above Threshold';
    end
    
    if ~isempty(hNominal)
        legendHandles(end+1) = hNominal;
        legendEntries{end+1} = 'Below Threshold';
    end
    
    if ~isempty(legendHandles)
        lh = legend(ax, legendHandles, legendEntries, 'Location', 'northeast');
        set(lh, 'Box', 'on', 'FontSize', 24);
    end
    
    hold(ax, 'off');
    
    % Export to PDF
    % (Keeping your original export call exactly as requested)
    journal.exportTight3DScatterPlots(gcf, ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name]);
end
function [] = show2DScatterOfLeakyBucketMap(sorted_x, sorted_y, hot_pixel_mask, name, show)
% SHOW2DSCATTEROFLEAKYBUCKETMAP  2D spatial scatter plot of the leaky-bucket accumulator.
%
%   SHOW2DSCATTEROFLEAKYBUCKETMAP(MAP) converts a 2D accumulator map to a 2D 
%   spatial point cloud. Points exceeding the persistent hot-pixel threshold (> 2.0) 
%   are enlarged and colored red. Nominal points are plotted faintly to show the 
%   structure of the scene without cluttering the defects.

    imgSz                       = [640, 480]; 

    filtered_x = sorted_x;
    filtered_y = sorted_y;

    % Strip the defective pixels globally
    hot_pixel_idx = find(hot_pixel_mask);

    if ~isempty(hot_pixel_idx)
        % Strip from the coordinate arrays
        sorted_lin_idx = sub2ind(imgSz, filtered_x, filtered_y);
        remove_mask = ismember(sorted_lin_idx, hot_pixel_idx);

        filtered_x(~remove_mask) = [];
        filtered_y(~remove_mask) = [];

    end

    
    if show
        fig = figure();
    else
        fig = figure('Visible', 'off'); 
    end
    
    ax  = axes('Parent', fig);
    set(fig, 'Renderer', 'opengl'); 
    hold(ax, 'on');
    
    % % --- Plot Scatter Data ---
    % % Plot nominal points first (in the background, small and faint gray)
    % scatter(ax, sorted_x./640, sorted_y./480, 5, [0.7 0.7 0.7], '.');
    % 
    % % Plot hot pixels on top (large and solid red)
    % scatter(ax, filtered_x./640, filtered_y./480, 40, [0.85 0.15 0.15], 'filled', 'MarkerEdgeColor', 'k');

    % Plot nominal points (Red)
    scatter(ax, sorted_x./640, sorted_y./480, 100, [0.7 0.7 0.7], '.');



    scatter(ax, filtered_x./640, filtered_y./480, 250, [0.15 0.7 0.15], '.');

    
    % --- Formatting ---
    xlabel(ax, 'X_{norm} [-]')
    ylabel(ax, 'Y_{norm} [-]')
    pbaspect(ax, [640 480 1]);
    grid(ax, 'on');
    box(ax, 'on');
    
    % Reverse Y-axis so it matches the top-left origin of an image/camera frame
    set(ax, 'YDir', 'reverse');
    
    % Adjust axis limits tightly to the normalized window
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    
    set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    
    hold(ax, 'off');
    
    % Export to PDF
    journal.exportTight3DScatterPlots(gcf, ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name]);
end
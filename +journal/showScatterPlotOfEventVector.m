function [] = showScatterPlotOfEventVector(x, y, z, name, show)
% VECTORSTOSCATTERPLOT  3D scatter plot from raw coordinate vectors.
%
%   VECTORSTOSCATTERPLOT(X, Y, Z, HOLDFIG) displays a 3D scatter
%   plot colored by z-value.
%
%   Inputs:
%     x, y, z - [N x 1] Coordinate vectors.
%     holdFig - Logical. If true, overlay on existing figure 16.
%
%   See also: plot.mapToScatterPlot

    x_trimmed = x./640;
    y_trimmed = y./480;
    % Normalize z to range [0,1]
    z_min = min(z(:));
    z_max = max(z(:));
    if z_max > z_min
        z_trimmed = (z - z_min) ./ (z_max - z_min);
    else
        z_trimmed = zeros(size(z)); % constant input -> map to zero
    end

    if show
        fig = figure();
    else
        fig = figure('Visible', 'off'); % Create a new figure if not holding
    end
    ax  = axes('Parent', fig);
    scatter3(x_trimmed, y_trimmed, z_trimmed, ...
        100, z_trimmed, '.');
    xlabel('X_{norm} [-]')
    ylabel('Y_{norm} [-]')
    zlabel('t_{norm} [-]')
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    zlim(ax, [0 1]);
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    colormap(jet);
    % colorbar;
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
    set(gcf, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    journal.exportTight3DScatterPlots(gcf, ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name]);

end
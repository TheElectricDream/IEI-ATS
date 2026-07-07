function [] = showScatterPlotOfRuleMaps3D(map, name, show)
% MAPTOSCATTERPLOT  3D scatter plot from a 2D value map.
%
%   MAPTOSCATTERPLOT(MAP, HOLDFIG) converts a 2D map to a point
%   cloud and displays it as a 3D scatter plot colored by z-value.
%
%   Inputs:
%     map     - [H x W] 2D value map.
%     holdFig - Logical. If true, overlay on existing figure 15.
%
%   See also: process.generateMeshFromFrame, plot.mapToSurfPlot

    [pointCloud] = process.generateMeshFromFrame(map');

    x_trimmed = pointCloud(~isnan(pointCloud(:,3)), 1)./640;
    y_trimmed = pointCloud(~isnan(pointCloud(:,3)), 2)./480;
    z_trimmed = pointCloud(~isnan(pointCloud(:,3)), 3);

    if show

        fig = figure();

    else

        fig = figure('Visible', 'off');
        
    end
    
    ax  = axes('Parent', fig);
    scatter3(x_trimmed, y_trimmed, z_trimmed, ...
        50, z_trimmed, '.');
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
    view(ax, 3);
    colormap(jet);
    %colorbar;
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
    set(gcf, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');

    journal.exportTight3DScatterPlots(gcf, ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name]);

end
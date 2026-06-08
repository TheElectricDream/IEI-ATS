function [] = showMetricLinePlot(metricData, name, show)
% SHOWMETRICLINEPLOT  2D line and scatter plot for sequential frame metrics.
%
%   SHOWMETRICLINEPLOT(METRICDATA, NAME, SHOW) plots a 1D array of metric
%   values (e.g., HotPixelCount) against their frame indices.
%
%   Inputs:
%     metricData - [1 x N] or [N x 1] array of metric values.
%     name       - String. Filename for the exported figure.
%     show       - Logical. If true, displays the figure; otherwise hidden.

% 1. Format the data 
y_data = metricData(:); % Ensure it is a column vector
x_data = (1:length(y_data))'; % X-axis corresponds to frame index

% 2. Initialize figure visibility
if show
    fig = figure();
else
    fig = figure('Visible', 'off');
end

ax = axes('Parent', fig);
hold(ax, 'on');

% 3. Plot the line connecting the points, then the scatter points on top
plot(ax, x_data, y_data, '-k', 'LineWidth', 1.5);
scatter(ax, x_data, y_data, 40, 'k', 'filled'); 

% 4. Format axes and labels (No title, per requirements)
xlabel(ax, 'Frame Index [-]');
ylabel(ax, 'Metric Count [-]'); 

grid(ax, 'on');
box(ax, 'on');

% 5. Set axis limits tightly to the data
if length(x_data) > 1
    xlim(ax, [min(x_data), max(x_data)]);
end

% Add a 5% padding to the Y-axis so the top/bottom scatter points don't clip the box
y_range = max(y_data) - min(y_data);
if y_range > 0
    ylim(ax, [min(y_data) - (0.05 * y_range), max(y_data) + (0.05 * y_range)]);
end

% 6. Apply publication typography
set(ax, 'FontSize', 16, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');

% 7. Export via the tight cropping routine to the manuscript directory
savePath = ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name];

journal.exportTight3DScatterPlots(fig, savePath);

hold(ax, 'off');
end
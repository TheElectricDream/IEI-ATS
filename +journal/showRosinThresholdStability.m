function [] = showRosinThresholdStability(theta_series, srr_series, name, show)
% SHOWROSINTHRESHOLDSTABILITY  Line plot for threshold and SRR over time.
%
%   SHOWROSINTHRESHOLDSTABILITY(THETA_SERIES, SRR_SERIES, NAME, SHOW) 
%
%   Inputs:
%     theta_series - [1 x N] array of threshold values per frame.
%     srr_series   - [1 x N] array of Signal Retention Rates per frame.
%     name         - String. Filename for the exported figure.
%     show         - Logical. If true, displays the figure; otherwise hidden.

% 1. Format the data 
y_theta = theta_series(:); 
y_srr   = srr_series(:);
x_data  = (1:length(y_theta))'; 

% 2. Initialize figure visibility
if show
    fig = figure();
else
    fig = figure('Visible', 'off');
end
ax = axes('Parent', fig);
hold(ax, 'on');

% 3. Plot the line connecting the points
%yyaxis(ax, 'left');
plot(ax, x_data, y_theta, '-k', 'LineWidth', 2);
ylabel(ax, 'Threshold [-]'); 
ax.YColor = 'k';

% yyaxis(ax, 'right');
% plot(ax, x_data, y_srr, '-k', 'LineWidth', 1.5, 'LineStyle', '--');
% ylabel(ax, 'Signal Retention Rate [-]');
% ax.YColor = 'k';

% 4. Format axes and labels (No title, per requirements)
xlabel(ax, 'Frame Index [-]');
grid(ax, 'on');
box(ax, 'on');

% 5. Set axis limits tightly to the data
if length(x_data) > 1
    xlim(ax, [min(x_data), max(x_data)]);
end

% y_range_theta = max(y_theta) - min(y_theta);
% if y_range_theta > 0
%     ylim(ax, [min(y_theta) - (0.05 * y_range_theta), max(y_theta) + (0.05 * y_range_theta)]);
% end

% yyaxis(ax, 'right');
% y_range_srr = max(y_srr) - min(y_srr);
% if y_range_srr > 0
%     ylim(ax, [min(y_srr) - (0.05 * y_range_srr), max(y_srr) + (0.05 * y_range_srr)]);
% end


% 6. Apply publication typography
set(ax, 'FontSize', 22, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');

pbaspect(ax, [640 640 1]);

% 2. Lock the axes tightly to the data box, ignoring outer labels during resize
% set(ax, 'ActivePositionProperty', 'position');

% 7. Export via the tight cropping routine to the manuscript directory
savePath = ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name];
journal.exportTight3DScatterPlots(fig, savePath);
hold(ax, 'off');
end
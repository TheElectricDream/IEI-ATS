function [] = showRosinThresholdConstructionNorm(diag, thresh_val, name, show)
% SHOWROSINTHRESHOLDCONSTRUCTIONNORM  Normalized geometric projection plot.
%
%   SHOWROSINTHRESHOLDCONSTRUCTIONNORM(DIAG, THRESH_VAL, NAME, SHOW) 
%
%   Inputs:
%     diag       - Struct output from testing.rosinThreshold.
%     thresh_val - Scalar threshold value.
%     name       - String. Filename for the exported figure.
%     show       - Logical. If true, displays the figure; otherwise hidden.

% 1. Format the data
idx_span = diag.peak_idx : diag.end_idx;
x_span = diag.bin_centers(idx_span);
y_span = diag.counts_smooth(idx_span);

x_n = (x_span - x_span(1)) / (x_span(end) - x_span(1));
y_n = (y_span - y_span(end)) / (y_span(1) - y_span(end));

% Find max deviation projection
t_rel_idx = find(idx_span == diag.thresh_idx, 1);
x0 = x_n(t_rel_idx);
y0 = y_n(t_rel_idx);
x_proj = (x0 - y0 + 1) / 2;
y_proj = (-x0 + y0 + 1) / 2;

% 2. Initialize figure visibility
if show
    fig = figure();
else
    fig = figure('Visible', 'off');
end
ax = axes('Parent', fig);
hold(ax, 'on');

% 3. Plot the curve, chord, and deviation
plot(ax, x_n, y_n, '-k', 'LineWidth', 2, 'DisplayName', 'Normalized Histogram');
plot(ax, [0 1], [1 0], '--k', 'LineWidth', 1.5, 'DisplayName', 'Chord');
plot(ax, [x0 x_proj], [y0 y_proj], '-r', 'LineWidth', 2, 'DisplayName', 'Max Deviation (d_b)');
scatter(ax, x0, y0, 60, 'r', 'filled', 'DisplayName', 'Optimal Threshold (\theta^*)');

% 4. Format axes and labels (No title, per requirements)
xlabel(ax, '$\tilde{x}_b$', 'Interpreter', 'latex');
ylabel(ax, '$\tilde{y}_b$', 'Interpreter', 'latex');
grid(ax, 'on');
box(ax, 'on');

% 5. Set axis limits tightly to the data
xlim(ax, [-0.05, 1.05]);
ylim(ax, [-0.05, 1.05]);
axis(ax, 'square'); % Maintain proper geometric aspect for the diagonal projection

pbaspect(ax, [640 640 1]);

% 6. Apply publication typography
set(ax, 'FontSize', 22, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
legend(ax, 'Location', 'northeast', 'FontSize', 18);

% 7. Export via the tight cropping routine to the manuscript directory
savePath = ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name];
journal.exportTight3DScatterPlots(fig, savePath);
hold(ax, 'off');
end
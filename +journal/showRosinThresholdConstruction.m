function [] = showRosinThresholdConstruction(diag, th, name, show)
% SHOWROSINTHRESHOLDCONSTRUCTION  Visualize Rosin unimodal-threshold selection.
%
%   Produces two figures from the diagnostics struct returned by
%   ROSINTHRESHOLD: (1) the (smoothed) coherence histogram with the
%   peak-to-tail chord and the selected threshold overlaid, and (2)
%   the per-bin perpendicular distance from the chord, whose maximum
%   defines the threshold.
%
%   Inputs:
%     diag  - Diagnostics struct from rosinThreshold (fields:
%             bin_centers, counts_smooth, chord_x, chord_y,
%             perp_dist, peak_idx, end_idx).
%     th    - Scalar selected threshold (bin-center value).
%     name  - Filename stem; '-Histogram.pdf' and '-PerpDist.pdf'
%             are appended.
%     show  - Logical. true displays figures; false renders offscreen.
%
%   See also: rosinThreshold, journal.exportTight3DScatterPlots

basePath = ['/home/alexandercrain/Dropbox/Graduate Documents/' ...
    'Doctor of Philosophy/Publications/Journals/AIAA Journal of ' ...
    'Spacecraft and Rockets/Event_Based_Spacecraft_Representation_' ...
    'Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/' ...
    'generated-figures/'];

c_chord = [0.85 0.10 0.10];   % chord / threshold accent
c_hist  = [0.80 0.80 0.80];   % histogram fill

% ----------------------------------------------------------------
% 1. Histogram with chord and selected threshold
% ----------------------------------------------------------------
if show, fig = figure(); else, fig = figure('Visible', 'off'); end
ax = axes('Parent', fig);
hold(ax, 'on');

bar(ax, diag.bin_centers, diag.counts_smooth, 1.0, ...
    'FaceColor', c_hist, 'EdgeColor', 'none');
plot(ax, diag.chord_x, diag.chord_y, '-', ...
    'Color', c_chord, 'LineWidth', 1.5);
yl = ylim(ax);
plot(ax, [th th], yl, '--k', 'LineWidth', 1.5);
ylim(ax, yl);

xlabel(ax, 'Coherence Score [-]');
ylabel(ax, 'Count [-]');
legend(ax, {'Histogram', 'Chord', '\theta^{*}'}, ...
    'Location', 'northeast', 'Box', 'off');
grid(ax, 'on'); box(ax, 'on');
pbaspect(ax, [4 3 1]);
set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultAxesFontName', 'Times New Roman');
journal.exportTight3DScatterPlots(fig, [basePath name '-Histogram.pdf']);

% ----------------------------------------------------------------
% 2. Perpendicular distance over the chord span
% ----------------------------------------------------------------
if show, fig = figure(); else, fig = figure('Visible', 'off'); end
ax = axes('Parent', fig);
hold(ax, 'on');

span = diag.peak_idx:diag.end_idx;
xs   = diag.bin_centers(span);
ds   = diag.perp_dist(span);
plot(ax, xs, ds, '-k', 'LineWidth', 1.5);

[d_max, mi] = max(ds);
plot(ax, xs(mi), d_max, 'o', 'MarkerSize', 9, ...
    'MarkerEdgeColor', c_chord, 'MarkerFaceColor', c_chord);
yl = ylim(ax);
plot(ax, [th th], yl, '--', 'Color', c_chord, 'LineWidth', 1.5);
ylim(ax, yl);

xlabel(ax, 'Coherence Score [-]');
ylabel(ax, 'd_{b} [-]');
grid(ax, 'on'); box(ax, 'on');
pbaspect(ax, [4 3 1]);
set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultAxesFontName', 'Times New Roman');
journal.exportTight3DScatterPlots(fig, [basePath name '-PerpDist.pdf']);
end
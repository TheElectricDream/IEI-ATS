function [] = showRosinThresholdConstructionNorm(diagStruct, th, name, show)
% SHOWROSINTHRESHOLDCONSTRUCTIONNORM  Scale-free Rosin construction figure.
%
%   Plots the coherence histogram in the [0,1]x[0,1] chord-span
%   coordinates used internally by ROSINTHRESHOLD, so the dominant
%   low-score noise spike no longer swamps the figure. Shows the
%   normalized histogram envelope, the peak-to-tail chord, and the
%   detected corner with its perpendicular drop to the chord.
%
%   Inputs:
%     diagStruct - Diagnostics struct (2nd output of rosinThreshold).
%     th         - Scalar selected threshold (actual score value),
%                  annotated in the legend.
%     name       - Filename stem; '-Norm.pdf' is appended.
%     show       - Logical. true displays; false renders offscreen.
%
%   See also: rosinThreshold, journal.showRosinThresholdConstruction

req = {'bin_centers', 'counts_smooth', 'perp_dist', ...
    'peak_idx', 'end_idx'};
if ~isstruct(diagStruct) || ~all(isfield(diagStruct, req))
    error('showRosinThresholdConstructionNorm:badInput', ...
        ['First arg must be the diagnostics struct (2nd output ' ...
        'of rosinThreshold); got a %s.'], class(diagStruct));
end

basePath = ['/home/alexandercrain/Dropbox/Graduate Documents/' ...
    'Doctor of Philosophy/Publications/Journals/AIAA Journal of ' ...
    'Spacecraft and Rockets/Event_Based_Spacecraft_Representation_' ...
    'Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/' ...
    'generated-figures/'];

c_chord  = [0.85 0.10 0.10];
c_curve  = [0.30 0.45 0.80];
c_corner = [0.10 0.10 0.10];

% ----------------------------------------------------------------
% 1. Reconstruct the algorithm's [0,1] chord-span normalization
% ----------------------------------------------------------------
span = diagStruct.peak_idx:diagStruct.end_idx;
xc   = diagStruct.bin_centers(span);
yc   = diagStruct.counts_smooth(span);

peak_count = diagStruct.counts_smooth(diagStruct.peak_idx);
x_range = xc(end) - xc(1);
y_range = peak_count - yc(end);

x_n = (xc - xc(1)) / x_range;     % peak -> 0,  tail anchor -> 1
y_n = (yc - yc(end)) / y_range;   % tail -> 0,  peak        -> 1

% Corner = max perpendicular distance (already computed in
% normalized space by rosinThreshold)
dperp = diagStruct.perp_dist(span);
[~, mi] = max(dperp);

% Foot of the perpendicular on the chord (for the drop line)
P0   = [x_n(1),  y_n(1)];                 % ~(0,1)
P1   = [x_n(end), y_n(end)];              % ~(1,0)
vhat = (P1 - P0) / norm(P1 - P0);
Pc   = [x_n(mi), y_n(mi)];
foot = P0 + dot(Pc - P0, vhat) * vhat;

% ----------------------------------------------------------------
% 2. Plot
% ----------------------------------------------------------------
if show, fig = figure(); else, fig = figure('Visible', 'off'); end
ax = axes('Parent', fig);
hold(ax, 'on');

area(ax, x_n, y_n, 'FaceColor', c_curve, 'FaceAlpha', 0.20, ...
    'EdgeColor', c_curve, 'LineWidth', 1.5);
plot(ax, [P0(1) P1(1)], [P0(2) P1(2)], '-', ...
    'Color', c_chord, 'LineWidth', 1.5);
plot(ax, [Pc(1) foot(1)], [Pc(2) foot(2)], ':', ...
    'Color', c_corner, 'LineWidth', 1.5);
plot(ax, Pc(1), Pc(2), 'o', 'MarkerSize', 9, ...
    'MarkerEdgeColor', c_corner, 'MarkerFaceColor', c_corner);

xlabel(ax, 'Normalized Score [-]');
ylabel(ax, 'Normalized Count [-]');
legend(ax, {'Histogram', 'Chord', 'd_{max}', ...
    sprintf('Corner (\\theta^{*} = %.4f)', th)}, ...
    'Location', 'northeast', 'Box', 'off');
xlim(ax, [0 1]); ylim(ax, [0 1]);
grid(ax, 'on'); box(ax, 'on');
pbaspect(ax, [4 3 1]);
set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultAxesFontName', 'Times New Roman');
journal.exportTight3DScatterPlots(fig, [basePath name '-Norm.pdf']);
end
function [] = showRosinSeparation(noise_scores, signal_scores, th, name, show)
% SHOWROSINSEPARATION  Validate Rosin threshold against known populations.
%
%   Overlays the noise and signal score histograms from a synthetic
%   map (where ground-truth membership is known) and marks the
%   selected threshold. The fraction of noise above and signal below
%   the threshold are reported as the false-positive and
%   false-negative rates.
%
%   Inputs:
%     noise_scores  - Vector of scores for ground-truth noise pixels.
%     signal_scores - Vector of scores for ground-truth signal pixels.
%     th            - Scalar selected threshold.
%     name          - Filename stem; '-Separation.pdf' is appended.
%     show          - Logical. true displays; false renders offscreen.
%
%   See also: rosinThreshold, demo_rosinThreshold

basePath = ['/home/alexandercrain/Dropbox/Graduate Documents/' ...
    'Doctor of Philosophy/Publications/Journals/AIAA Journal of ' ...
    'Spacecraft and Rockets/Event_Based_Spacecraft_Representation_' ...
    'Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/' ...
    'generated-figures/'];

c_noise  = [0.30 0.45 0.80];
c_signal = [0.90 0.55 0.10];
c_thresh = [0.10 0.10 0.10];

% ----------------------------------------------------------------
% 1. Error rates at the selected threshold
% ----------------------------------------------------------------
fpr = 100 * sum(noise_scores  > th) / max(numel(noise_scores), 1);
fnr = 100 * sum(signal_scores < th) / max(numel(signal_scores), 1);

% ----------------------------------------------------------------
% 2. Overlaid histograms (shared binning) on a log count axis
% ----------------------------------------------------------------
if show, fig = figure(); else, fig = figure('Visible', 'off'); end
ax = axes('Parent', fig);
hold(ax, 'on');

lo = min([noise_scores(:); signal_scores(:)]);
hi = max([noise_scores(:); signal_scores(:)]);
edges = linspace(lo, hi, 80);

histogram(ax, noise_scores,  edges, 'FaceColor', c_noise,  ...
    'EdgeColor', 'none', 'FaceAlpha', 0.55);
histogram(ax, signal_scores, edges, 'FaceColor', c_signal, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.55);
set(ax, 'YScale', 'log');

yl = ylim(ax);
plot(ax, [th th], yl, '--', 'Color', c_thresh, 'LineWidth', 1.5);
ylim(ax, yl);

xlabel(ax, 'Coherence Score [-]');
ylabel(ax, 'Count [-]');
legend(ax, { sprintf('Noise (%.2f%% above)',  fpr), ...
    sprintf('Signal (%.2f%% below)', fnr), ...
    '\theta^{*}' }, ...
    'Location', 'northeast', 'Box', 'off');
grid(ax, 'on'); box(ax, 'on');
pbaspect(ax, [4 3 1]);
set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultAxesFontName', 'Times New Roman');
journal.exportTight3DScatterPlots(fig, [basePath name '-Separation.pdf']);
end
function [] = showRosinThresholdStability(theta_series, srr_series, name, show)
% SHOWROSINTHRESHOLDSTABILITY  Per-frame threshold and retention plots.
%
%   Plots the Rosin-selected coherence threshold and the signal
%   retention rate (SRR) as a function of frame index over a full
%   sequence, demonstrating that the per-frame selector is stable
%   rather than jittery.
%
%   Inputs:
%     theta_series - Vector of per-frame thresholds (FilterThreshold).
%     srr_series   - Vector of per-frame signal retention rates.
%     name         - Filename stem; '-Threshold.pdf' and '-SRR.pdf'
%                    are appended.
%     show         - Logical. true displays; false renders offscreen.
%
%   See also: computeSignalRetentionRate, rosinThreshold

basePath = ['/home/alexandercrain/Dropbox/Graduate Documents/' ...
    'Doctor of Philosophy/Publications/Journals/AIAA Journal of ' ...
    'Spacecraft and Rockets/Event_Based_Spacecraft_Representation_' ...
    'Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/' ...
    'generated-figures/'];

% ----------------------------------------------------------------
% 1. Selected threshold per frame
% ----------------------------------------------------------------
if show, fig = figure(); else, fig = figure('Visible', 'off'); end
ax = axes('Parent', fig);
hold(ax, 'on');
x_data = linspace(1, numel(theta_series), numel(theta_series));
plot(ax, x_data, theta_series, '-k', 'LineWidth', 1.5);
xlabel(ax, 'Frame Index [-]');
ylabel(ax, '\theta^{*} [-]');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0 numel(theta_series)]);
pbaspect(ax, [4 3 1]);
set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultAxesFontName', 'Times New Roman');
journal.exportTight3DScatterPlots(fig, [basePath name '-Threshold.pdf']);

% ----------------------------------------------------------------
% 2. Signal retention rate per frame
% ----------------------------------------------------------------
if show, fig = figure(); else, fig = figure('Visible', 'off'); end
ax = axes('Parent', fig);
hold(ax, 'on');
x_data = linspace(1, numel(srr_series), numel(srr_series));
plot(ax, x_data, srr_series, '-k', 'LineWidth', 1.5);
xlabel(ax, 'Frame Index [-]');
ylabel(ax, 'SRR [-]');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0 numel(srr_series)]);
ylim(ax, [0 1]);
pbaspect(ax, [4 3 1]);
set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', ...
    'DefaultAxesFontName', 'Times New Roman');
journal.exportTight3DScatterPlots(fig, [basePath name '-SRR.pdf']);
end
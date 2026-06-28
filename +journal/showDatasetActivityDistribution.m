function [] = showDatasetActivityDistribution(xk, yk, tk, name, show)
% SHOWDATASETACTIVITYDISTRIBUTION  Per-pixel event-count distribution.
%
%   Aggregates the entire event stream into a per-pixel event count
%   and visualizes its distribution to assess long-/heavy-tailedness.
%   Two panels are exported:
%     (1) log-binned histogram of events-per-pixel on log-log axes,
%     (2) complementary CDF P(X >= x) on log-log axes, the standard
%         heavy-tail diagnostic. A straight CCDF suggests a power-law
%         tail; downward curvature a lighter (e.g. lognormal) tail.
%
%   Only active pixels (those firing at least once) are included.
%
%   Inputs:
%     xk, yk - Event pixel coordinates over the whole dataset.
%     tk     - Event timestamps [s] (used to report duration).
%     name   - Filename stem; '-EventCountHist.pdf' and
%              '-EventCountCCDF.pdf' are appended.
%     show   - Logical. true displays; false renders offscreen.
%
%   See also: rosinThreshold, journal.showRosinThresholdConstruction

    basePath = ['/home/alexandercrain/Dropbox/Graduate Documents/' ...
        'Doctor of Philosophy/Publications/Journals/AIAA Journal of ' ...
        'Spacecraft and Rockets/Event_Based_Spacecraft_Representation_' ...
        'Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/' ...
        'generated-figures/'];

    n_bins = 50;
    c_hist = [0.80 0.80 0.80];
    c_line = [0.10 0.10 0.10];

    % ----------------------------------------------------------------
    % 1. Per-pixel event count (active pixels only, no imgSz needed)
    % ----------------------------------------------------------------
    [~, ~, g] = unique([double(xk(:)), double(yk(:))], 'rows');
    counts = accumarray(g, 1);          % events per active pixel
    counts = counts(counts > 0);

    if numel(counts) < 10
        warning('showDatasetActivityDistribution:tooFewPixels', ...
            'Fewer than 10 active pixels; nothing to plot.');
        return;
    end

    % ----------------------------------------------------------------
    % 2. Quantitative tail indicators (printed for the manuscript)
    % ----------------------------------------------------------------
    T        = max(tk) - min(tk);
    N_pix    = numel(counts);
    sorted_c = sort(counts, 'descend');
    top1pct  = sorted_c(1:max(1, round(0.01 * N_pix)));
    conc     = 100 * sum(top1pct) / sum(counts);
    fprintf('Dataset duration:            %.1f s\n', T);
    fprintf('Active pixels:               %d\n', N_pix);
    fprintf('Events/pixel  median | max:  %d | %d\n', ...
        median(counts), max(counts));
    fprintf('Top 1%% of pixels hold:       %.1f%% of all events\n', conc);

    % ----------------------------------------------------------------
    % 3. Panel 1 — log-binned histogram, log-log axes
    % ----------------------------------------------------------------
    if show, fig = figure(); else, fig = figure('Visible', 'off'); end
    ax = axes('Parent', fig);
    hold(ax, 'on');

    edges = logspace(0, log10(double(max(counts))), n_bins + 1);
    histogram(ax, counts, edges, 'FaceColor', c_hist, 'EdgeColor', 'none');
    set(ax, 'XScale', 'log', 'YScale', 'log');

    xlabel(ax, 'Events per Pixel [-]');
    ylabel(ax, 'Number of Pixels [-]');
    grid(ax, 'on'); box(ax, 'on');
    pbaspect(ax, [4 3 1]);
    set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', ...
        'DefaultAxesFontName', 'Times New Roman');
    journal.exportTight3DScatterPlots(fig, [basePath name '-EventCountHist.pdf']);

    % ----------------------------------------------------------------
    % 4. Panel 2 — complementary CDF P(X >= x), log-log axes
    % ----------------------------------------------------------------
    if show, fig = figure(); else, fig = figure('Visible', 'off'); end
    ax = axes('Parent', fig);
    hold(ax, 'on');

    cs = sort(counts, 'ascend');
    N  = numel(cs);
    [xv, ifirst] = unique(cs, 'first');     % light: unique values only
    ccdf = (N - ifirst + 1) / N;            % P(X >= xv)
    plot(ax, xv, ccdf, '-', 'Color', c_line, 'LineWidth', 1.5);
    set(ax, 'XScale', 'log', 'YScale', 'log');

    xlabel(ax, 'Events per Pixel, x [-]');
    ylabel(ax, 'P(X \geq x) [-]');
    grid(ax, 'on'); box(ax, 'on');
    pbaspect(ax, [4 3 1]);
    set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', ...
        'DefaultAxesFontName', 'Times New Roman');
    journal.exportTight3DScatterPlots(fig, [basePath name '-EventCountCCDF.pdf']);
end
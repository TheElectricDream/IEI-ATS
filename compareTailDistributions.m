function [] = compareTailDistributions(user_data, user_label, name, show)
% COMPARETAILDISTRIBUTIONS  Reference plate for assessing heavy-tailedness.
%
%   Overlays a user-supplied empirical sample against reference
%   distributions of KNOWN tail class, in two views: a log-binned
%   histogram and a complementary CDF, both on log-log axes. The CCDF
%   is the decisive panel: a power law is a straight line, an
%   exponential drops in a concave cliff, and lognormal /
%   stretched-exponential curve gently between them.
%
%   All samples are normalized to unit median so their bodies align
%   and only the tails diverge (scaling x is harmless on a log-log
%   CCDF -- it shifts the curve without changing its slope).
%
%   Inputs:
%     user_data  - Empirical sample to assess (e.g. per-pixel event
%                  counts from showDatasetActivityDistribution). Pass
%                  [] to plot the references alone.
%     user_label - Legend label for the user data (e.g. 'EVOS-NOM').
%     name       - Filename stem; '-TailHist.pdf' and '-TailCCDF.pdf'
%                  are appended.
%     show       - Logical. true displays; false renders offscreen.
%
%   See also: showDatasetActivityDistribution, rosinThreshold
%
%   Ref: Clauset, Shalizi & Newman, "Power-law distributions in
%        empirical data," SIAM Review 51(4):661-703, 2009.
%        doi:10.1137/070710111

    basePath = ['/home/alexandercrain/Dropbox/Graduate Documents/' ...
        'Doctor of Philosophy/Publications/Journals/AIAA Journal of ' ...
        'Spacecraft and Rockets/Event_Based_Spacecraft_Representation_' ...
        'Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/' ...
        'generated-figures/'];

    rng(7);                       % reproducible references
    N = 200000;
    n_bins = 45;

    % ----------------------------------------------------------------
    % 1. Reference samples (inverse transform -> no Stats toolbox)
    % ----------------------------------------------------------------
    U = rand(N, 1);
    Z = randn(N, 1);

    refs(1).label = 'Exponential (Poisson null, light)';
    refs(1).data  = -log(U);
    refs(1).color = [0.17 0.50 0.72];
    refs(1).style = '-';

    refs(2).label = 'Stretched-exp (Weibull k=0.5)';
    refs(2).data  = (-log(U)).^2;
    refs(2).color = [0.25 0.67 0.36];
    refs(2).style = '-';

    refs(3).label = 'Lognormal (\sigma=1)';
    refs(3).data  = exp(1.0 * Z);
    refs(3).color = [0.85 0.37 0.05];
    refs(3).style = '-';

    refs(4).label = 'Power-law (\alpha=1.5, heavy)';
    refs(4).data  = U.^(-1/1.5);
    refs(4).color = [0.80 0.07 0.07];
    refs(4).style = '-';

    % Append user data as a thick black series, if provided
    if ~isempty(user_data)
        k = numel(refs) + 1;
        refs(k).label = user_label;
        refs(k).data  = double(user_data(:));
        refs(k).color = [0 0 0];
        refs(k).style = '-';
    end

    % Normalize every series to unit median
    for i = 1:numel(refs)
        d = refs(i).data;
        d = d(d > 0);
        refs(i).data = d / median(d);
    end

    % ----------------------------------------------------------------
    % 2. Panel 1 -- log-binned density as lines, log-log
    % ----------------------------------------------------------------
    if show, fig = figure(); else, fig = figure('Visible', 'off'); end
    ax = axes('Parent', fig);
    hold(ax, 'on');

    for i = 1:numel(refs)
        d = refs(i).data;
        edges = logspace(log10(min(d)), log10(max(d)), n_bins + 1);
        h = histcounts(d, edges, 'Normalization', 'pdf');
        c = sqrt(edges(1:end-1) .* edges(2:end));   % geometric centers
        m = h > 0;
        lw = 1.5 + 1.5 * (i == numel(refs) && ~isempty(user_data));
        plot(ax, c(m), h(m), refs(i).style, ...
            'Color', refs(i).color, 'LineWidth', lw);
    end
    set(ax, 'XScale', 'log', 'YScale', 'log');

    xlabel(ax, 'Value / median [-]');
    ylabel(ax, 'Probability Density [-]');
    legend(ax, {refs.label}, 'Location', 'southwest', 'Box', 'off', ...
        'FontSize', 14);
    grid(ax, 'on'); box(ax, 'on');
    pbaspect(ax, [4 3 1]);
    set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', ...
        'DefaultAxesFontName', 'Times New Roman');
    journal.exportTight3DScatterPlots(fig, [basePath name '-TailHist.pdf']);

    % ----------------------------------------------------------------
    % 3. Panel 2 -- CCDF P(X >= x), log-log (the decisive panel)
    % ----------------------------------------------------------------
    if show, fig = figure(); else, fig = figure('Visible', 'off'); end
    ax = axes('Parent', fig);
    hold(ax, 'on');

    for i = 1:numel(refs)
        s = sort(refs(i).data, 'ascend');
        n = numel(s);
        ccdf = (n - (1:n)' + 1) / n;     % P(X >= s)
        lw = 1.5 + 1.5 * (i == numel(refs) && ~isempty(user_data));
        plot(ax, s, ccdf, refs(i).style, ...
            'Color', refs(i).color, 'LineWidth', lw);
    end
    set(ax, 'XScale', 'log', 'YScale', 'log');
    ylim(ax, [1/N 1.5]);

    xlabel(ax, 'x / median [-]');
    ylabel(ax, 'P(X \geq x) [-]');
    legend(ax, {refs.label}, 'Location', 'southwest', 'Box', 'off', ...
        'FontSize', 14);
    grid(ax, 'on'); box(ax, 'on');
    pbaspect(ax, [4 3 1]);
    set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', ...
        'DefaultAxesFontName', 'Times New Roman');
    journal.exportTight3DScatterPlots(fig, [basePath name '-TailCCDF.pdf']);
end
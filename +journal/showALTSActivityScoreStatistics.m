function [] = showALTSActivityScoreStatistics(alts_activity_score, name, show)

    if show
        fig = figure();
    else
        fig = figure('Visible', 'off'); 
    end
    ax = axes('Parent', fig);
    hold(ax, 'on'); 
    x_data = linspace(1,length(alts_activity_score.mean),length(alts_activity_score.mean));
    plot(ax, x_data, alts_activity_score.mean, '-k', 'LineWidth', 1.5);
    % scatter(ax, x_data, alts_activity_score.mean, 40, 'k', 'filled'); 
    xlabel(ax, 'Frame Index [-]')
    ylabel(ax, '\mu_{ALTS} [-]')
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [0 650]);
    pbaspect(ax, [4 3 1]);
    set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    
    % Export to PDF
    journal.exportTight3DScatterPlots(gcf, ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name '-Mean.pdf']);

    if show
        fig = figure();
    else
        fig = figure('Visible', 'off'); 
    end
    ax = axes('Parent', fig);
    hold(ax, 'on'); 
    x_data = linspace(1,length(alts_activity_score.std),length(alts_activity_score.std));
    plot(ax, x_data, alts_activity_score.std, '-k', 'LineWidth', 1.5);
    % scatter(ax, x_data, alts_activity_score.std, 40, 'k', 'filled'); 
    xlabel(ax, 'Frame Index [-]')
    ylabel(ax, '\sigma_{ALTS} [-]')
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [0 650]);
    pbaspect(ax, [4 3 1]);
    set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');

    % Export to PDF
    journal.exportTight3DScatterPlots(gcf, ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name '-StandDev.pdf']);

    if show
        fig = figure();
    else
        fig = figure('Visible', 'off'); 
    end
    ax = axes('Parent', fig);
    hold(ax, 'on'); 
    x_data = linspace(1,length(alts_activity_score.median),length(alts_activity_score.median));
    plot(ax, x_data, alts_activity_score.median, '-k', 'LineWidth', 1.5);
    % scatter(ax, x_data, alts_activity_score.median, 40, 'k', 'filled'); 
    xlabel(ax, 'Frame Index [-]')
    ylabel(ax, 'Median_{ALTS} [-]')
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [0 650]);
    pbaspect(ax, [4 3 1]);
    set(ax, 'FontSize', 24, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');

    % Export to PDF
    journal.exportTight3DScatterPlots(gcf, ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name '-Median.pdf']);




end
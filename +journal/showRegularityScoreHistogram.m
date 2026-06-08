function [] = showRegularityScoreHistogram(regularity_map, name, show)
% SHOWREGULARITYSCOREHISTOGRAM  1D histogram of the regularity scores S(x,y).
%
%   SHOWREGULARITYSCOREHISTOGRAM(REGULARITY_MAP) plots the relative frequency
%   of regularity scores across the frame, highlighting the threshold that
%   separates structured spacecraft edges from aperiodic noise.

    % Flatten the map and ignore empty pixels (zeros or NaNs)
    scores = regularity_map(:);
    valid_scores = scores(scores > 0 & ~isnan(scores));
    if show
        fig = figure();
    else
        fig = figure('Visible', 'off'); % Create a figure without displaying it
    end
        
    ax  = axes('Parent', fig);
    set(fig, 'Renderer', 'opengl'); 
    hold(ax, 'on');

    % Plot the histogram
    histogram(ax, valid_scores, 100, 'Normalization', 'probability', ...
        'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

    % Draw the Aperiodic Threshold line
    threshold = 0.9;
    xline(ax, threshold, '--r', 'Aperiodic Threshold (S = 0.9)', ...
        'LineWidth', 2.5, 'LabelVerticalAlignment', 'top', ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 14, 'FontName', 'Times New Roman');

    % --- Formatting ---
    xlabel(ax, 'Regularity Score, S(x,y) [-]')
    ylabel(ax, 'Relative Frequency [-]')
    grid(ax, 'on');
    box(ax, 'on');
    axis(ax, 'square')
    
    xlim(ax, [0 1]);
    grid(ax, 'on');   
    
    set(ax, 'FontSize', 16, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    
    hold(ax, 'off');
    
    % Export to PDF
    journal.exportTight3DScatterPlots(gcf, ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name]);
end
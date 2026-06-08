function [] = showLeakyBucketEvolution(frames, hot_history, edge_history, noise_history)
% SHOWLEAKYBUCKETEVOLUTION  Line plot showing A_k progression over time.
%
%   Requires 1D arrays representing the frame indices, and the A_k value 
%   at each frame for three representative pixels: a real hot pixel, a 
%   genuine spacecraft edge, and random thermal noise.

    fig = figure();
    ax  = axes('Parent', fig);
    set(fig, 'Renderer', 'opengl'); 
    hold(ax, 'on');

    % Plot the three representative pixel histories
    p1 = plot(ax, frames, hot_history, '-r', 'LineWidth', 2.5);
    p2 = plot(ax, frames, edge_history, '-b', 'LineWidth', 2.0);
    p3 = plot(ax, frames, noise_history, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

    % Draw the Defect Threshold line
    threshold = 2.0;
    yline(ax, threshold, '--k', 'Defect Threshold (A_k = 2.0)', ...
        'LineWidth', 2.0, 'LabelHorizontalAlignment', 'left', ...
        'FontSize', 14, 'FontName', 'Times New Roman');

    % --- Formatting ---
    xlabel(ax, 'Frame Index, k [-]')
    ylabel(ax, 'Accumulator Value, A_k [-]')
    grid(ax, 'on');
    box(ax, 'on');
    
    % Add a legend
    legend(ax, [p1, p2, p3], {'Defective Pixel', 'Spacecraft Edge', 'Thermal Noise'}, ...
        'Location', 'northwest', 'FontSize', 14);
    
    % Dynamically set axis limits
    xlim(ax, [min(frames) max(frames)]);
    ylim(ax, [0 max(max(hot_history), threshold + 1.0)]);
    
    set(ax, 'FontSize', 16, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    
    hold(ax, 'off');
    
    % Export to PDF
    exportgraphics(fig, '/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/Leaky-Bucket-Time-Evolution.pdf');
end
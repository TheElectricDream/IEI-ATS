function [] = showSpatialMaskEvolution(raw_counts, aperiodic_mask, final_hot_mask)
% SHOWSPATIALMASKEVOLUTION  Generates three independent figures for the spatial mask.
%
%   Displays and exports (1) Raw Event Counts (Log-Normalized), 
%   (2) The Aperiodic Mask (S < 0.8), and (3) The final thresholded hot-pixel mask.
%   These are output as separate PDFs to be combined via LaTeX subfigures.

    % Custom colormaps for binary masks
    bw_cmap  = [1 1 1; 0 0 0];          % 0=White, 1=Black
    red_cmap = [1 1 1; 0.85 0.15 0.15]; % 0=White, 1=Red

    % Common Formatting Options
    fontName = 'Times New Roman';
    fontSize = 18;
    exportDir = '/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/';

    % =========================================================================
    % FIGURE 1: RAW EVENT COUNTS (Log-Normalized for Visibility)
    % =========================================================================
    fig1 = figure('Name', 'Raw Event Counts');
    ax1  = axes('Parent', fig1);
    set(fig1, 'Renderer', 'opengl', 'DefaultTextFontName', fontName, 'DefaultAxesFontName', fontName); 
    
    % Apply log1p to compress extreme hot-pixel outliers and reveal faint structure
    log_counts = log1p(raw_counts');
    
    imagesc(ax1, log_counts); 
    colormap(ax1, 'parula');
    axis(ax1, 'image'); 
    axis(ax1, 'off');
    
    journal.exportTight3DScatterPlots(fig1, [exportDir, 'Spatial-Mask-1-Raw-Counts.png']);
    
    % =========================================================================
    % FIGURE 2: APERIODIC MASK
    % =========================================================================
    fig2 = figure('Name', 'Aperiodic Mask');
    ax2  = axes('Parent', fig2);
    set(fig2, 'Renderer', 'opengl', 'DefaultTextFontName', fontName, 'DefaultAxesFontName', fontName); 
    
    imagesc(ax2, aperiodic_mask'); 
    colormap(ax2, bw_cmap);
    axis(ax2, 'image'); 
    axis(ax2, 'off');
    
    journal.exportTight3DScatterPlots(fig2, [exportDir, 'Spatial-Mask-2-Aperiodic.png']);

    % =========================================================================
    % FIGURE 3: FINAL HOT-PIXEL MASK
    % =========================================================================
    fig3 = figure('Name', 'Final Defect Mask');
    ax3  = axes('Parent', fig3);
    set(fig3, 'Renderer', 'opengl', 'DefaultTextFontName', fontName, 'DefaultAxesFontName', fontName); 
    
    imagesc(ax3, final_hot_mask'); 
    colormap(ax3, red_cmap);
    axis(ax3, 'image'); 
    axis(ax3, 'off');
    
    journal.exportTight3DScatterPlots(fig3, [exportDir, 'Spatial-Mask-3-Final-Defect.png']);

    % Clean up by closing the independent figures if running in a loop
    % close(fig1); close(fig2); close(fig3);
end

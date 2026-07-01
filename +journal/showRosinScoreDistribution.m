function [] = showRosinScoreDistribution(diag, thresh_val, name, show)
% SHOWROSINSCOREDISTRIBUTION  Log-scale bar plot of normalized active scores.
%
%   SHOWROSINSCOREDISTRIBUTION(DIAG, THRESH_VAL, NAME, SHOW) 

% 1. Format the data and normalize the Y-axis to [0, 1]
x_data = diag.bin_centers(:);
norm_factor = max(diag.counts(:)); 

y_data = diag.counts(:) / norm_factor;
y_smooth = diag.counts_smooth(:) / norm_factor;

% 2. Initialize figure visibility
if show
    fig = figure();
else
    fig = figure('Visible', 'off');
end
ax = axes('Parent', fig);
hold(ax, 'on');

% 3. Plot the histogram, smoothed curve, noise peak, and threshold
bar(ax, x_data, y_data, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none', 'DisplayName', 'Raw Counts');
plot(ax, x_data, y_smooth, '-k', 'LineWidth', 1.5, 'DisplayName', 'Smoothed Counts');
%xline(ax, diag.bin_centers(diag.peak_idx), '--b', 'LineWidth', 1.5, 'DisplayName', 'Noise Peak (b_p)');
xline(ax, thresh_val, '-r', 'LineWidth', 2, 'DisplayName', 'Threshold (\theta^*)');

% 4. Format axes and labels (No title, per requirements)
xlabel(ax, 'Normalized Coherence Score [-]');
ylabel(ax, 'Log10 of Normalized Pixel Count [-]'); 
grid(ax, 'on');
box(ax, 'on');

% 5. Apply Logarithmic Scale and handle limits
set(ax, 'YScale', 'log');

if length(x_data) > 1
    xlim(ax, [min(x_data), max(x_data)]);
end

% Find the smallest non-zero value to set a clean lower bound for the log axis
min_nonzero = min(y_data(y_data > 0));
if isempty(min_nonzero)
    min_nonzero = 10^-6; % Fallback 
end
% Set limit from slightly below the minimum non-zero value up to 1.5 (to leave room for the legend)
ylim(ax, [min_nonzero * 0.5, 2.0]); 

% 6. Apply publication typography & Aspect Ratio Lock
set(ax, 'FontSize', 26, 'FontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
legend(ax, 'Location', 'northeast', 'FontSize', 22);

% Enforce 4:3 aspect ratio and lock it against outer dimension changes
pbaspect(ax, [640 480 1]);
set(ax, 'ActivePositionProperty', 'position'); 

% 7. Export via the tight cropping routine to the manuscript directory
savePath = ['/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/' name];
journal.exportTight3DScatterPlots(fig, savePath);
hold(ax, 'off');
end
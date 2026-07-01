function showAdaptiveThresholdFigures(diag, thresh_val, theta_series, varargin)
% SHOWADAPTIVETHRESHOLDFIGURES Generates publication-ready figures for the 
% adaptive thresholding section using real pipeline data.
%
%   showAdaptiveThresholdFigures(diag, thresh_val, theta_series)
% 
% Inputs:
%   diag         : The 'd' diagnostics struct from testing.rosinThreshold
%                  for a representative frame (used for Figs Xa and Xb).
%   thresh_val   : The scalar threshold value for that specific frame.
%   theta_series : An array (1 x N_frames) containing the calculated 
%                  thresholds over time (used for Fig Y).
%
% Optional Inputs:
%   'ScatterData': A cell array of size (1 x N_frames) containing raw active 
%                  scores per frame. If provided, it plots the grey background 
%                  scatter in Figure Y. If omitted, plots only the trend line.

    % Parse optional scatter data
    p = inputParser;
    addOptional(p, 'ScatterData', {});
    parse(p, varargin{:});
    scatter_data = p.Results.ScatterData;

    %% --------------------------------------------------------------------
    % FIGURE Xa: The Justification (Raw Histogram of Real Data)
    % --------------------------------------------------------------------
    if ~isempty(diag)
        fig_xa = figure('Name', 'Figure Xa: Score Distribution', 'Position', [100, 100, 600, 400], 'Color', 'w');
        
        % Plot the raw counts from the diagnostics struct
        bar(diag.bin_centers, diag.counts, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
        hold on; grid on;

        % Highlight the peak and the chosen threshold
        xline(diag.bin_centers(diag.peak_idx), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Noise Peak ($b_p$)');
        xline(thresh_val, 'r-', 'LineWidth', 2, 'DisplayName', 'Threshold ($\theta^{\star}$)');

        xlabel('Normalized Coherence Score $\bar{C}(x,y)$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('Pixel Count $h(b)$', 'Interpreter', 'latex', 'FontSize', 12);
        title('Figure Xa: Unimodal Score Distribution', 'FontSize', 14);
        
        lgd = legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);
        lgd.ItemTokenSize = [20, 18]; % Clean up legend line lengths
        set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
    end

    %% --------------------------------------------------------------------
    % FIGURE Xb: The Method (Normalized Geometric Plot)
    % --------------------------------------------------------------------
    if ~isempty(diag)
        fig_xb = figure('Name', 'Figure Xb: Geometric Method', 'Position', [750, 100, 600, 400], 'Color', 'w');
        hold on; grid on;

        % 1. Extract the span of interest from the peak to the far anchor
        idx_span = diag.peak_idx : diag.end_idx;
        x_span = diag.bin_centers(idx_span);
        y_span = diag.counts_smooth(idx_span);

        % 2. Normalize both axes to [0, 1] exactly as described in the paper
        x_n = (x_span - x_span(1)) / (x_span(end) - x_span(1));
        y_n = (y_span - y_span(end)) / (y_span(1) - y_span(end));

        % 3. Plot the normalized histogram curve and the chord
        plot(x_n, y_n, 'k-', 'LineWidth', 2, 'DisplayName', 'Normalized Histogram');
        plot([0 1], [1 0], 'b--', 'LineWidth', 1.5, 'DisplayName', 'Chord');

        % 4. Find the max deviation coordinate and compute the projection
        t_rel_idx = find(idx_span == diag.thresh_idx, 1);
        if ~isempty(t_rel_idx)
            x0 = x_n(t_rel_idx);
            y0 = y_n(t_rel_idx);

            % Perpendicular projection of (x0, y0) onto the line y = -x + 1
            x_proj = (x0 - y0 + 1) / 2;
            y_proj = (-x0 + y0 + 1) / 2;

            % Plot the perpendicular deviation line and threshold point
            plot([x0 x_proj], [y0 y_proj], 'r-', 'LineWidth', 2, 'DisplayName', 'Max Deviation ($d_b$)');
            scatter(x0, y0, 60, 'r', 'filled', 'DisplayName', 'Optimal Threshold ($\theta^{\star}$)');
        end

        xlabel('Normalized Score Bin ($\tilde{x}_b$)', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('Normalized Count ($\tilde{y}_b$)', 'Interpreter', 'latex', 'FontSize', 12);
        title('Figure Xb: Corner Detection via Max Deviation', 'FontSize', 14);
        axis([-0.05 1.05 -0.05 1.05]);
        axis square;
        
        lgd2 = legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);
        lgd2.ItemTokenSize = [20, 18];
        set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
    end

    %% --------------------------------------------------------------------
    % FIGURE Y: The Results (Time-Series Thresholding)
    % --------------------------------------------------------------------
    if ~isempty(theta_series)
        fig_y = figure('Name', 'Figure Y: Dynamic Threshold', 'Position', [100, 550, 800, 350], 'Color', 'w');
        hold on; grid on;

        num_frames = length(theta_series);
        frame_axis = 1:num_frames;

        % Plot background scatter if provided (to show the noise floor changing)
        if ~isempty(scatter_data) && length(scatter_data) == num_frames
            for f = 1:num_frames
                if ~isempty(scatter_data{f})
                    % Downsample heavily for plotting speed/cleanliness
                    f_scores = scatter_data{f};
                    n_pts = min(150, length(f_scores));
                    ds_scores = f_scores(randperm(length(f_scores), n_pts));
                    
                    scatter(repmat(f, n_pts, 1), ds_scores, 10, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.3, 'HandleVisibility','off');
                end
            end
            % Add a dummy scatter for the legend
            scatter(NaN, NaN, 10, [0.7 0.7 0.7], 'filled', 'DisplayName', 'Active Scores (Sampled)');
        end

        % Plot the threshold trend
        plot(frame_axis, theta_series, 'r-', 'LineWidth', 2, 'DisplayName', 'Dynamic Threshold ($\theta^{\star}$)');

        xlabel('Time [Frame]', 'FontSize', 12);
        ylabel('Threshold ($\theta^{\star}$)', 'Interpreter', 'latex', 'FontSize', 12);
        title('Figure Y: Adaptive Threshold Tracking Scene Dynamics', 'FontSize', 14);
        ylim([0 1]);
        xlim([1 num_frames]);
        
        lgd3 = legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);
        set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
    end

end
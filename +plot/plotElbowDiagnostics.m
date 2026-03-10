function fig = plotElbowDiagnostics(auto_threshold, diagnostics, varargin)
% PLOTELBOWDIAGNOSTICS  Two-panel diagnostic figure for elbow threshold.
%
%   FIG = PLOTELBOWDIAGNOSTICS(AUTO_THRESHOLD, DIAGNOSTICS) produces a
%   two-panel figure showing (left) the raw mean NN distance vs.
%   threshold curve with the selected elbow marked, and (right) the
%   normalized curve with the chord and perpendicular distance peak.
%
%   FIG = PLOTELBOWDIAGNOSTICS(..., 'Name', Value) accepts optional
%   name-value arguments:
%     'Title'       - ('') String prepended to subplot titles.
%     'FigPosition' - ([100 100 900 350]) Figure position vector.
%
%   Inputs:
%     auto_threshold - Scalar elbow threshold from findElbowThreshold.
%     diagnostics    - Struct returned by findElbowThreshold (second
%                      output). Must contain fields: th_vec, mean_nn,
%                      th_norm, nn_norm, perp_dist, elbow_idx.
%
%   Outputs:
%     fig - Figure handle.
%
%   See also: process.findElbowThreshold

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'Title', '', @ischar);
    addParameter(p, 'FigPosition', [100 100 900 350], ...
        @(x) isnumeric(x) && numel(x) == 4);
    parse(p, varargin{:});
    label_prefix = p.Results.Title;
    fig_pos      = p.Results.FigPosition;

    % Unpack
    th_g = diagnostics.th_vec;
    nn_g = diagnostics.mean_nn;
    th_n = diagnostics.th_norm;
    nn_n = diagnostics.nn_norm;
    ei   = diagnostics.elbow_idx;

    if ~isempty(label_prefix)
        label_prefix = [label_prefix ' — '];
    end

    fig = figure('Position', fig_pos);

    % --- Left panel: raw curve ---
    subplot(1, 2, 1);
    plot(th_g, nn_g, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
    hold on;
    xline(auto_threshold, 'r--', 'LineWidth', 1.5);
    plot(auto_threshold, nn_g(ei), 'rp', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'r');
    xlabel('Threshold');
    ylabel('Mean NN distance [px]');
    title(sprintf('%sRaw curve — elbow @ %.4f', label_prefix, ...
        auto_threshold));
    grid on;

    % --- Right panel: normalized + chord ---
    subplot(1, 2, 2);
    plot(th_n, nn_n, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
    hold on;
    plot([th_n(1) th_n(end)], [nn_n(1) nn_n(end)], 'k--', ...
        'LineWidth', 1.0);
    plot(th_n(ei), nn_n(ei), 'rp', ...
        'MarkerSize', 14, 'MarkerFaceColor', 'r');

    % Draw perpendicular line from elbow to chord
    p1 = [th_n(1), nn_n(1)];
    v  = [th_n(end), nn_n(end)] - p1;
    w  = [th_n(ei), nn_n(ei)] - p1;
    proj_len = dot(w, v) / dot(v, v);
    foot = p1 + proj_len * v;
    plot([th_n(ei) foot(1)], [nn_n(ei) foot(2)], 'r-', ...
        'LineWidth', 1.0);

    xlabel('Threshold (normalized)');
    ylabel('Mean NN dist (normalized)');
    title(sprintf('%sNormalized — max \\perp dist', label_prefix));
    grid on;

end
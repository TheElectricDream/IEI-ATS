function [] = showPersistenceMatches3D(rowsA, colsA, valsA, ...
        rowsB, colsB, valsB, matchIdx, name, show)
% SHOWPERSISTENCEMATCHES3D  Feature-matching view of cross-frame persistence.
%
%   SHOWPERSISTENCEMATCHES3D(ROWSA, COLSA, VALSA, ROWSB, COLSB, VALSB,
%   MATCHIDX, NAME, SHOW) plots two event point clouds side by side in a
%   single 3D axes and connects a random subset of matched point pairs with
%   lines, in the style of a feature-matching montage. The left cloud is
%   frame A; the right cloud is frame B, shifted along the X axis.
%
%   Inputs:
%     rowsA, colsA, valsA - [Na x 1] row, column and value of frame-A points.
%     rowsB, colsB, valsB - [Nb x 1] row, column and value of frame-B points.
%     matchIdx            - [Na x 1] index into the frame-B arrays giving the
%                           match for each frame-A point; NaN where a point
%                           has no persistence match.
%     name                - Filename stem for export ('' / [] skips export).
%     show                - Logical; if false the figure is built invisibly.
%
%   Outputs:
%     None. Produces a figure and, when NAME is non-empty, an exported PDF.
%
%   Notes:
%     - Both clouds live in one axes so the connecting lines share the 3D
%       coordinate frame and stay valid under rotation; frame B is offset by
%       (1 + GAP) in normalised X to read as two plots side by side.
%     - Up to 10 matched pairs are drawn, sampled with a local RandStream so
%       the figure is reproducible without touching the global RNG.
%     - X and Y are normalised by the 640x480 sensor resolution to match
%       plot.showScatterPlotOfRuleMaps3D; values are assumed pre-normalised.
%
%   See also: plot.showScatterPlotOfRuleMaps3D, journal.exportTight3DScatterPlots

    % --- 1. Style and layout constants --------------------------------
    colorContext = [0.65, 0.65, 0.65];   % subtle grey for non-matched points
    colorA       = [0.00, 0.45, 0.74];   % deep blue    - frame A
    colorB       = [0.85, 0.33, 0.10];   % burnt orange - frame B
    colorLine    = [0.30, 0.30, 0.30];   % connecting lines
    lineAlpha    = 0.60;
    nShow        = 10;                    % number of matches to draw
    gap          = 0.12;                  % horizontal gap between clouds
    W = 640;  H = 480;                    % sensor resolution

    xOffset = 1 + gap;                    % frame-B shift in normalised X

    % --- 2. Normalise coordinates -------------------------------------
    xA = colsA(:) / W;   yA = rowsA(:) / H;   zA = valsA(:);
    xB = colsB(:) / W;   yB = rowsB(:) / H;   zB = valsB(:);
    xBoff = xB + xOffset;                 % shifted copy for display

    % --- 3. Select a random subset of valid matches -------------------
    validA = find(~isnan(matchIdx(:)));
    k      = min(nShow, numel(validA));
    s      = RandStream('twister', 'Seed', 0);
    sel    = validA(randperm(s, numel(validA), k));
    selB   = matchIdx(sel);               % matched indices into frame B

    % --- 4. Figure -----------------------------------------------------
    if show
        fig = figure();
    else
        fig = figure('Visible', 'off');
    end
    ax = axes('Parent', fig);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); rotate3d(ax, 'on');

    % Context points (all points, faint)
    scatter3(ax, xA,    yA, zA, 15, colorContext, 'filled', 'MarkerFaceAlpha', 0.15);
    scatter3(ax, xBoff, yB, zB, 15, colorContext, 'filled', 'MarkerFaceAlpha', 0.15);

    % Connecting lines between matched pairs
    hL = gobjects(1);
    for i = 1:k
        a = sel(i);  b = selB(i);
        h = plot3(ax, [xA(a), xBoff(b)], [yA(a), yB(b)], [zA(a), zB(b)], ...
            '-', 'Color', [colorLine, lineAlpha], 'LineWidth', 1.2);
        if i == 1, hL = h; end
    end

    % Matched endpoints (drawn last so they sit on top)
    hA = scatter3(ax, xA(sel),    yA(sel),   zA(sel),   70, colorA, 'd', ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    hB = scatter3(ax, xBoff(selB), yB(selB),  zB(selB),  70, colorB, 'o', ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);

    % --- 5. Axes cosmetics --------------------------------------------
    xlabel(ax, 'X_{norm} [-]');
    ylabel(ax, 'Y_{norm} [-]');
    zlabel(ax, '\rho_{norm} [-]');
    view(ax, 3);
    axis(ax, 'tight');
    legend([hA, hB, hL], {'Frame k', 'Frame k+1', 'Matches'}, ...
        'Location', 'northeast', 'Box', 'off');
    set(ax,  'FontSize', 16, 'FontName', 'Times New Roman');
    set(fig, 'DefaultTextFontName', 'Times New Roman', ...
             'DefaultAxesFontName', 'Times New Roman', 'Color', 'white');

    % --- 6. Export -----------------------------------------------------
    if nargin >= 8 && ~isempty(name)
        basePath = ['/home/alexandercrain/Dropbox/Graduate Documents/' ...
            'Doctor of Philosophy/Publications/Journals/' ...
            'AIAA Journal of Spacecraft and Rockets/' ...
            'Event_Based_Spacecraft_Representation_Using_Inter_Event_' ...
            'Interval_Adaptive_Time_Surfaces/results/generated-figures/'];
        journal.exportTight3DScatterPlots(gcf, [basePath name]);
    end
end
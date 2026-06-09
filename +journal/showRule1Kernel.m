function fig = showRule1Kernel(sorted_x, sorted_y, sorted_t, ...
    t_interval, imgSz, r_s, varargin)
% SHOWRULE1KERNEL  3D event scatter with spatial-density kernel overlay.
%
%   FIG = SHOWRULE1KERNEL(SORTED_X, SORTED_Y, SORTED_T, T_INTERVAL,
%   IMGSZ, R_S) renders the events of a single frame window in
%   normalized [0, 1]^3 coordinates and overlays a transparent ball
%   of radius R_S centred on a randomly selected event. Events that
%   fall inside the ball are highlighted; the remaining events are
%   shown as a faint background cloud. Intended as a methods figure
%   illustrating Rule 1 of the three-rule coherence filter.
%
%   FIG = SHOWRULE1KERNEL(..., 'Name', Value) accepts optional
%   name-value arguments:
%     'CenterIdx'   - ([]) Linear index into the sorted event arrays
%                     to centre the kernel on. When empty, a random
%                     event is selected.
%     'NumOutside'  - (5000) Maximum number of background ("outside")
%                     events scattered for legibility. The "inside"
%                     events are never subsampled.
%     'RandomSeed'  - ([]) RNG seed for reproducible figures. When
%                     empty, the current RNG state is used.
%     'FigPosition' - ([100 100 720 640]) Figure position vector.
%     'ExportPath'  - ('') If non-empty, the figure is saved to this
%                     path via exportgraphics at the end of the call.
%
%   Inputs:
%     sorted_x   - [N x 1] Row coordinates of events (1-indexed).
%     sorted_y   - [N x 1] Column coordinates of events (1-indexed).
%     sorted_t   - [N x 1] Event timestamps [s], any monotonic offset.
%     t_interval - Scalar frame window duration [s].
%     imgSz      - [1 x 2] Image dimensions [nRows, nCols].
%     r_s        - Scalar normalized search radius, matching
%                  coh_params.r_s in coherence.computeCoherenceMask.
%
%   Outputs:
%     fig - Figure handle.
%
%   Notes:
%     - Coordinates are normalized by dividing rows by imgSz(1),
%       columns by imgSz(2), and (t - min(t)) by t_interval, matching
%       the [0, 1]^3 space used by Rule 1 of the coherence filter.
%     - 'axis equal' is enforced so the kernel renders as a ball.
%       Without it, MATLAB stretches axes independently and the
%       sphere appears as an ellipsoid even though the underlying
%       geometric object is spherical in normalized coordinates.
%     - Only the background ("outside") cloud is thinned by
%       'NumOutside'. The "inside" events are always plotted in full
%       since they are the quantity the rule integrates over.
%     - Coordinates: x = row, y = col, following repository
%       convention.
%
%   See also: coherence.computeCoherenceMask,
%             filters.findSpatialNeighbours

    % ----------------------------------------------------------------
    % 0. Parse optional arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'CenterIdx',   [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'NumOutside',  5000, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'RandomSeed',  [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'FigPosition', [100 100 720 640], ...
        @(x) isnumeric(x) && numel(x) == 4);
    addParameter(p, 'ExportPath',  '', ...
        @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    centerIdx   = p.Results.CenterIdx;
    numOutside  = p.Results.NumOutside;
    randomSeed  = p.Results.RandomSeed;
    figPosition = p.Results.FigPosition;
    exportPath  = char(p.Results.ExportPath);

    if ~isempty(randomSeed)
        rng(10000);
    end

    % ----------------------------------------------------------------
    % 1. Normalize event coordinates to [0, 1]^3
    % ----------------------------------------------------------------
    N  = numel(sorted_x);
    if N == 0
        error('showRule1Kernel:emptyFrame', ...
            'No events provided for visualization.');
    end

    t0 = min(sorted_t);

    ex = double(sorted_x) / imgSz(1);
    ey = double(sorted_y) / imgSz(2);
    et = (double(sorted_t) - t0) / t_interval;

    % ----------------------------------------------------------------
    % 2. Select kernel centre event
    % ----------------------------------------------------------------
    if isempty(centerIdx)
        centerIdx = randi(N);
    elseif centerIdx > N
        error('showRule1Kernel:invalidCenter', ...
            'CenterIdx (%d) exceeds the number of events (%d).', ...
            centerIdx, N);
    end
    c = [ex(centerIdx), ey(centerIdx), et(centerIdx)];

    % ----------------------------------------------------------------
    % 3. Partition events by kernel membership
    % ----------------------------------------------------------------
    d                  = sqrt((ex - c(1)).^2 + (ey - c(2)).^2 + ...
                              (et - c(3)).^2);
    inside             = d <= r_s;
    inside(centerIdx)  = false;             % exclude centre from "inside"
    isCenter           = false(N, 1);
    isCenter(centerIdx) = true;

    % Thin the outside cloud for legibility (inside events unchanged)
    outsideMask = ~inside & ~isCenter;
    outsideIdx  = find(outsideMask);
    if numel(outsideIdx) > numOutside
        keep       = randperm(numel(outsideIdx), numOutside);
        outsideIdx = outsideIdx(keep);
    end

    % ----------------------------------------------------------------
    % 4. Build figure
    % ----------------------------------------------------------------
    fig = figure('Color', 'w', 'Position', figPosition);
    ax  = axes('Parent', fig);
    hold(ax, 'on');

    % Outside events (faint grey background)
    scatter3(ax, ex(outsideIdx), ey(outsideIdx), et(outsideIdx), ...
        6, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.8);

    % Inside events (blue, fully drawn)
    scatter3(ax, ex(inside), ey(inside), et(inside), ...
        30, 'b', 'filled');

    % Centre event (red, outlined)
    scatter3(ax, c(1), c(2), c(3), 100, 'r', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);

    % Transparent kernel ball of radius r_s
    [sx, sy, sz] = sphere(60);
    surf(ax, r_s*sx + c(1), r_s*sy + c(2), r_s*sz + c(3), ...
        'FaceColor', 'r', 'FaceAlpha', 0.12, 'EdgeColor', 'none');

    % ----------------------------------------------------------------
    % 5. Format axes and lighting
    % ----------------------------------------------------------------
    axis(ax, 'equal');
    grid(ax, 'on');
    box(ax, 'on');
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');

    xlabel(ax, 'X_{norm} [-]');
    ylabel(ax, 'Y_{norm} [-]');
    zlabel(ax, 't_{norm} [-]');
    % title(ax, sprintf(['Rule 1 kernel at event %d: ' ...
    %     'r_s = %.3f, %d neighbours inside'], ...
    %     centerIdx, r_s, nnz(inside)));
    legend('Excluded Events', 'Included Events', 'Seed Event','Kernel Sphere', 'FontSize', 16, 'FontName', 'Times New Roman');
    view( 3); % Forces the 3D perspective to prevent locking to 2D
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
    set(gcf, 'DefaultTextFontName', 'Times New Roman', 'DefaultAxesFontName', 'Times New Roman');
    journal.exportTight3DScatterPlots(gcf,'/home/alexandercrain/Dropbox/Graduate Documents/Doctor of Philosophy/Publications/Journals/AIAA Journal of Spacecraft and Rockets/Event_Based_Spacecraft_Representation_Using_Inter_Event_Interval_Adaptive_Time_Surfaces/results/generated-figures/Spherical-Kernel-Sample-Nominal-Rot.pdf');


end
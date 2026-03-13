function R = computeClarkEvans(map, varargin)
% CALCULATECLARKEVANS  Clark-Evans nearest-neighbour aggregation index.
%
%   R = CALCULATECLARKEVANS(MAP) extracts active pixel coordinates from
%   MAP (map > 0) and computes the Clark-Evans R statistic — the ratio of
%   the observed mean nearest-neighbour distance to the distance expected
%   under complete spatial randomness (CSR).
%
%   R = CALCULATECLARKEVANS(MAP, 'Name', Value) accepts optional
%   name-value arguments:
%     'ImageSize' - ([H W], derived from MAP) Area of the observation
%                   window in pixels. Defaults to size(map), so override
%                   only when MAP is a cropped ROI of a larger sensor.
%
%   Inputs:
%     map - [H x W] 2D value map. Active pixels are those with map > 0.
%
%   Outputs:
%     R - Scalar Clark-Evans index.
%           R ≈ 1 : randomly distributed (CSR)
%           R < 1 : clustered
%           R > 1 : over-dispersed / regular
%
%   Notes:
%     - Nearest-neighbour search uses KNNSEARCH (k = 2, self-excluded)
%       rather than an O(N^2) explicit loop, giving O(N log N) scaling.
%     - Returns NaN and emits a warning when fewer than 2 active pixels
%       are present, as the statistic is undefined in that case.
%     - Reference: P. J. Clark and F. C. Evans, "Distance to Nearest
%       Neighbor as a Measure of Spatial Relationships in Populations,"
%       Ecology, vol. 35, no. 4, pp. 445-453, Oct. 1954,
%       doi: 10.2307/1931034.
%
%   See also: plot.showSpreadDistribution, stats.findSpatialNeighbours

    % ----------------------------------------------------------------
    % 0. Parse optional arguments
    % ----------------------------------------------------------------
    [H, W] = size(map);
    p = inputParser;
    addParameter(p, 'ImageSize', [H, W], ...
        @(x) isnumeric(x) && numel(x) == 2 && all(x > 0));
    parse(p, varargin{:});
    img_sz = p.Results.ImageSize;

    % ----------------------------------------------------------------
    % 1. Extract active pixel coordinates [col, row]
    % ----------------------------------------------------------------
    [rows, cols] = find(map > 0);
    pts = [cols, rows];
    N   = size(pts, 1);

    if N < 2
        warning('calculateClarkEvans:insufficientPoints', ...
            'Fewer than 2 active pixels — R statistic is undefined.');
        R = NaN;
        return;
    end

    % ----------------------------------------------------------------
    % 2. Nearest-neighbour distances  (O(N log N) via kd-tree)
    % ----------------------------------------------------------------
    % k = 2: first neighbour is the point itself; second is the nearest
    % other point.
    [~, D_nn] = knnsearch(pts, pts, 'K', 2);
    r_obs = mean(D_nn(:, 2));

    % ----------------------------------------------------------------
    % 3. Clark-Evans statistic
    % ----------------------------------------------------------------
    A     = img_sz(1) * img_sz(2);     % observation window area [px^2]
    r_exp = 0.5 / sqrt(N / A);        % expected NND under CSR
    R     = r_obs / r_exp;

end
function [splatted_map] = splatTimestampMap(x, y, t, imgSz, sigma)
% SPLATTIMESTAMPMAP  Gaussian normalized convolution of event timestamps.
%
%   SPLATTED_MAP = SPLATTIMESTAMPMAP(X, Y, T, IMGSZ, SIGMA) spreads
%   event timestamps spatially using Gaussian-weighted normalized
%   convolution. The result is a smooth per-pixel timestamp estimate
%   that accounts for unobserved pixels.
%
%   Reference:
%       Knutsson, H. and Westin, C.-F. (1993), "Normalized and
%       Differential Convolution," Proc. IEEE CVPR, pp. 515-523.
%
%   Inputs:
%     x, y   - [N x 1] Pixel coordinates (row, col).
%     t      - [N x 1] Timestamp values [s].
%     imgSz  - [1 x 2] Image dimensions [nRows, nCols].
%     sigma  - Scalar Gaussian kernel standard deviation [pixels].
%
%   Outputs:
%     splatted_map - [imgSz] Spatially smoothed timestamp map.
%
%   Algorithm:
%     1. Scatter timestamps onto a raw map (latest event wins).
%     2. Build a binary observation mask.
%     3. Gaussian-blur both the raw map and the mask.
%     4. Divide: splatted = blurred_values / blurred_weights.
%
%   Notes:
%     - Regions with no nearby events are set to zero.
%     - Coordinates: x = row, y = col, sub2ind(imgSz, x, y).
%
%   See also: stats.spreadEventsSpatially

    % ----------------------------------------------------------------
    % 0. Initialize
    % ----------------------------------------------------------------
    if ~isa(t, 'single'), t = single(t); end

    raw_map  = zeros(imgSz, 'single');
    mask_map = zeros(imgSz, 'single');

    % ----------------------------------------------------------------
    % 1. Scatter timestamps (latest overwrites at collisions)
    % ----------------------------------------------------------------
    linear_idx = sub2ind(imgSz, x, y);
    raw_map(linear_idx)  = t;
    mask_map(linear_idx) = 1.0;

    % ----------------------------------------------------------------
    % 2. Gaussian normalized convolution
    % ----------------------------------------------------------------
    k_size = 2 * ceil(2 * sigma) + 1;

    blurred_vals = imgaussfilt(raw_map, sigma, ...
        'FilterSize', k_size, 'Padding', 0);
    blurred_weights = imgaussfilt(mask_map, sigma, ...
        'FilterSize', k_size, 'Padding', 0);

    % ----------------------------------------------------------------
    % 3. Normalize and clean
    % ----------------------------------------------------------------
    splatted_map = blurred_vals ./ (blurred_weights + eps('single'));
    splatted_map(blurred_weights < 1e-4) = 0;

end
function [voxelOccupancy, xEdges, yEdges, tEdges] = ...
    discretizeEventsToVoxels(x, y, t, opts)
% DISCRETIZEEVENTSTOVOXELS  Bin events into a 3D spatiotemporal voxel grid.
%
%   [VOXELOCCUPANCY, XEDGES, YEDGES, TEDGES] =
%   DISCRETIZEEVENTSTOVOXELS(X, Y, T, OPTS) assigns each event to a
%   voxel in a regular 3D grid and returns a binary occupancy volume.
%
%   Inputs:
%     x, y - [N x 1] Pixel coordinates (row, col).
%     t    - [N x 1] Timestamps [s].
%     opts - Struct with fields:
%              .dx    - Spatial bin width in x [pixels] (default: 1)
%              .dy    - Spatial bin width in y [pixels] (default: 1)
%              .dt    - Temporal bin width [s]
%              .imgSz - [1 x 2] Image dimensions [nRows, nCols].
%                       Used to define spatial bin edges.
%
%   Outputs:
%     voxelOccupancy - [Nx x Ny x Nt] Logical 3D array. True where
%                      at least one event falls in the voxel.
%     xEdges         - [Nx+1 x 1] Bin edges along x dimension.
%     yEdges         - [Ny+1 x 1] Bin edges along y dimension.
%     tEdges         - [Nt+1 x 1] Bin edges along t dimension.
%
%   Algorithm:
%     1. Define spatial bin edges from imgSz (or from opts.dx/dy).
%     2. Define temporal bin edges from the data range and opts.dt.
%     3. Discretize each event into (ix, iy, it) voxel indices.
%     4. Accumulate counts and convert to binary occupancy.
%
%   Notes:
%     - Spatial bins default to 1-pixel resolution when imgSz is
%       provided. For coarser binning, increase opts.dx / opts.dy.
%     - Events outside the bin edges are dropped (NaN from
%       discretize).
%
%   See also: discretize, accumarray

    % ----------------------------------------------------------------
    % 0. Define bin edges
    % ----------------------------------------------------------------
    % Spatial edges from image dimensions
    if isfield(opts, 'imgSz')
        xEdges = 1:(opts.imgSz(1) + 1);
        yEdges = 1:(opts.imgSz(2) + 1);
    else
        xEdges = (floor(min(x) / opts.dx) * opts.dx) : opts.dx : ...
            (ceil(max(x) / opts.dx) * opts.dx + opts.dx);
        yEdges = (floor(min(y) / opts.dy) * opts.dy) : opts.dy : ...
            (ceil(max(y) / opts.dy) * opts.dy + opts.dy);
    end

    % Temporal edges from data range
    tEdges = (floor(min(t) / opts.dt) * opts.dt) : opts.dt : ...
        (ceil(max(t) / opts.dt) * opts.dt + opts.dt);

    % ----------------------------------------------------------------
    % 1. Discretize events into voxel indices
    % ----------------------------------------------------------------
    ix = discretize(x, xEdges);
    iy = discretize(y, yEdges);
    it = discretize(t, tEdges);

    % Drop events outside bin edges
    keep = ~isnan(ix) & ~isnan(iy) & ~isnan(it);
    subs = [ix(keep), iy(keep), it(keep)];

    % ----------------------------------------------------------------
    % 2. Accumulate and convert to occupancy
    % ----------------------------------------------------------------
    sz = [numel(xEdges)-1, numel(yEdges)-1, numel(tEdges)-1];
    voxelCount = accumarray(subs, 1, sz, @sum, 0);
    voxelOccupancy = voxelCount > 0;

end
function [pointCloud] = generateMeshFromFrame(map)
% GENERATEMESHFROMFRAME Convert a 2D value map to an [H*W x 3] point cloud.
%
%     pointCloud = generateMeshFromFrame(map)
%
%     Inputs:
%       map - [H x W] Matrix of values (e.g., a time surface).
%
%     Outputs:
%       pointCloud - [H*W x 3] Array of [x, y, z] points, where x is
%                    the column index (1..W), y is the row index
%                    (1..H), and z is the map value. Rows are ordered
%                    column-major (y varies fastest). Zero-valued
%                    pixels yield NaN in the z column only; their x
%                    and y entries remain valid.
%
%     See also MESHGRID, PLOT.MAPTOSCATTERPLOT, PLOT.MAPTOSURFPLOT.

    map(map == 0) = nan;
    [x, y] = meshgrid(1:size(map, 2), 1:size(map, 1));
    pointCloud = [x(:), y(:), map(:)];
end

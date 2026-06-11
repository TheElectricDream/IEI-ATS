function Map_Filled = localized_mls_fill(Map, SearchRadius)
% Map: 640x480 matrix with NaNs
if nargin < 2, SearchRadius = 5; end % Pixel neighborhood radius

[Rows, Cols] = size(Map);
Map_Filled = Map;
Hole_Indices = find(Map==0);

% Convert grid coordinates into explicit (X, Y, Val) structures
[X, Y] = meshgrid(1:Cols, 1:Rows);

% Loop over each hole pixel to evaluate a localized robust local plane
for k = 1:length(Hole_Indices)
    idx = Hole_Indices(k);
    [r, c] = ind2sub([Rows, Cols], idx);

    % Define local bounding box around the hole point
    r_min = max(1, r - SearchRadius); r_max = min(Rows, r + SearchRadius);
    c_min = max(1, c - SearchRadius); c_max = min(Cols, c + SearchRadius);

    % Extract localized neighbor patches
    Sub_X = X(r_min:r_max, c_min:c_max);
    Sub_Y = Y(r_min:r_max, c_min:c_max);
    Sub_V = Map(r_min:r_max, c_min:c_max);

    % Filter out other NaNs in this window to gather valid anchor points
    Valid_Mask = ~isnan(Sub_V);
    if sum(Valid_Mask(:)) < 6
        continue; % Skip if there aren't enough points to fit a model
    end

    Pts_X = Sub_X(Valid_Mask);
    Pts_Y = Sub_Y(Valid_Mask);
    Pts_V = Sub_V(Valid_Mask);

    % Fit an Edge-Preserving Robust Linear Model (M-estimator using Huber weights)
    % This rejects points on the far side of a sharp boundary step as outliers.
    Model = robustfit([Pts_X, Pts_Y], Pts_V, 'huber', 1.5);

    % Evaluate the fitted local surface exactly at the missing coordinate (c, r)
    Map_Filled(r, c) = Model(1) + Model(2)*c + Model(3)*r;
end
end
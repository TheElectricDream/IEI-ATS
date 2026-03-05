function x_out = linearRemap(x, tau_min, tau_max)
% LINEARREMAP  Linear mapping between specified bounds.
%
%   X_OUT = LINEARREMAP(X, TAU_MIN, TAU_MAX) remaps values in X to
%   [TAU_MIN, TAU_MAX] using min-max linear scaling.
%
%   Inputs:
%     x       - Numeric array of values to remap.
%     tau_min - Scalar lower bound of the output range.
%     tau_max - Scalar upper bound of the output range.
%
%   Outputs:
%     x_out   - Numeric array mapped to [tau_min, tau_max].
%
%   Notes:
%     - Returns midpoint if input is constant (avoids division by 0).
%
%   See also: process.sigmoidRemap

    x = double(x);

    if isempty(x)
        x_out = zeros(size(x));
        return;
    end

    xmin = min(x(:));
    xmax = max(x(:));

    if xmax == xmin
        x_out = repmat((tau_min + tau_max) / 2, size(x));
        return;
    end

    x_norm = (x - xmin) / (xmax - xmin);
    x_out = tau_min + x_norm .* (tau_max - tau_min);

end
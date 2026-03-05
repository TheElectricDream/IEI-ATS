function x_out = sigmoidRemap(x, tau_min, tau_max, k)
% SIGMOIDREMAP  Sigmoid mapping between specified bounds.
%
%   X_OUT = SIGMOIDREMAP(X, TAU_MIN, TAU_MAX) remaps values in X to
%   [TAU_MIN, TAU_MAX] using a sigmoid-shaped transfer function.
%
%   X_OUT = SIGMOIDREMAP(X, TAU_MIN, TAU_MAX, K) specifies the
%   steepness parameter K (default: 6). Larger K = steeper sigmoid.
%
%   Inputs:
%     x       - Numeric array of values to remap.
%     tau_min - Scalar lower bound of the output range.
%     tau_max - Scalar upper bound of the output range.
%     k       - (Optional) Steepness parameter. (default: 6)
%
%   Outputs:
%     x_out   - Numeric array mapped to [tau_min, tau_max].
%
%   Notes:
%     - Returns midpoint if input is constant (avoids division by 0).
%
%   See also: process.linearRemap

    if nargin < 4
        k = 6;
    end

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

    x_norm = (x - xmin) / (xmax - xmin) * 2 * k - k;
    s = 1 ./ (1 + exp(-x_norm));
    x_out = tau_min + s .* (tau_max - tau_min);

end
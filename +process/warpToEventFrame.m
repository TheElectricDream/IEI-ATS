function [warped_frame, H] = warpToEventFrame(frame, K_zed, K_event, t, varargin)
% WARPTOEVENTFRAME  Warp an RGB frame to align with the event camera viewpoint.
%
%   [WARPED_FRAME, H] = WARPTOEVENTFRAME(FRAME, K_ZED, K_EVENT, T)
%   applies a planar homography to transform a frame captured by the
%   ZED camera so that it appears as if captured from the event
%   camera's optical centre. The two cameras are assumed to be
%   parallel (no relative rotation).
%
%   [WARPED_FRAME, H] = WARPTOEVENTFRAME(..., 'Name', Value) accepts
%   optional name-value arguments:
%     'OutputSize'  - [H W] Output dimensions [rows, cols]. Default
%                     is derived from K_event: [2*cy, 2*cx] where
%                     (cx, cy) is the principal point.
%     'SceneDepth'  - Scalar assumed scene depth Z [metres] for
%                     parallax correction. Default is Inf, which
%                     applies the infinite-depth (pure intrinsics)
%                     homography. For SPOT testbed distances of
%                     1–3 m, set this to the approximate target
%                     range to correct for the baseline offset.
%     'Interpolation' - ('bilinear') Interpolation method passed
%                       to imwarp. Options: 'nearest', 'bilinear',
%                       'bicubic'.
%     'FillValue'   - (0) Value for pixels outside the source frame.
%
%   Inputs:
%     frame   - [M x N] or [M x N x 3] Image from the ZED camera.
%     K_zed   - [3 x 3] Intrinsic matrix of the ZED camera.
%     K_event - [3 x 3] Intrinsic matrix of the event camera.
%     t       - [3 x 1] Translation vector FROM the event camera
%               TO the ZED camera, expressed in the event camera's
%               coordinate frame [metres]. A zero vector means the
%               two optical centres coincide.
%
%   Outputs:
%     warped_frame - [H x W] or [H x W x 3] Warped image at the
%                    event camera's resolution and viewpoint.
%     H            - [3 x 3] Homography matrix mapping ZED pixel
%                    coordinates to event camera pixel coordinates
%                    (p_event ~ H * p_zed).
%
%   Algorithm:
%     For parallel cameras (R = I) separated by translation t,
%     a 3-D point at depth Z in the ZED frame projects to:
%
%       p_event ~ K_event * ( Z * inv(K_zed) * p_zed  +  t )
%
%     Exploiting the unit third component of homogeneous pixel
%     coordinates (e3' * p_zed = 1), this factors as:
%
%       H = K_event * ( Z * inv(K_zed)  +  t * e3' )
%
%     At infinite depth (Z -> Inf) the translation term vanishes
%     and the mapping reduces to the intrinsics-only homography:
%
%       H_inf = K_event * inv(K_zed)
%
%   Notes:
%     - The homography is exact for a fronto-parallel plane at
%       depth Z. For non-planar scenes, it is an approximation
%       that improves as the scene depth increases relative to
%       the baseline ||t||.
%     - The homography is defined up to scale; it is normalized
%       so that H(3,3) = 1 before being passed to projective2d.
%     - MATLAB's projective2d uses the TRANSPOSE convention:
%       [u' v' w'] = [u v 1] * T, where T = H'. This function
%       handles the transposition internally.
%
%   Example:
%     % ZED2 left camera intrinsics (1280x720)
%     K_zed = [527.5  0     640.0;
%               0    527.5  360.0;
%               0      0      1  ];
%
%     % DVXplorer Micro intrinsics (640x480)
%     K_event = [320.0  0     320.0;
%                  0   320.0  240.0;
%                  0     0      1  ];
%
%     % Translation: event cam is 5 cm left, 2 cm above ZED
%     t = [-0.05; 0.02; 0.0];
%
%     warped = warpToEventFrame(rgb_frame, K_zed, K_event, t, ...
%                  'SceneDepth', 2.0, 'OutputSize', [480, 640]);
%
%   See also: projective2d, imwarp, imref2d

    % ----------------------------------------------------------------
    % 0. Parse optional arguments
    % ----------------------------------------------------------------
    p = inputParser;
    addRequired(p, 'frame',   @(x) isnumeric(x) && ndims(x) >= 2);
    addRequired(p, 'K_zed',   @(x) isnumeric(x) && isequal(size(x), [3 3]));
    addRequired(p, 'K_event', @(x) isnumeric(x) && isequal(size(x), [3 3]));
    addRequired(p, 't',       @(x) isnumeric(x) && numel(x) == 3);

    % Derive default output size from K_event principal point
    cx = K_event(1,3);
    cy = K_event(2,3);
    default_size = round([2*cy, 2*cx]);

    addParameter(p, 'OutputSize',    default_size, ...
        @(x) isnumeric(x) && numel(x) == 2 && all(x > 0));
    addParameter(p, 'SceneDepth',    Inf,          ...
        @(x) isscalar(x) && isnumeric(x) && x > 0);
    addParameter(p, 'Interpolation', 'bilinear',   ...
        @(x) ismember(x, {'nearest', 'bilinear', 'bicubic'}));
    addParameter(p, 'FillValue',     0,            ...
        @(x) isnumeric(x) && isscalar(x));
    parse(p, frame, K_zed, K_event, t, varargin{:});

    out_sz  = p.Results.OutputSize;
    Z       = p.Results.SceneDepth;
    interp  = p.Results.Interpolation;
    fill_val = p.Results.FillValue;

    % Ensure column vector
    t = t(:);

    % ----------------------------------------------------------------
    % 1. Compute the homography  H = K_event * (Z * inv(K_zed) + t * e3')
    % ----------------------------------------------------------------
    e3 = [0; 0; 1];

    if isinf(Z)
        % Infinite-depth: pure intrinsics mapping
        H = K_event * (K_zed \ eye(3));
    else
        % Finite-depth: includes parallax correction
        H = K_event * (Z * (K_zed \ eye(3)) + t * e3');
    end

    % Normalize so H(3,3) = 1
    H = H ./ H(3,3);

    % ----------------------------------------------------------------
    % 2. Apply the warp using imwarp
    % ----------------------------------------------------------------
    %   projective2d expects the TRANSPOSE: [u' v' w'] = [u v 1] * T
    tform = projective2d(H');

    output_ref = imref2d(out_sz);

    warped_frame = imwarp(frame, tform, interp, ...
        'OutputView', output_ref, ...
        'FillValues', fill_val);

end
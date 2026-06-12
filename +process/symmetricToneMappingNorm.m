function norm_S = symmetricToneMappingNorm(S, scale)
% SYMMETRICTONEMAPPINGNORM Symmetric tanh tone mapping of a signed surface to (0,1).
%
%     norm_S = symmetricToneMappingNorm(S, scale)
%
%     Inputs:
%       S      - Array of signed surface values (any size).
%       scale  - Positive scalar contrast parameter. Larger values give a
%                softer curve with more headroom for extremes; smaller
%                values give a steeper curve with more contrast near zero.
%                Optional; default is 3.0.
%
%     Outputs:
%       norm_S - Array of the same size as S with values in the open
%                interval (0,1). Zero maps to 0.5 (mid-gray); extremes
%                approach 0 and 1 asymptotically.
%
%     Notes:
%       Fully vectorized; operates elementwise on S.
%
%       Sigmoidal compression is a standard tone-mapping operator for
%       high-dynamic-range imagery; see Reinhard et al. (2002),
%       "Photographic Tone Reproduction for Digital Images," ACM Trans.
%       Graphics (SIGGRAPH), Eq. 4.
%
%     See also TANH, MAT2GRAY, RESCALE.

    % --- 1. Defaults ---
    if nargin < 2
        scale = 3.0;
    end
    
    % --- 2. Symmetric tanh mapping ---
    norm_S = 0.5 + 0.5 * tanh(S / scale);
    
end
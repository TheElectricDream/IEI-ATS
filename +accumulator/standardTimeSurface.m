function [visual_surface, last_t_map, last_p_map] = ...
    standardTimeSurface(last_t_map, last_p_map, x, y, t, p, ...
    imgSz, t_current, tau)
% standardTimeSurface Generates a Time Surface (SAE) with motion blur.
%
%   INPUTS:
%   last_t_map  - (Height x Width) Map of the timestamps of the last events.
%   last_p_map  - (Height x Width) Map of the polarity of the last events.
%   x, y, t, p  - Vectors of event data for the current packet.
%   imgSz       - Size of the image [Height, Width].
%   t_current   - The current simulation time (usually max(t) of packet).
%   tau         - Time constant [seconds] for the visual decay (motion blur).
%
%   OUTPUTS:
%   visual_surface - The visualization frame [-1 to 1] with exponential decay.
%   last_t_map     - Updated timestamp map.
%   last_p_map     - Updated polarity map.

    % 1. Update the Global State Maps
    % We only care about the *last* event at each pixel in this packet.
    % Since x, y, t are sorted by time, we can just write them sequentially.
    % However, using linear indices is faster and handles overwrites automatically.
    
    linear_idx = sub2ind(imgSz, x, y);
    
    % Update timestamps (latest event overwrites previous ones at same pixel)
    last_t_map(linear_idx) = t;
    
    % Update polarities (preserve the polarity of the last event)
    % Ensure polarity is signed (-1, 1)
    p_signed = double(p);
    p_signed(p_signed == 0) = -1;
    last_p_map(linear_idx) = p_signed;

    % 2. Calculate the Time Surface (Motion Blur)
    % Formula: S(x,y) = Polarity * exp( - (t_now - t_last) / tau )
    
    % Calculate time delta
    dt_map = t_current - last_t_map;
    
    % Handle potential negative times (if t_current < t_last due to slice issues)
    dt_map(dt_map < 0) = 0;
    
    % Calculate exponential decay value (0 to 1)
    decay_value = exp(-dt_map / tau);
    
    % Apply polarity
    visual_surface = last_p_map .* decay_value;

end
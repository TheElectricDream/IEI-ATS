function [visual_surface, last_t_map, last_p_map] = ...
    standardTimeSurface(last_t_map, last_p_map, x, y, t, p, ...
    imgSz, t_current, tau)
% STANDARDTIMESURFACE Generates a Time Surface using Linear Decay (SAE).
%
%   Unlike exponential decay (which never truly hits zero), this method
%   linearly maps the event age to intensity.
%   - Event happening NOW = 1.0
%   - Event happening 'tau' seconds ago = 0.0
%   - Events older than 'tau' = 0.0 (Hard cutoff)
%
%   INPUTS:
%   tau - Acts as the "Window Size" or "Memory Horizon" (not a half-life).

    % Update the Global State Maps
    linear_idx = sub2ind(imgSz, x, y);
    
    % Update timestamps (latest event overwrites previous ones)
    last_t_map(linear_idx) = t;
    
    % Update polarities
    p_signed = double(p);
    p_signed(p_signed == 0) = -1;
    last_p_map(linear_idx) = p_signed;

    % Calculate Linear Decay Surface
    
    % Calculate Age: How long ago did the event happen?
    age_map = t_current - last_t_map;
    
    % Handle potential negative times (safety for sorting issues)
    age_map(age_map < 0) = 0;
    
    % LINEAR DECAY FORMULA
    % Value = 1.0 - (Age / Horizon)
    % This creates a straight line from White (now) to Black (tau seconds ago)
    decay_value = 1 - (age_map / tau);
    
    % Clip values: Anything older than 'tau' becomes exactly 0 (Black)
    decay_value(decay_value < 0) = 0;
    
    % Apply polarity (-1 to 1 range)
    visual_surface = last_p_map .* decay_value;
    
end
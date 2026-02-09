function [visual_surface, last_t_map, tau_map, last_p_map] = ...
    speedInvariantTimeSurface(last_t_map, tau_map, last_p_map, ...
    x, y, t, p, imgSz, t_current, sits_scale)
% SPEEDINVARIANTTIMESURFACE Generates a surface with velocity-dependent decay.
%
%   Implements the core concept of Manderscheid et al. (2019):
%   The decay constant 'tau' is updated per-pixel based on the inter-event
%   time difference (dt). This ensures the visual "trail" of an object
%   has a consistent spatial length regardless of how fast the object moves.
%
%   INPUTS:
%   last_t_map  - (HxW) State: Timestamp of the previous event.
%   tau_map     - (HxW) State: Current decay constant for each pixel.
%   last_p_map  - (HxW) State: Polarity of last event.
%   x, y, t, p  - Vectors of NEW events.
%   imgSz       - [Height, Width].
%   t_current   - Current simulation time (for decay calculation).
%   sits_scale  - Scaling factor 'K'. Standard choice is between 1 and 3.
%                 (tau = dt * K). Controls the "length" of the trail.
%
%   OUTPUTS:
%   visual_surface - The visualization frame [-1 to 1].
%   last_t_map     - Updated timestamps.
%   tau_map        - Updated local decay map.
%   last_p_map     - Updated polarity.

    %% 1. UPDATE STATE (Per-Pixel Tau Calculation)
    
    linear_idx = sub2ind(imgSz, x, y);
    
    % 1. Retrieve the LAST timestamp at these specific pixels
    %    (Before we overwrite them with the new time)
    t_prev = last_t_map(linear_idx);
    
    % 2. Calculate dt (Inter-event time)
    %    This is the inverse proxy for speed. 
    %    Small dt = Fast speed; Large dt = Slow speed.
    dt_events = t - t_prev;
    
    % 3. Filter Invalid dt
    %    If dt is Inf (first event ever) or too large (>0.5s), it's likely 
    %    background noise or a new object, not part of a continuous flow. 
    %    We clamp these to a reasonable max to prevent stuck pixels.
    max_valid_dt = 0.5; % [seconds]
    valid_mask = (dt_events < max_valid_dt) & (dt_events > 0);
    
    % 4. Update Tau Map
    %    Rule: tau = dt * K
    %    We only update tau where we have a valid flow measurement.
    %    For invalid pixels, we keep the old tau (or could reset to default).
    
    new_tau_values = dt_events(valid_mask) * sits_scale;
    
    % Update the global tau map at the valid indices
    tau_map(linear_idx(valid_mask)) = new_tau_values;
    
    % Optional: Enforce a minimum tau to prevent instant disappearance
    % tau_map(tau_map < 1e-4) = 1e-4; 

    %% 2. UPDATE TIMESTAMPS & POLARITY
    % Standard SAE updates
    last_t_map(linear_idx) = t;
    
    p_signed = double(p);
    p_signed(p_signed == 0) = -1;
    last_p_map(linear_idx) = p_signed;
    
    %% 3. COMPUTE VISUALIZATION
    % Formula: S = P * exp( - (t_now - t_last) / tau(x,y) )
    
    % Age of the event at every pixel
    age_map = t_current - last_t_map;
    age_map(age_map < 0) = 0; % Sanity check
    
    % Apply spatially varying decay
    % Note: Pixels with no events (tau=inf or age=inf) will naturally handle 
    % themselves if initialized correctly, but explicit handling is safer.
    
    decay_factor = exp( -age_map ./ tau_map );
    
    % Handle potential NaNs (0/0) or Infs
    decay_factor(isnan(decay_factor)) = 0;
    
    visual_surface = last_p_map .* decay_factor;

end
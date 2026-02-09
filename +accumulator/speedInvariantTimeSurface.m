function [visual_surface, last_t_map, tau_map, last_p_map] = ...
    speedInvariantTimeSurface(last_t_map, tau_map, last_p_map, ...
    x, y, t, p, imgSz, t_current, sits_scale)

    % Prepare indices
    linear_idx = sub2ind(imgSz, x, y);
    num_events = length(t);
    
    % Parameter: Max dt to consider part of a "flow" (0.5s)
    max_valid_dt = 0.5;
    
    % --- STEP 1: SEQUENTIAL STATE UPDATE ---
    % We must loop to ensure that if a pixel appears twice in the batch,
    % the second appearance references the updated time of the first.
    
    for i = 1:num_events
        idx = linear_idx(i);
        t_now = t(i);
        
        % 1. Get previous time from the map (which might have just been updated!)
        t_prev = last_t_map(idx);
        
        % 2. Update Timestamp immediately
        last_t_map(idx) = t_now;
        
        % 3. Update Polarity
        p_val = double(p(i));
        if p_val == 0, p_val = -1; end
        last_p_map(idx) = p_val;
        
        % 4. Calculate dt and Update Tau
        dt = t_now - t_prev;
        
        % Only update tau if this looks like valid motion (positive and recent)
        if dt > 0 && dt < max_valid_dt
            % SITS Rule: tau scales with dt
            new_tau = dt * sits_scale;
            tau_map(idx) = new_tau;
        end
    end
    
    % --- STEP 2: VISUALIZATION (Vectorized) ---
    % This part is safe to vectorize because it relies on the final state maps
    
    % Calculate age of the last event at every pixel
    age_map = t_current - last_t_map;
    age_map(age_map < 0) = 0; 
    
    % Safety: Ensure no division by zero or tiny taus causing NaNs
    % Clamp tau to a small minimum (e.g., 1 microsecond)
    safe_tau_map = max(tau_map, 1e-6);
    
    % Apply Decay
    decay_factor = exp( -age_map ./ safe_tau_map );
    
    visual_surface = last_p_map .* decay_factor;
    
end
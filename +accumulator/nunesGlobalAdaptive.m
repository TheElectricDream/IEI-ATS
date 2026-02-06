function [nunes_activity_map, G_activity, G_last_t] = nunesGlobalAdaptive(sorted_t, group_ends, ...
    unique_idx, pos, G_activity, G_last_t)

    % Iterate over each active pixel
    for k = 1:length(unique_idx)
        idx = unique_idx(k);
        
        % Get the sequence of timestamps for this specific pixel
        pixel_events_t = sorted_t(pos(k):group_ends(k));
        
        % Retrieve previous state from Globals
        a_curr = G_activity(idx);
        t_curr = G_last_t(idx);
        
        % Apply recursive update for every event in the sequence
        % Formula: a_i = beta_i * a_{i-1} + 1
        %          beta_i = 1 / (1 + a_{i-1} * (t_i - t_{i-1}))
        for j = 1:length(pixel_events_t)
            t_new = pixel_events_t(j);
            dt = t_new - t_curr;
            
            % Eq (15): Adaptive Decay
            % Protect against potential dt=0 if events are simultaneous
            beta = 1.0 / (1.0 + a_curr * dt);
            
            % Eq (10): Event Activity Update
            a_curr = beta * a_curr + 1;
            
            % Update time
            t_curr = t_new;
        end
        
        % Store updated state back to Globals
        G_activity(idx) = a_curr;
        G_last_t(idx) = t_curr;
    end
    
    % For visualization: Use the updated global activity map
    nunes_activity_map = G_activity;

end
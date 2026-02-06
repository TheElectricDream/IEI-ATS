function [time_surface_map, tau_filtered, decayed_surface] = ...
    localAdaptiveTimeSurface(t_mean, last_event_timestamp,...
    filtered_coherence_map, time_surface_map_prev, frameIndex, surface_k_tau, ...
    surface_tau_max, surface_tau_min, recency_T, filter_mask, norm_trace_map)

    t_mean(isnan(t_mean)) = 0; 
    filtered_coherence_map(isnan(filtered_coherence_map)) = 0;

    % Calculate the candidate decay time
    candidate_tau_n = surface_k_tau.*max(t_mean, eps);

    % Bound the candidate map
    candidate_tau_n(candidate_tau_n>surface_tau_max) = surface_tau_max;
    candidate_tau_n(candidate_tau_n<surface_tau_min) = surface_tau_min;

    % Calculate recency 
    recency_function = zeros(size(candidate_tau_n));
    last_event_timestamp(isnan(last_event_timestamp)) = 0;
    mask_recency = (last_event_timestamp ~= 0);
    recency_function(mask_recency) = ...
        exp(-last_event_timestamp(mask_recency) ./ recency_T);
    recency_weighted_tau = candidate_tau_n .* recency_function + ...
        surface_tau_min .* (1 - recency_function);

    % Accumulate the time surface
    if frameIndex == 1
        time_surface_map = filtered_coherence_map; % Initialize with input, not time constants

        % % Initialize previous maps
        % time_surface_map_prev = time_surface_map;
        % tau_map_prev = recency_weighted_tau;

        % For visualization
        tau_filtered = recency_weighted_tau;
        decayed_surface = time_surface_map_prev .* exp(-0.3 ./ tau_filtered); 
    else

        tau_current = recency_weighted_tau;
        tau_current(isnan(tau_current)) = surface_tau_min; 

        tau_filtered = imgaussfilt(tau_current, 3.0, "FilterSize", 3); 

        % Logic correction from previous turn (Accumulate, don't overwrite)
        decayed_surface = time_surface_map_prev .* exp(-0.3 ./ tau_filtered);

        % 3. FINAL ACCUMULATION: Add new clean input to decayed history
        time_surface_map = time_surface_map_prev*0.8 + (1-0.8).*imgaussfilt(filter_mask.*1, 3.0, "FilterSize", 9).*norm_trace_map;

    end

end
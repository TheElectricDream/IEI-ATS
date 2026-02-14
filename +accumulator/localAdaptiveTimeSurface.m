function [normalized_output_frame, time_surface_map, tau_filtered, decayed_surface] = ...
    localAdaptiveTimeSurface(t_mean, last_event_timestamp, ...
    time_surface_map_prev, alts_params, filter_mask, polarity_map)

    % Extract parameters
    surface_tau_max      = alts_params.surface_tau_max;
    surface_tau_min      = alts_params.surface_tau_min;
    dt                   = alts_params.dt;
    recency_filter_sigma = alts_params.recency_filter_sigma;
    recency_filter_size  = alts_params.recency_filter_size;

    % Set any NaN values to 0 for computation
    t_mean(isnan(t_mean)) = 0; 

    % Calculate the candidate decay time
    candidate_tau_n = max(t_mean, eps);

    % Bound the candidate map
    candidate_tau_n(candidate_tau_n>surface_tau_max) = surface_tau_max;
    candidate_tau_n(candidate_tau_n<surface_tau_min) = surface_tau_min;

    % Calculate recency 
    recency_function = zeros(size(candidate_tau_n));
    mask_recency = (last_event_timestamp ~= 0);
    recency_function(mask_recency) = ...
        exp(-last_event_timestamp(mask_recency));
    recency_weighted_tau = candidate_tau_n .* recency_function + ...
        surface_tau_max .* (1 - recency_function);

    % Smooth the tau map (spatial consistency)
    tau_current = recency_weighted_tau;
    tau_current(isnan(tau_current)) = surface_tau_min; 
    tau_filtered = imgaussfilt(tau_current, recency_filter_sigma,...
        "FilterSize", recency_filter_size); 

    % A. Apply Adaptive Decay
    decay_factor = exp(-dt ./ tau_current);
    decayed_surface = time_surface_map_prev .* decay_factor;

    % B. Add Weighted Polarity Input
    % We use the polarity_map instead of norm_trace_map.
    % This adds +1 (ON) or -1 (OFF) to the surface.
    
    % Apply the filter_mask to the polarity map to remove noise
    masked_input = polarity_map .* filter_mask;
    
    % ACCUMULATE
    time_surface_map = masked_input+decayed_surface;

    normalized_output_frame = normalizeSurface(time_surface_map);

    % % Hot Pixel Map
    % normalized_output_frame(normalized_output_frame==1) = 0;

end

function norm_S = normalizeSurface(S)

    % Clip values to a reasonable contrast integration range    
    % Robust Auto-scale (Ignore outliers)
    mean_val = mean(S(:));
    std_val = std(S(:));
    min_v = mean_val - 3*std_val;
    max_v = mean_val + 3*std_val;
    S_clipped = max(min(S, max_v), min_v);
    norm_S = (S_clipped - min_v) / (max_v - min_v);

end
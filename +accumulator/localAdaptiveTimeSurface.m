function [normalized_output_frame, time_surface_map, tau_filtered, decayed_surface] = ...
    localAdaptiveTimeSurface(t_mean, time_surface_map_prev, alts_params,...
    filter_mask, polarity_map)

    % Extract parameters
    surface_tau_max      = alts_params.surface_tau_max;
    surface_tau_min      = alts_params.surface_tau_min;
    surface_tau_release  = alts_params.surface_tau_release;
    dt                   = alts_params.dt;
    recency_filter_sigma = alts_params.recency_filter_sigma;
    recency_filter_size  = alts_params.recency_filter_size;

    % Set any NaN values to 0 for computation
    t_mean(isnan(t_mean)) = 0; 

    % Calculate the candidate decay time
    tau_active = max(t_mean, eps);
    
    % Sigmoid-like mapping: small IEI -> tau_min, large IEI -> tau_max
    tau_active = process.sigmoidRemap(tau_active, surface_tau_min, surface_tau_max);

    % Smooth the tau map (spatial consistency)
    % tau_active = recency_weighted_tau;
    % tau_active(isnan(tau_active)) = surface_tau_min; 
    tau_filtered = imgaussfilt(tau_active, recency_filter_sigma,...
        "FilterSize", recency_filter_size); 

    % Build a current-activity indicator: 1 where events arrive, 
    % 0 where idle
    activity_indicator = single(polarity_map ~= 0);

    % Smooth it spatially so the transition isn't pixel-sharp
    activity_blurred = imgaussfilt(activity_indicator, recency_filter_sigma, ...
        "FilterSize", recency_filter_size);
    activity_blurred = min(activity_blurred ./ max(activity_blurred(:) + eps), 1.0);

    % Active pixels keep tau_active, idle pixels with residual get tau_release
    tau_effective = tau_active .* activity_blurred + ...
                    surface_tau_release .* (1 - activity_blurred); 

    % Apply Adaptive Decay
    decay_factor = exp(-dt ./ tau_effective);
    decayed_surface = time_surface_map_prev .* decay_factor;

    % Add Weighted Polarity Input
    % We use the polarity_map instead of norm_trace_map.
    % This adds +1 (ON) or -1 (OFF) to the surface.
    
    % Apply the filter_mask to the polarity map to remove noise
    masked_input = polarity_map .* filter_mask;
    
    % Accumulate
    time_surface_map = masked_input + decayed_surface;

    % Normalize the output frame
    normalized_output_frame = normalizeSurface(time_surface_map);

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
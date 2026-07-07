function [lats_outputs] = ...
    localAdaptiveTimeSurface(t_mean, time_surface_map_prev, ...
    alts_params, filter_mask, polarity_map, counts)
    % LOCALADAPTIVETIMESURFACE  IEI-ATS adaptive local time surface.
    
    % Extract parameters
    surface_tau_release  = alts_params.surface_tau_release;
    dt                   = alts_params.dt;
    div_norm_exp         = alts_params.div_norm_exp;
    symmetric_tone_scale = alts_params.symmetric_tone_scale;
    sigma_base           = alts_params.sigma_base;
    sigma_outlier        = alts_params.sigma_outlier;

    % Sanitize IEI input
    t_mean(isnan(t_mean)) = 0;

    % ----------------------------------------------------------------
    % 1. Time constant mapping (identity: tau = IEI)
    % ----------------------------------------------------------------
    tau_filtered = max(t_mean, eps);

    % ----------------------------------------------------------------
    % 2. Activity indicator (morphological dilation)
    % ----------------------------------------------------------------
    % Apply coherence mask BEFORE computing activity so the envelope
    % sees only events that survive filtering. Without this, filtered
    % pixels appear "active" and get tau_active instead of tau_release.
    masked_input = polarity_map .* filter_mask;
    
    % Apply outlier rejection to filter out any points that are way too big
    [masked_input, ~, outlier_mask] =...
        stats.rejectPolarityOutliers(masked_input, masked_input, sigma_outlier);

    activity_indicator = single(masked_input ~= 0);
    activity_indicator(isnan(activity_indicator)) = 0;

    % ----------------------------------------------------------------
    % 3. Asymmetric attack-release envelope
    % ----------------------------------------------------------------
    tau_effective = tau_filtered .* activity_indicator + ...
                    surface_tau_release .* (1 - activity_indicator);

    % ----------------------------------------------------------------
    % 4. Coupled decay factor (Inverted Gain)
    % ----------------------------------------------------------------
    gamma_min = 0.15; 
    
    % Map the exponential decay so it bottoms out at gamma_min instead of 0
    raw_decay = exp(-dt ./ tau_effective);
    adaptive_gains = gamma_min + (1 - gamma_min) .* raw_decay;
    
    % ----------------------------------------------------------------
    % 5. Amplitude-Modulated Leaky Integrator
    % ----------------------------------------------------------------

    clamped_input = tanh(masked_input);

    time_surface_map_raw = (adaptive_gains .* clamped_input) ...
        + (adaptive_gains .* time_surface_map_prev);

    % ----------------------------------------------------------------
    % 6. Divisive normalization (Carandini-Heeger)
    % ----------------------------------------------------------------
    magnitude = abs(time_surface_map_raw);
    time_surface_map = zeros(size(masked_input));

    % Sanitize auxiliary data arrays to prevent denominator corruption
    if any(outlier_mask(:))
        % Calculate safe fallbacks (ignoring the outliers and empty space)
        valid_counts = counts(~outlier_mask & counts > 0);
        if ~isempty(valid_counts)
            safe_count = median(valid_counts);
            counts(outlier_mask) = safe_count;
        end

        % If magnitude is tracking separately prior to the surface update:
        valid_mags = magnitude(~outlier_mask & magnitude > 0);
        if ~isempty(valid_mags)
            safe_mag = median(valid_mags);
            magnitude(outlier_mask) = safe_mag;
        end
    end

    % Calculate the time surface map with divisive normalization
    time_surface_map((magnitude > 0)) = ...
        time_surface_map_raw((magnitude > 0)) ./ ...
        (sigma_base + counts((magnitude > 0)) .^ div_norm_exp);

    % Bipolar surface mapped to [0, 1] via symmetric sigmoid.
    % Zero maps to 0.5 (mid-gray).
    normalized_output_frame = ...
        process.symmetricToneMappingNorm(time_surface_map, ...
        symmetric_tone_scale);

    % Return LATS outputs
    lats_outputs.normalized_output_frame = normalized_output_frame;
    lats_outputs.time_surface_map = time_surface_map;
    lats_outputs.magnitude = magnitude;
    lats_outputs.adaptive_gains = adaptive_gains;
    lats_outputs.time_surface_map_raw = time_surface_map_raw;
    lats_outputs.clamped_input = clamped_input;
    lats_outputs.tau_effective = tau_effective;
    lats_outputs.tau_filtered = tau_filtered;

end


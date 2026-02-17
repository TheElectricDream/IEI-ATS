function [normalized_output_frame, time_surface_map_raw, tau_filtered, decayed_surface, adaptive_gains] = ...
    localAdaptiveTimeSurface(t_mean, time_surface_map_prev, alts_params,...
    filter_mask, polarity_map, counts)

    % Extract parameters
    surface_tau_max      = alts_params.surface_tau_max;
    surface_tau_min      = alts_params.surface_tau_min;
    surface_tau_release  = alts_params.surface_tau_release;
    dt                   = alts_params.dt;
    recency_filter_sigma = alts_params.recency_filter_sigma;
    recency_filter_size  = alts_params.recency_filter_size;
    pool_sigma = recency_filter_sigma;
    pool_filter_size = recency_filter_size;

    % Set any NaN values to 0 for computation
    t_mean(isnan(t_mean)) = 0; 

    % Calculate the candidate decay time
    tau_active = max(t_mean, eps);
    
    % Sigmoid-like mapping: small IEI -> tau_min, large IEI -> tau_max
    tau_active = process.sigmoidRemap(tau_active, surface_tau_min, surface_tau_max);
    tau_filtered = imgaussfilt(tau_active, recency_filter_sigma,...
        "FilterSize", recency_filter_size); 

    % Apply the filter_mask BEFORE computing activity so the
    % attack-release envelope sees only events that survive filtering.
    % Without this, filtered pixels appear "active" and get the slow
    % tau_active decay instead of the fast tau_release, causing trails.
    masked_input = polarity_map .* filter_mask;

    % Build activity indicator from the FILTERED input
    activity_indicator = single(masked_input ~= 0);
    activity_indicator(isnan(activity_indicator))=0;

    % Smooth it spatially so the transition isn't pixel-sharp
    % activity_blurred = imgaussfilt(activity_indicator, recency_filter_sigma, ...
    %     "FilterSize", recency_filter_size);
    % activity_blurred = min(activity_blurred ./ max(activity_blurred(:) + eps), 1.0);
    se = strel('disk', 1);  % fills 1-pixel gaps, no boundary smear
    activity_blurred = imdilate(activity_indicator, se);

    % Active pixels keep tau_active, idle pixels with residual get tau_release
    tau_effective = tau_active .* activity_blurred + ...
                    surface_tau_release .* (1 - activity_blurred);

    % % Apply Adaptive Decay
    decay_factor = exp(-dt./tau_effective);
    decayed_surface = time_surface_map_prev .* decay_factor;
    
    % % Accumulate
    % time_surface_map = masked_input + decayed_surface;

    % Compute per-pixel blending coefficient (coupled gain + decay)
    adaptive_gains = 1 - exp(-dt ./ tau_effective);

    adaptive_gains = imgaussfilt(adaptive_gains, recency_filter_sigma, ...
         "FilterSize", recency_filter_size);
   
    % Remap the results
    %masked_input_remapped = process.linearRemap(masked_input,-0.5, 0.5)+0.5;
    
    % % EMA Update
    % time_surface_map = adaptive_gains .* masked_input + (1 - adaptive_gains) .* time_surface_map_prev;

    % 1. EMA update (unchanged — controls temporal tracking rate)
    time_surface_map_raw = adaptive_gains .* masked_input ...
                         + (1 - adaptive_gains) .* time_surface_map_prev;
    
    % Decompose into magnitude and sign
    magnitude = abs(time_surface_map_raw);

    % 2. Build the normalization pool (smoothed local energy)
    %    abs() because polarity_map is signed
    activity_pool = imgaussfilt(magnitude, recency_filter_sigma, 'FilterSize', recency_filter_size);

    % Build normalization pool from EVENT COUNTS (unsigned activity)
    counts_smooth = imgaussfilt(single(counts), pool_sigma, 'FilterSize', pool_filter_size);
    
    % 3. Divisive normalization
    %    sigma controls the crossover: regions with activity >> sigma get
    %    compressed toward 1/activity_pool; regions with activity << sigma
    %    pass through nearly unchanged.
    
    time_surface_map = zeros(size(masked_input));
    good_mask = (abs(masked_input)>0);
    sigma = median(activity_pool(abs(activity_pool)>0));  
    n = 1.0;        

    % Divisive normalization: signed surface / unsigned activity
    time_surface_map(good_mask) = time_surface_map_raw(good_mask) ./ (sigma + counts_smooth(good_mask) .^ n);
    
    % Simple outlier rejection
    mean_value_pos = mean(time_surface_map(time_surface_map>0));
    mean_value_neg = mean(time_surface_map(time_surface_map<0));
    std_value_pos = std(time_surface_map(time_surface_map>0));
    std_value_neg = std(time_surface_map(time_surface_map<0));

    % Reject points which are still 3\sigma outside the mean
    pos_threshold = mean_value_pos + 4*std_value_pos;
    neg_threshold = mean_value_neg - 4*std_value_neg;

    time_surface_map(time_surface_map>pos_threshold)= median(time_surface_map(time_surface_map>0));
    time_surface_map(time_surface_map<neg_threshold)= median(time_surface_map(time_surface_map<0));

    % Normalize the output frame
    normalized_output_frame = normalizeSurface(time_surface_map, 8, 1.0);
    % remapped_time_surface = process.linearRemap(time_surface_map, -0.5, 0.5)+0.5;
    % shift_value = mean(remapped_time_surface((abs(masked_input)==0)))-0.5;
    % normalized_output_frame = remapped_time_surface-shift_value;
    % bad_pixels = (normalized_output_frame>0.78 | normalized_output_frame<0.2);
    % time_surface_map_raw(bad_pixels) = sigma;
    % 
    % normalized_output_frame(bad_pixels)=0.5;
    % normalized_output_frame(bad_pixels)=0.5;

end

% function norm_S = normalizeSurface(S)
% 
%     % Clip values to a reasonable contrast integration range    
%     % Robust Auto-scale (Ignore outliers)
%     S(isnan(S))=0;
%     mean_val = mean(S(:));
%     std_val = std(S(:));
%     min_v = mean_val - 3*std_val;
%     max_v = mean_val + 3*std_val;
%     S_clipped = max(min(S, max_v), min_v);
%     norm_S = (S_clipped - min_v) / (max_v - min_v);
% 
% end

% function norm_S = normalizeSurface(S)
%     % 1. Log-compress to reduce extreme dynamic range before CLAHE
%     sign_S = sign(S);
%     abs_S  = abs(S);
%     compressed = sign_S .* log1p(abs_S);
% 
%     % 2. Rescale to [0, 1] for adapthisteq input
%     c_min = min(compressed(:));
%     c_max = max(compressed(:));
%     if c_max - c_min < eps
%         norm_S = 0.5 * ones(size(S));
%         return;
%     end
%     prescaled = (compressed - c_min) / (c_max - c_min);
% 
%     % 3. CLAHE: locally adaptive equalization
%     %    - NumTiles controls the spatial granularity of adaptation
%     %    - ClipLimit controls how much contrast enhancement is allowed
%     %      (lower = more uniform, higher = more local contrast)
%     norm_S = adapthisteq(prescaled, ...
%         'NumTiles',  [16 16], ...
%         'ClipLimit', 0.02, ...
%         'Distribution', 'uniform');
% end

% function norm_S = normalizeSurface(S, scale)
%     % Fixed symmetric tone mapping via hyperbolic tangent.
%     % Maps S -> [0, 1] with midpoint at 0.5 (zero surface = gray).
%     %
%     %   scale controls contrast:
%     %     - Larger scale  = softer curve, more headroom for extremes
%     %     - Smaller scale = steeper curve, more contrast in quiet regions
%     %
%     % The tanh function is a standard sigmoidal tone-mapping operator;
%     % see Reinhard et al. (2002), "Photographic Tone Reproduction for
%     % Digital Images," ACM SIGGRAPH, Eq. 4 and discussion of sigmoid
%     % compression for high dynamic range imagery.
% 
%     if nargin < 2
%         scale = 3.0;
%     end
% 
%     norm_S = 0.5 + 0.5 * tanh(S / scale);
% end

function norm_S = normalizeSurface(S, scale, detail_boost)

    if nargin < 2, scale = 3.0; end
    if nargin < 3, detail_boost = 1.0; end

    S = double(S);

    % 1. Edge-aware decomposition
    base = imbilatfilt(S, 2.0, 8);
    detail = S - base;

    % 2. Compress only the base
    base_compressed = tanh(base / scale);

    % 3. Recombine — detail is NOT divided by scale
    %    detail_boost directly controls local contrast strength
    recombined = base_compressed + detail_boost * detail;

    % 4. Fixed mapping to [0, 1]
    norm_S = 0.5 + 0.5 * tanh(recombined);

end
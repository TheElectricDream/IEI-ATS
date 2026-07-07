function [agd] = initializeNunesAGD(img_size, frame_total)

    % ADAPTIVE GLOBAL DECAY (AGD) TIME SURFACE PARAMETERS
    % REF: https://ieeexplore.ieee.org/document/10205486/

    % Tuning Guide:
    %   The Adaptive Global Decay Time Surface constructs a
    %   spatiotemporal representation in which past event activity decays
    %   at a rate that automatically adapts to the global scene dynamics.
    %   Rather than applying a fixed exponential time constant (which
    %   must be manually tuned and cannot simultaneously handle both
    %   fast and slow motion), the decay rate is governed by the event
    %   activity: a running estimate of the global event rate. When the
    %   scene produces events rapidly (fast motion or dense texture),
    %   the activity is high, the decay is steep, and only very recent
    %   events contribute to the surface. When the scene is slow or
    %   sparse, the activity is low, the decay is gentle, and events
    %   persist longer. The adaptive decay for each pixel takes the
    %   form beta = 1 / (1 + a * dt), where a is the current event
    %   activity and dt is the elapsed time since the pixel was last
    %   updated.
    %
    %   alpha (activity smoothing factor): Controls how quickly the
    %   event activity estimate responds to changes in the global event
    %   rate. Increasing alpha makes the activity estimate more
    %   reactive, causing the decay rate to track rapid fluctuations in
    %   scene dynamics more closely. This suits scenes with abrupt
    %   transitions between fast and slow motion. Decreasing alpha
    %   smooths the activity estimate over a longer history, producing
    %   more stable but slower-adapting decay behaviour. This suits
    %   scenes with gradually varying dynamics where stability of the
    %   representation is preferred over responsiveness.
    %
    %   K (activity scaling constant): Scales the event activity
    %   relative to the sensor resolution and expected event rates.
    %   Increasing K amplifies the effective activity, which steepens
    %   the decay and causes past events to fade more rapidly. This may
    %   be necessary for high-resolution sensors or scenes that produce
    %   very high event rates. Decreasing K attenuates the effective
    %   activity, producing a gentler decay that retains past events
    %   longer. This suits lower-resolution sensors or scenes with
    %   moderate event rates. As a starting point, K can be set on the
    %   order of the total number of pixels in the sensor array, and
    %   then adjusted based on whether the resulting surface appears
    %   overly decayed (reduce K) or overly persistent (increase K).
    
    agd.agd_params.alpha            = 0.001;   
    agd.agd_params.K                = 50000.0;  
    agd.agd_surface                 = zeros(img_size); 
    agd.agd_state.last_t_map        = zeros(img_size);
    agd.agd_state.activity          = 0; 
    agd.agd_state.last_update_time  = 0; 
    agd.agd_activity_store          = zeros(length(frame_total),1);
    agd.normalized_output_frame     = zeros(img_size);

end
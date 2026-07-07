function [hots] = initializeLagorceHOTS(img_size)

    % TIME SURFACE (HOTS) PARAMETERS
    % REF: http://ieeexplore.ieee.org/document/7508476/

    % Tuning Guide:
    %   The Time Surface accumulator constructs a spatiotemporal
    %   representation of recent event activity by maintaining a map of
    %   the most recent event timestamp at each pixel. When a new event
    %   arrives, the time surface at its location is computed by applying
    %   an exponential decay kernel to the elapsed time since the last
    %   event at each pixel: S(u) = exp(-(t_now - t_last(u)) / tau).
    %   The resulting surface takes values between 0 and 1, where values
    %   near 1 indicate very recent activity and values near 0 indicate
    %   stale or absent activity. This provides a continuous, polarity-
    %   independent descriptor of the local spatiotemporal context around
    %   each event.
    %
    %   ts_time_constant (exponential decay time constant, tau): Controls
    %   how rapidly past event activity fades in the time surface, in
    %   seconds. Increasing tau causes past events to persist longer in
    %   the representation, producing smoother and more temporally
    %   extended surfaces. This is appropriate for slow-moving objects or
    %   low event rates where the temporal separation between related
    %   events may be large. Decreasing tau causes past activity to decay
    %   more quickly, producing surfaces that are sharply peaked around
    %   only the most recent events. This suits fast motion or high event
    %   rates where only the immediate temporal context is relevant.
    %   Lagorce et al. report typical first-layer values of 20-50 ms,
    %   with subsequent layers in a hierarchical architecture using
    %   progressively larger time constants (scaled by a factor K_tau
    %   per layer) to integrate activity over longer temporal windows.
    
    hots.ts_time_constant = 2;  
    hots.ts_t_map = -inf(img_size);
    hots.normalized_output_frame = zeros(img_size);

end
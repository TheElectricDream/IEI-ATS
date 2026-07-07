function [mcf] = initializeWangMCF(img_size, frame_total)

    % MOTION CONSISTENCY FILTER
    % REF: https://ieeexplore.ieee.org/document/8953966/

    % Tuning Guide:
    %   The Motion Consistency Filter classifies each event as signal or
    %   noise by fitting a plane to neighbouring events in the
    %   spatiotemporal (x, y, t) domain via least squares. The slope of
    %   the fitted plane yields a local optical flow velocity estimate.
    %   Events produced by genuine object motion should lie on a plane
    %   whose velocity is non-zero and physically plausible, whereas noise
    %   events lack the spatial coherence to produce a consistent plane.
    %
    %   delta_t (temporal neighbourhood half-width): Defines how far
    %   forward and backward in time, in seconds, the filter searches for
    %   neighbouring events to include in the plane fit. Increasing this
    %   value draws on a longer temporal window, which is appropriate for
    %   slow-moving objects or low event rates where nearby supporting
    %   events may be temporally spread out. Decreasing it restricts the
    %   fit to very recent events, suiting fast motion where temporally
    %   distant events are unlikely to belong to the same local motion.
    %
    %   V_max (maximum admissible velocity): The upper bound on the
    %   optical flow speed, in pixels per second, that the filter will
    %   accept as plausible motion. Events whose fitted plane yields a
    %   velocity magnitude exceeding this value are rejected as noise.
    %   Raising V_max permits faster apparent motion to pass through,
    %   which may be necessary for high-speed scenes or close-range
    %   objects but increases the risk of admitting noise that
    %   coincidentally forms a steep plane. Lowering it enforces stricter
    %   motion plausibility at the risk of discarding valid events from
    %   fast-moving edges. Note that events with a fitted velocity of
    %   exactly zero are also rejected, as a stationary plane implies no
    %   genuine motion.
    %
    %   min_neighbours (minimum support count): The minimum number of
    %   neighbouring events required within the spatiotemporal
    %   neighbourhood (3x3 spatial, +/- delta_t temporal) before
    %   attempting a plane fit. Increasing this value demands more
    %   evidence before evaluating an event, improving the reliability of
    %   the fitted velocity at the cost of rejecting isolated valid
    %   events in sparse regions. Decreasing it allows plane fits with
    %   very few supporting events, which may retain more signal in
    %   sparse scenes but produces less stable velocity estimates.
    
    mcf.mcf_params.delta_t         = 0.1;   
    mcf.mcf_params.V_max           = 500;    
    mcf.mcf_params.min_neighbours  = 2;       
    mcf.mcf_eventBuffer.x          = [];
    mcf.mcf_eventBuffer.y          = [];
    mcf.mcf_eventBuffer.t          = [];
    mcf.mcf_n_passed_store         = zeros(frame_total, 1);
    mcf.mcf_n_total_store          = zeros(frame_total, 1);

end
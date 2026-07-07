function [mets] = initializeXuMETS(img_size)

    % MOTION-ENCODED TIME-SURFACE (METS) PARAMETERS
    % REF: https://link.springer.com/10.1007/s11263-025-02379-6

    % Tuning Guide:
    %   The Motion-Encoded Time-Surface (METS) constructs an event
    %   representation in which the exponential decay rate at each pixel
    %   is dynamically modulated by the locally estimated edge velocity,
    %   rendering the surface invariant to motion speed. When an event
    %   activates a pixel, METS estimates the instantaneous velocity of
    %   the scene edge at that location by extracting a temporal sequence
    %   from a polarity-separated timestamp memory within a local
    %   observation window. The pixel is set to 1 and then decays
    %   exponentially at a rate tied to the estimated velocity: the
    %   surface value drops by a factor of e for every d pixels of
    %   subsequent edge displacement. Once the edge has moved d_th
    %   pixels (inferred from elapsed time and velocity), the pixel is
    %   reset to zero. Because the decay is expressed in units of edge
    %   displacement rather than absolute time, the resulting surface
    %   profile is identical regardless of whether the edge is moving
    %   slowly or rapidly. The authors report that parameter settings
    %   remain fixed across all scenes and speeds with no tuning
    %   required.
    %
    %   R (observation window half-width): Defines the spatial extent
    %   of the square window used to estimate local edge velocity. The
    %   full window is (2R+1) x (2R+1) pixels. Increasing R draws on
    %   a larger spatial neighbourhood for the velocity estimate,
    %   improving robustness in noisy or texturally complex scenes but
    %   increasing per-event computation. Decreasing R reduces the
    %   spatial context, which may suffice for clean data with simple
    %   edge structures and lowers computational cost.
    %
    %   n (velocity estimation depth): The number of edge-width steps
    %   sampled from the temporal memory to estimate instantaneous edge
    %   velocity. The temporal sequence extracted has length
    %   n*(2R+1) + 1 entries, and its span is interpreted as the time
    %   for the edge to traverse n pixels. Increasing n integrates
    %   velocity over a longer motion history, producing a smoother
    %   estimate that is more robust to timestamp noise. Decreasing it
    %   makes the estimate more responsive to instantaneous changes in
    %   motion but more susceptible to noise.
    %
    %   d (decay step): The edge displacement, in pixels, over which
    %   the surface value decays by a factor of e (approximately 2.72).
    %   Increasing d produces a gentler gradient behind the moving
    %   edge, yielding a wider and smoother trailing slope on the
    %   surface. Decreasing it steepens the gradient, concentrating
    %   high surface values closer to the most recent edge position.
    %
    %   d_th (decay distance threshold): The edge displacement, in
    %   pixels, at which a pixel is reset to zero. This controls the
    %   effective thickness of edges in the representation and the
    %   signal-to-noise ratio of the surface. Increasing d_th allows
    %   past activity to persist over a longer trailing distance,
    %   which may improve alignment in downstream tasks but risks
    %   retaining stale edges. Decreasing it produces thinner, crisper
    %   edges at the cost of a narrower temporal support region.
    
    mets.mets_params.R    = 4;   
    mets.mets_params.n    = 3;    
    mets.mets_params.d    = 5;   
    mets.mets_params.d_th = 8;   
    mets.mets_state.t_last_pos = zeros(img_size);   
    mets.mets_state.t_last_neg = zeros(img_size);  
    mets.mets_state.t_last_any = zeros(img_size);   
    mets.mets_state.p_last     = zeros(img_size); 
    mets.mets_surface = zeros(img_size);
    mets.normalized_output_frame = zeros(img_size);

end
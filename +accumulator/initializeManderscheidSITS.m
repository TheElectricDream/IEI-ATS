function [sits] = initializeManderscheidSITS(img_size)

    % SPEED INVARIANT TIME SURFACE (SITS) PARAMETERS
    % REF: https://arxiv.org/pdf/1903.11332

    % Tuning Guide:
    %   The Speed Invariant Time Surface (SITS) constructs a
    %   spatiotemporal representation of recent event activity that is
    %   independent of the apparent speed of moving objects. Unlike the
    %   standard exponential-decay time surface, which produces steeper
    %   or shallower gradients depending on how fast an edge traverses
    %   the pixel array, the SITS maintains an identical profile
    %   regardless of speed. When a new event arrives at pixel (x, y)
    %   with polarity p, all neighbouring pixels within a (2R+1) x (2R+1)
    %   window whose stored value is greater than or equal to the current
    %   pixel's value are decremented by one. The event pixel is then set
    %   to (2R+1)^2. The resulting surface takes integer values between 0
    %   and (2R+1)^2, with the highest value at the most recent event and
    %   a fixed-slope ramp behind a moving edge whose length is exactly R
    %   pixels, independent of the edge's speed.
    %
    %   sits_R (neighbourhood half-width): Defines the spatial extent of
    %   the update region around each incoming event. The full
    %   neighbourhood is (2R+1) x (2R+1) pixels. Increasing R produces
    %   a longer and smoother ramp behind moving edges, which captures
    %   more spatial context and may improve robustness for downstream
    %   tasks such as feature detection. However, larger values increase
    %   the per-event computation cost (which scales as (2R+1)^2) and
    %   may blur fine spatial structure. Decreasing R shortens the ramp,
    %   preserving sharper spatial detail at the cost of less context.
    %   Manderscheid et al. report using R = 6 in combination with an
    %   8 x 8 classifier input patch for corner detection.
    
    sits.sits_t_map = zeros(img_size);
    sits.sits_R = 3;
    sits.normalized_output_frame = zeros(img_size);

end
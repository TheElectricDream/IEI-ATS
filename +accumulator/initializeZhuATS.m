function [zhu] = initializeZhuATS(img_size)

    % ZHU ADAPTIVE TIME SURFACE (EVO-ATS) PARAMETERS
    % REF: https://ieeexplore.ieee.org/document/10342048/

    % Tuning Guide:
    %   The Zhu Adaptive Time Surface (EVO-ATS) addresses the whiteout
    %   and blackout problems of constant-decay time surfaces by
    %   computing a pixel-wise decay rate that adapts to both local
    %   motion speed and scene texture complexity. For each pixel, the
    %   decay rate is derived from the timestamps of a fixed number of
    %   neighbouring pixels selected by a sparse spatial pattern. The
    %   mean temporal gap between the current time and the n most recent
    %   neighbouring events is subtracted from an upper bound to yield
    %   the local decay rate: tau(x) = max(tau_u - mean(t - t_last,i),
    %   tau_l). In regions of fast motion or dense texture, neighbouring
    %   events are recent, the mean gap is small, and the decay rate
    %   remains close to tau_u, causing rapid decay that prevents
    %   overlapping trails. In regions of slow motion or sparse texture,
    %   the mean gap is large, the subtraction lowers the decay rate
    %   toward tau_l, and past activity persists longer to retain
    %   sufficient information.
    %
    %   tau_u (upper bound on decay rate): The maximum decay rate, in
    %   seconds, assigned to any pixel. This value is used when the
    %   surrounding neighbourhood has seen very recent activity,
    %   indicating fast motion or dense texture. Increasing tau_u
    %   allows the surface to retain activity longer even in high-
    %   activity regions, which may cause trail overlap. Decreasing it
    %   produces faster decay in active regions, yielding thinner and
    %   sharper edge trails.
    %
    %   tau_l (lower bound on decay rate): The minimum decay rate, in
    %   seconds, that prevents the surface from decaying too rapidly in
    %   any region. This floor ensures that even in very active areas,
    %   the exponential decay does not collapse to zero before the
    %   surface can be sampled. Increasing tau_l raises the minimum
    %   persistence of all pixels. Decreasing it permits more
    %   aggressive decay in the most active regions but risks losing
    %   information before it can be used.
    %
    %   n_neighbors (number of neighbouring timestamps): The number of
    %   most recent neighbouring events, selected from a sparse spatial
    %   pattern surrounding the target pixel, used to compute the mean
    %   temporal gap. Increasing this value draws on more spatial
    %   context, producing a smoother and more representative estimate
    %   of local activity at the cost of additional lookups. Decreasing
    %   it makes the decay rate more sensitive to individual
    %   neighbouring pixels, which may be appropriate for very fine
    %   spatial structure but is more susceptible to noise.
    %
    %   blur_sigma (Gaussian blur standard deviation): The standard
    %   deviation, in pixels, of a Gaussian filter applied to the
    %   completed ATS to smooth out sharp discontinuities between
    %   adjacent pixels with different decay rates. Increasing this
    %   value produces a smoother surface with gentler gradients,
    %   which may aid convergence in downstream optimisation tasks.
    %   Decreasing it preserves sharper spatial detail at the risk of
    %   gradient discontinuities.
    %
    %   median_sz (median filter kernel size): The side length, in
    %   pixels, of a median filter applied to the ATS after Gaussian
    %   smoothing. This filter suppresses impulsive noise from isolated
    %   hot pixels or spurious events. Increasing the kernel size
    %   provides stronger noise suppression but may erode fine edge
    %   detail. Decreasing it preserves detail at the cost of reduced
    %   noise robustness.
    
    zhu.zhu_state.t_last       = zeros(img_size);
    zhu.zhu_params.tau_u       = 0.8;    
    zhu.zhu_params.tau_l       = 0.001;  
    zhu.zhu_params.n_neighbors = 16;     
    zhu.zhu_params.blur_sigma  = 1.0;   
    zhu.zhu_params.median_sz   = 3;  
    zhu.zhu_surface            = zeros(img_size);
    zhu.zhu_tau_map            = zeros(img_size);
    zhu.normalized_output_frame = zeros(img_size);

end
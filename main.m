%% Loop control
% Set this to "true" to run this code in a loop across all available
% filters
isLooping = true;

% If-else logic for the loop
if isLooping == false
    
    % Clear out all variables
    clear;
    clc;
    close all;

    % Use buffered data
    useBuffer = false;
    
    % Select algorithms for filtering and accumulation
    % Set 'None' for filter selection to skip filtering entirely
    
    % Set filtering
    % Options: 'NONE', 'STC', 'BAF', 'EDF', 'STCC', 'MCF', 'COH'
    filterSelection = 'EDF';
    
    % Set accumulator
    % Options: 'HOTS', 'SITS', 'METS', 'IEI-ATS', 'AGD', 'EVO-ATS'
    accumulatorSelection = 'IEI-ATS';

    % Set dataset name
    %fileName = 'recording_20260127_145247.hdf5';  % Jack W. (LED Cont)
    %fileName = 'recording_20251029_131131.hdf5';  % EVOS - NOM - ROT
    fileName = 'recording_20251029_135047.hdf5';  % EVOS - SG - ROT
    %fileName = 'recording_20251029_134602.hdf5';  % EVOS - DARK - ROT
    %fileName = 'ZED-ROT-NOM.h5';
    %fileName = 'ZED-ROT-NOM-SLOWMO.h5';

else

    % Use buffered data
    useBuffer = false;

end

%% Define processing range
% Define start and end time to process [seconds]
t_start_process = 0; 
t_end_process   = 1000; 

%% Import events for inspection

% Set path to datasets
hdf5Path = ['/home/alexandercrain/Dropbox/Graduate Documents' ...
    '/Doctor of Philosophy/Thesis Research/Datasets/SPOT/HDF5/'];

% Set output video name
videoOutPath = fullfile('/home/alexandercrain/Videos/Research', ...
    sprintf('Normalized-Output-Fltr-%s-Acmtr-%s-%s.avi', ...
    filterSelection, accumulatorSelection, datestr(now,'yyyymmdd-HHMMSS'))); %#ok<TNOW1,DATST>

if useBuffer == false

    % Load the data (actual event data)
    tk = double(h5read([hdf5Path fileName], '/timestamp'));
    xk = single(h5read([hdf5Path fileName], '/x'));
    yk = single(h5read([hdf5Path fileName], '/y'));
    pk = single(h5read([hdf5Path fileName], '/polarity'));

    % % Load the data (v2e simulated data)
    % data = h5read([hdf5Path fileName], '/events');
    % 
    % % data is [4 x N] in MATLAB (HDF5 transposes row/column order)
    % tk  = double(data(1, :))';   % timestamp [s]
    % xk  = single(data(2, :))';   % x
    % yk  = single(data(3, :))';   % y
    % pk  = single(data(4, :))';   % polarity (0 or 1)

    % Convert time to seconds
    tk = (tk - tk(1))./1e6;
    
    % Convert to single data type to use less memory
    tk = single(tk);
    
    % Find indices within the valid range
    valid_idx = tk >= t_start_process & tk <= t_end_process;
    
    % Filter the data vectors
    tk = tk(valid_idx);
    xk = xk(valid_idx);
    yk = yk(valid_idx);
    pk = pk(valid_idx);
    
    % Shift time to start at 0 for the new window
    % This ensures your frame loop starts correctly at frame 1
    tk = tk - t_start_process; 
    
    % Clear unused variables for memory
    clearvars valid_idx;

    % Set the time interval to accumulate over
    t_total                     = max(tk);  % [s]

else

    % Create buffered event reader
    buf = process.EventBuffer(fullfile(hdf5Path, fileName), ...
        t_start_process, t_end_process);

    t_total                     = buf.t_total;

end


% Set the time interval to accumulate over
t_interval                  = 0.33;     % [s]
frame_total                 = floor(t_total / t_interval);

% Set the image size
imgSz                       = [640, 480]; 

%% Initialize Filter Parameters

% BACKGROUND ACTIVITY FILTER PARAMETERS (Delbruck 2008)
% -------------------------------------------------------------------------
% REF: https://www.zora.uzh.ch/handle/20.500.14742/41289
% -------------------------------------------------------------------------
% Tuning Guide:
%   The Background Activity Filter removes uncorrelated noise events
%   caused by thermal noise and junction leakage currents in the
%   sensor hardware. For each incoming event, its timestamp is written
%   into the timestamp memory of all eight immediately neighbouring
%   pixels. The event is then accepted only if the timestamp already
%   stored at its own pixel location indicates that a nearby event
%   occurred within the support time T. Because genuine activity from
%   moving objects produces spatiotemporally correlated events in
%   adjacent pixels, these events mutually support one another and
%   pass through. Isolated noise events, which lack recent
%   neighbouring activity, are rejected.
%
%   T (support time): The maximum allowable time, in seconds, between
%   the current event and the most recent neighbouring event for the
%   current event to be accepted as correlated activity. Increasing T
%   extends the support window, which is appropriate for slow-moving
%   objects or low event rates where the time between correlated
%   events at adjacent pixels may be long. Decreasing it tightens the
%   window, improving rejection of background activity noise at the
%   risk of discarding valid events from slow or sparse motion. The
%   optimal value depends on the expected speed of objects in the
%   scene relative to the sensor resolution: faster motion produces
%   shorter inter-event intervals at neighbouring pixels, permitting
%   a smaller T.

baf_params.T                = 10e-1;   
baf_lastTimesMap            = -Inf(imgSz);
baf_n_passed_store          = zeros(1, frame_total);
baf_n_total_store           = zeros(1, frame_total);

% SPATIOTEMPORAL CORRELATION FILTER PARAMETERS 
% -------------------------------------------------------------------------
% REF: http://ieeexplore.ieee.org/document/7168735/
% -------------------------------------------------------------------------
% Tuning Guide:
%   The Spatiotemporal Correlation Filter removes background activity
%   (BA) noise by exploiting the fact that real events arising from
%   moving features are temporally correlated with nearby events,
%   whereas BA noise events occur randomly and in isolation. The
%   sensor pixel array is spatially subsampled so that blocks of
%   pixels map onto a single filter cell. When an event arrives at a
%   cell, it opens a time window during which any subsequent event
%   mapping to the same cell is accepted as correlated activity. If
%   no prior event has occurred in the cell within the time window,
%   the incoming event is rejected as uncorrelated noise.
%
%   dT (temporal correlation window): The duration, in seconds, for
%   which a filter cell remains "open" after receiving an event.
%   Any subsequent event that maps to the same cell within this
%   window is classified as correlated and allowed to pass.
%   Increasing dT extends the window, which is appropriate for
%   slow-moving objects or low event rates where the inter-event
%   interval within a cell may be long. Decreasing it tightens the
%   window, improving noise rejection at the risk of discarding
%   valid events from slow or sparse activity. The original hardware
%   design supports a range of 5 ns to 1 s, with typical values of
%   1-10 ms.
%
%   subsample_rate (spatial neighbourhood exponent): Controls the
%   size of the pixel block that maps onto a single filter cell.
%   A value of S causes each block of 2^S x 2^S sensor pixels to
%   share one filter cell (e.g. S=0 gives 1x1, S=1 gives 2x2, S=2
%   gives 4x4). Increasing this value widens the spatial
%   neighbourhood, so that events from a larger region of the sensor
%   provide temporal support for one another. This improves
%   robustness for sparse or slow-moving scenes but reduces spatial
%   selectivity, potentially allowing spatially distant noise events
%   to support each other. Decreasing it restricts support to fewer
%   pixels, preserving spatial precision at the cost of requiring
%   denser local activity for events to pass.

stc_params.dT               = 10e-3;    
stc_params.subsample_rate   = 2;       
stc_block_size              = 2^stc_params.subsample_rate;
stc_cellSz                  = ceil(imgSz ./ stc_block_size);
stc_lastTimesMap            = -Inf(stc_cellSz);
stc_n_passed_store          = zeros(1, frame_total);
stc_n_total_store           = zeros(1, frame_total);

% MOTION CONSISTENCY FILTER
% -------------------------------------------------------------------------
% REF: https://ieeexplore.ieee.org/document/8953966/
% -------------------------------------------------------------------------
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

mcf_params.delta_t         = 0.1;   
mcf_params.V_max           = 500;    
mcf_params.min_neighbours  = 2;       
mcf_eventBuffer.x          = [];
mcf_eventBuffer.y          = [];
mcf_eventBuffer.t          = [];
mcf_n_passed_store         = zeros(frame_total, 1);
mcf_n_total_store          = zeros(frame_total, 1);

% EVENT DENSITY FILTER PARAMETERS
% -------------------------------------------------------------------------
% REF: https://doi.org/10.3390/app10062024
% -------------------------------------------------------------------------
% Tuning Guide:
%   The Event Density Filter classifies each incoming event as signal or
%   noise by counting the number of events that have occurred within a
%   spatiotemporal neighbourhood. If the total count (the event density)
%   falls below a threshold, the event is rejected as background
%   activity. A second pass removes hot-pixel noise by checking whether
%   the surviving events in a small 3x3 neighbourhood have any spatial
%   support from surrounding pixels.
%
%   L (spatial neighbourhood size): Defines the side length, in pixels,
%   of the square region centred on the candidate event within which
%   supporting events are counted. Increasing L draws support from a
%   wider area, improving robustness in sparse scenes but potentially
%   allowing spatially distant noise events to support one another.
%   Decreasing L restricts support to the immediate vicinity, suiting
%   dense scenes with fine spatial structure.
%
%   dt (temporal neighbourhood size): Defines how far back in time, in
%   seconds, the filter looks for supporting events. Increasing dt
%   extends the temporal window, which is appropriate for slow-moving
%   objects or low event rates where supporting events are spread over
%   longer intervals. Decreasing it tightens the window, suiting fast
%   motion where only recent events are likely to be correlated with
%   the candidate.
%
%   Psi (event density threshold): The minimum number of events that
%   must be present in the spatiotemporal neighbourhood for the
%   candidate event to be retained. Raising Psi demands more supporting
%   evidence before accepting an event, producing more aggressive noise
%   removal at the risk of discarding valid signal events at the edges
%   of moving objects. Lowering it retains more events at the risk of
%   admitting background activity noise.

edf_params.L                = 5;        
edf_params.dt               = 10e-1;    
edf_params.Psi              = 3;      
edf_eventBuffer.x           = [];
edf_eventBuffer.y           = [];
edf_eventBuffer.t           = [];
edf_n_passed_store          = zeros(1, frame_total);
edf_n_total_store           = zeros(1, frame_total);

% STCC-FILTER PARAMETERS
% -------------------------------------------------------------------------
% REF: https://doi.org/10.1016/j.image.2024.117136
% -------------------------------------------------------------------------
% Tuning Guide:
%   The STCC-Filter computes a probability-of-signal (POS) score for
%   each incoming event by evaluating its spatiotemporal correlation
%   with N neighbouring events. Events whose POS falls below the
%   discrimination threshold TH are classified as noise and rejected.
%
%   N (context window size): The number of recent events used to
%   evaluate each candidate. Larger values provide more statistical
%   support and improve robustness, but increase computation time.
%   Smaller values are faster but may produce unreliable POS estimates.
%
%   sigma_d (spatial scale): Controls how rapidly the influence of
%   neighbouring events decays with spatial distance, in pixels.
%   Increasing this value extends the spatial support region, which is
%   appropriate for sparse scenes or broad edge structures. Decreasing
%   it restricts support to the immediate vicinity, suiting dense or
%   fine-featured scenes.
%
%   sigma_t (temporal scale): Controls how rapidly the influence of
%   neighbouring events decays with temporal distance, in seconds.
%   Increasing this value extends the temporal support window, which is
%   appropriate for slow-moving objects or low event rates. Decreasing
%   it tightens the window, suiting fast motion where temporally
%   distant events are unlikely to be correlated.
%
%   sigma_p (polarity scale): Controls how strongly polarity agreement
%   is weighted. Decreasing this value enforces stricter same-polarity
%   agreement between neighbouring events. Increasing it relaxes this
%   constraint, permitting greater mixed-polarity support.
%
%   TH (base discrimination threshold): The minimum POS value an event
%   must exceed to be retained. Raising TH produces more aggressive
%   filtering at the risk of discarding valid signal events. Lowering
%   it retains more events at the risk of admitting noise. Note that
%   the content-correlation mechanism adjusts TH at run-time, so this
%   value serves as the initial operating point for that adaptation.

stcc_params.N               = 100;      
stcc_params.sigma_d         = 6;       
stcc_params.sigma_t         = 10e-1;   
stcc_params.sigma_p         = 5;        
stcc_params.TH              = 0.1;     
stcc_n_passed_store         = zeros(1, frame_total);
stcc_n_total_store          = zeros(1, frame_total);

% COHERENCE PARAMETERS (Crain 2026)
% -------------------------------------------------------------------------
% REF: TBD
% -------------------------------------------------------------------------
% Tuning Guide:

coh_params.iei_alpha        = 0.8; 
coh_params.r_s              = 50/imgSz(1);
filter_mask                 = ones(imgSz);
iei_map                     = zeros(imgSz);

%% Initialize Accumulator Parameters

% TIME SURFACE (HOTS) PARAMETERS
% -------------------------------------------------------------------------
% REF: http://ieeexplore.ieee.org/document/7508476/
% -------------------------------------------------------------------------
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

ts_time_constant = 2;  
ts_t_map = -inf(imgSz); 

% SPEED INVARIANT TIME SURFACE (SITS) PARAMETERS
% -------------------------------------------------------------------------
% REF: https://arxiv.org/pdf/1903.11332
% -------------------------------------------------------------------------
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

sits_t_map = zeros(imgSz);
sits_R = 3;

% ADAPTIVE GLOBAL DECAY (AGD) TIME SURFACE PARAMETERS
% -------------------------------------------------------------------------
% REF: https://ieeexplore.ieee.org/document/10205486/
% -------------------------------------------------------------------------
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

agd_params.alpha            = 0.001;   
agd_params.K                = 50000.0;  
agd_surface                 = zeros(imgSz); 
agd_state.last_t_map        = zeros(imgSz);
agd_state.activity          = 0; 
agd_state.last_update_time  = t_start_process; 
agd_activity_store          = zeros(length(frame_total),1);

% MOTION-ENCODED TIME-SURFACE (METS) PARAMETERS
% -------------------------------------------------------------------------
% REF: https://link.springer.com/10.1007/s11263-025-02379-6
% -------------------------------------------------------------------------
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

mets_params.R    = 4;   
mets_params.n    = 3;    
mets_params.d    = 5;   
mets_params.d_th = 8;   
mets_state.t_last_pos = zeros(imgSz);   
mets_state.t_last_neg = zeros(imgSz);  
mets_state.t_last_any = zeros(imgSz);   
mets_state.p_last     = zeros(imgSz);   

% ZHU ADAPTIVE TIME SURFACE (EVO-ATS) PARAMETERS
% -------------------------------------------------------------------------
% REF: https://ieeexplore.ieee.org/document/10342048/
% -------------------------------------------------------------------------
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

zhu_state.t_last = zeros(imgSz);
zhu_params.tau_u       = 0.8;    
zhu_params.tau_l       = 0.001;  
zhu_params.n_neighbors = 16;     
zhu_params.blur_sigma  = 1.0;   
zhu_params.median_sz   = 3;      

% ADAPTIVE LOCAL TIME-SURFACE PARAMETERS (Crain, 2026)
% -------------------------------------------------------------------------
% REF: TBD
% -------------------------------------------------------------------------

alts_params.dt                   = t_interval;
alts_params.filter_size          = 11;
alts_params.filter_sigma         = 9.0;
alts_params.surface_tau_release  = 3.0;
alts_params.div_norm_exp         = 1.0;
alts_params.symmetric_tone_scale = 1.0;
alts_activity_score.mean         = zeros(frame_total, 1);
alts_activity_score.median       = zeros(frame_total, 1);
alts_activity_score.std          = zeros(frame_total, 1);

%% Initialize all figure code for video output
% Indicate which videos should be saved
cohOut = false;
atsOut = true;

% Initialize the videos
[hFigs, hAxs, hImgs, videoWriters] = plot.initializeEventVideos(cohOut,...
    atsOut, imgSz, videoOutPath);

%% Initialize data storage and perform data optimizations
% Identify number of events
current_idx     = 1;

if useBuffer == false
    n_events        = length(tk);
    frame_time      = zeros(frame_total, 1);
end

% Initialize per-pixel timestamp tracking
last_event_timestamp    = zeros(imgSz);
norm_trace_map_prev     = zeros(imgSz);
time_surface_map_prev   = zeros(imgSz);

% Pre-allocate metrics storage.
frame_metrics.SRR                  = zeros(1, frame_total);
frame_metrics.ClarkEvansRemoved    = zeros(1, frame_total);
frame_metrics.ClarkEvansRemaining  = zeros(1, frame_total);
frame_metrics.ComputeTimeFilter    = zeros(1, frame_total);
frame_metrics.ComputeTimeAccumulator = zeros(1, frame_total);
frame_metrics.EventsInFrame        = zeros(1, frame_total);
frame_metrics.FilteredEvents       = zeros(1, frame_total);
frame_metrics.FilteringMEVs        = zeros(1, frame_total);
frame_metrics.AccumulatorMEVs      = zeros(1, frame_total);

elbowDiagnostics = {};

% Track the previous frame's output for temporal SSIM computation.
% This is separate from time_surface_map_prev (which is the raw
% unnormalized EMA state for feedback). This stores the normalized
% display-ready output for metrics comparison.
prev_output_for_metrics = zeros(imgSz);

%% Data processing starting point
% Loop through the figures to capture each frame
for frameIndex = 1:frame_total 

    if frameIndex == 1
        fprintf('\n')
        fprintf('------------------ OPTIONS -----------------\n')
        fprintf('--------------------------------------------\n')
        fprintf('\n')
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    % ============================= SETUP ================================%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    
    % Increment the interval
    t_range_c = (frameIndex - 1) * t_interval;
    t_range_n = (t_range_c+t_interval);

    % Slice the events to a valid range
    if useBuffer == false
        [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
        process.sliceToValidRange(t_range_n, xk, yk, tk, pk, imgSz, current_idx);
    else
        [x_valid, y_valid, t_valid, p_valid] = buf.nextWindow(t_range_n, imgSz);
    end

    % Confirm the presence of valid events in the packet
    % If no events are present, we skip this frame
    if isempty(t_valid)
        
        fprintf('There are no events in this slice, skipping... \n');
        continue;
        
    end   

    % ------------------------ EVENT PREPERATION--------------------------%
    % --------------------------------------------------------------------%
    
    % Convert 2D subscripts (x,y) to 1D linear indices
    % Imagine you are a post-man with a disorganized stack of letters. 
    % Instead of dealing with letters for
    % 3rd Avenue, 5th street, the sub2ind function assigns and "ID" for
    % each house. So going forward, (3,5) might just be "House #1".
    % Programatically this just means that (1,1) is "1", (1,2) is "2".
    % So you will have a list at the end which is of size x*y. In this
    % case with a frame of size 640 by 480, the MAXIMUM size of the list
    % is 307,200. However the practical size of the list will change as
    % not every pixel is active.
    linear_idx = sub2ind(imgSz, x_valid, y_valid);
    
    % So now that we have a linear index, we sort them by "ID". So
    % ultimately what we get here is a sorted list of event times where we
    % have all event time differences for each pixel coordinate groupped. 
    % So the time intervals may not actually be monotonic. So
    % sorted_idx is basically the list of "houses" grouped together.
    [sorted_idx, sort_order] = sort(linear_idx);
    sorted_t = t_valid(sort_order);
    sorted_x = x_valid(sort_order);
    sorted_y = y_valid(sort_order);
    sorted_p = p_valid(sort_order);

    % Reset the frames for the current loop
    counts = zeros(imgSz);
    
    % With the sorted index list, aka the houses, we can extract where
    % in that list each "house" starts. That would be the "pos" output.
    % The unique_idx list is the complete list of all houses. 
    [unique_idx, pos, ~] = unique(sorted_idx);
           
    % Now we need to find where each "house" starts and ends. The start
    % is easy because we get it from "pos". The end can be inferred
    % from the start because the start of the NEXT "house" group is
    % always one more than the end of the previous group. So,
    % pos(2:end)-1 gives the list of ends. We do not know when the last
    % group will end, but we do know that it MUST end by the end of the
    % dataset. So this is the length(sorted_idx).
    group_ends = [pos(2:end)-1; length(sorted_idx)];

    % This is a bit of cheeky MATLAB code. Because "unique_idx" is a
    % linear index list, MATLAB knows to map it to the size of the 2D
    % array "counts". So although group_end and pos are vectors, the
    % output is a 640x480 array. Additionally, the ending position of
    % each "house" minus the starting position (i.e. the difference)
    % gives you the total number of "things" assigned to that "ID". 
    counts(unique_idx) = group_ends - pos + 1;

    % ------------------------- STATISTICS -------------------------------%
    % --------------------------------------------------------------------%
    [t_mean, t_std, t_max, t_min, t_mean_diff, t_std_diff] = ...
        stats.computeNeighborhoodStats(sorted_t, unique_idx, pos, ...
        group_ends, imgSz);

    % Update a persistant map of the inter-event interval so that sparse
    % data is retained
    new_obs_mask = (t_mean_diff > 0); 
    iei_map(new_obs_mask) = (1 - coh_params.iei_alpha) .* iei_map(new_obs_mask) +...
        coh_params.iei_alpha .* t_mean_diff(new_obs_mask);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    % =========================== FILTERING ==============================%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    
    % Start filter timer
    filtering_start = tic;

    if frameIndex == 1
        fprintf(['FILTER SELECTION: ' filterSelection '\n']);
    end

    if strcmp(filterSelection, 'STC') == 1

        % -------------------- STC CORRELATION FILTER --------------------%
        % ----------------------------------------------------------------%

        [filter_mask, stc_lastTimesMap, stc_n_pass, stc_n_tot] = ...
            filters.spatiotemporalCorrelation(sorted_x, sorted_y, ...
            sorted_t, imgSz, stc_lastTimesMap, stc_params);

        % Store metrics for post-processing analysis
        stc_n_passed_store(frameIndex) = stc_n_pass;
        stc_n_total_store(frameIndex)  = stc_n_tot;

    elseif strcmp(filterSelection, 'BAF') == 1

        % ------------- BACKGROUND ACTIVITY FILTER (BAF) -----------------%
        % ----------------------------------------------------------------%

        [filter_mask, baf_lastTimesMap, baf_n_pass, baf_n_tot] = ...
            filters.backgroundActivityFilter(sorted_x, sorted_y, ...
            sorted_t, imgSz, baf_lastTimesMap, baf_params);

        % Store metrics for post-processing analysis
        baf_n_passed_store(frameIndex) = baf_n_pass;
        baf_n_total_store(frameIndex)  = baf_n_tot;

    elseif strcmp(filterSelection, 'EDF') == 1

        % ------------ EVENT DENSITY FILTER (EDF) ------------------------%
        % ----------------------------------------------------------------%

        [filter_mask, edf_eventBuffer, edf_n_pass, edf_n_tot] = ...
            filters.eventDensityFilter(sorted_x, sorted_y, ...
            sorted_t, imgSz, edf_eventBuffer, edf_params);

        % Store metrics for post-processing analysis
        edf_n_passed_store(frameIndex) = edf_n_pass;
        edf_n_total_store(frameIndex)  = edf_n_tot;

    elseif strcmp(filterSelection, 'STCC') == 1

        % ------- STCC-FILTER (Space-Time-Content Correlation) -----------%
        % ----------------------------------------------------------------%

        % Compute signed polarity for this filter
        p_signed_for_stcc = sorted_p * 2 - 1;

        [filter_mask, stcc_n_pass, stcc_n_tot, stcc_threshold] = ...
            filters.stccFilter(sorted_x, sorted_y, ...
            sorted_t, p_signed_for_stcc, imgSz, stcc_params);

        % Store metrics for post-processing analysis
        stcc_n_passed_store(frameIndex) = stcc_n_pass;
        stcc_n_total_store(frameIndex)  = stcc_n_tot;
        frame_metrics.STCCThreshold(frameIndex) = stcc_threshold;

    elseif strcmp(filterSelection, 'MCF') == 1
 
        % ----------- MOTION CONSISTENCY FILTER (MCF) ----------------%
        % ------------------------------------------------------------%

        [filter_mask, mcf_eventBuffer, mcf_n_pass, mcf_n_tot] = ...
            filters.motionConsistencyFilter(sorted_x, sorted_y, ...
            sorted_t, imgSz, mcf_eventBuffer, mcf_params);
 
        % Store metrics for post-processing analysis
        mcf_n_passed_store(frameIndex) = mcf_n_pass;
        mcf_n_total_store(frameIndex)  = mcf_n_tot;
        

    elseif strcmp(filterSelection, 'COH') == 1

        % --------------------- COHERENCE FILTER -------------------------%
        % ----------------------------------------------------------------%
        [norm_trace_map, norm_similarity_map, norm_persist_map,...
            filtered_coherence_map] = filters.computeCoherenceMask(sorted_x,...
            sorted_y, sorted_t, imgSz, t_interval, unique_idx, pos, ...
            group_ends, coh_params, frameIndex, norm_trace_map_prev, t_std_diff,...
            t_mean_diff, counts);

        % Set any retention variables
        norm_trace_map_prev = norm_trace_map;
        
        % Create the filter
        filter_mask = filtered_coherence_map;
        filter_mask(isnan(filter_mask)) = 0;
        filter_mask = single(imgaussfilt(single(filter_mask), 5.0, "FilterSize", 9));

        % Use the Kneedle algorithm to find a suitable threshold
        % Downsample filter_mask from 640x480 to 320x240 before thresholding
        % Assume imgSz corresponds to original size; compute target size half in each dim
        targetSz = ceil(imgSz / 2);
        % Use imresize to downsample (preserve range); convert to single for processing
        ds_mask = imresize(single(filter_mask), targetSz, 'bilinear');
        [th_lo, diag] = stats.findElbowThreshold(ds_mask);
        % [th_lo, diag] = stats.findElbowThreshold(norm_trace_map);
        % norm_trace_map(norm_trace_map<th_lo)=0;
        % [th_lo, diag] = stats.findElbowThreshold(norm_persist_map);
        % norm_persist_map(norm_persist_map<th_lo)=0;
        % [th_lo, diag] = stats.findElbowThreshold(norm_similarity_map);
        % norm_similarity_map(norm_similarity_map<th_lo)=0;

        % filter_mask = norm_similarity_map + norm_persist_map;
          
        % Use the threshold on the mask & log the result
        filter_mask(filter_mask < th_lo) = 0;
        % last_pass = bwareaopen(filter_mask,100);
        % filter_mask = filter_mask.*last_pass.*1.0;

        % Compute event-level pass/total counts for the COH filter.
        % The other filters return these directly; COH needs them
        % derived from the pixel-level mask and event locations.
        coh_linear_idx = sub2ind(imgSz, sorted_x(:), sorted_y(:));
        coh_n_total  = numel(sorted_x);
        coh_n_passed = sum(filter_mask(coh_linear_idx) > 0);

        frame_metrics.FilterThreshold(frameIndex) = th_lo;
        elbowDiagnostics{frameIndex} = diag;

    elseif strcmp(filterSelection, 'NONE') == 1

        % Use a unity mask instead
        filter_mask = (counts>0).*1.0;

    end

    % Stop filter timer
    filtering_stop = toc(filtering_start);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    % ========================= ACCUMULATION =============================%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    
    % Start accumulator timer
    accumulator_start = tic;

    if frameIndex == 1
        fprintf(['ACCUMULATOR SELECTION: ' accumulatorSelection '\n']);
    end

    % Any events which fall within the filter mask should not be included
    % in the accumulation
    event_linear_idx = sub2ind(imgSz, sorted_x, sorted_y);
    keep = filter_mask(event_linear_idx) > 0;
    
    % Apply to all event vectors simultaneously
    filtered_x = sorted_x(keep);
    filtered_y = sorted_y(keep);
    filtered_t = sorted_t(keep);
    filtered_p = sorted_p(keep);

    % Ensure polarity is -1 and 1 (if it's 0 and 1)
    p_signed = double(filtered_p);
    p_signed(p_signed == 0) = -1;

    if strcmp(accumulatorSelection, 'IEI-ATS') == 1

        % ------------ ADAPTIVE LOCAL TIME-SURFACE UPDATE ----------------%
        % ----------------------------------------------------------------%

        % Accumulate polarity into a 2D grid
        % If multiple events land on one pixel, we sum their polarities 
        polarity_map = accumarray([filtered_x, filtered_y], p_signed, imgSz, @sum, 0);
    
        [normalized_output_frame, time_surface_map_raw, tau_filtered, adaptive_gains] = ...
        accumulator.localAdaptiveTimeSurface(iei_map,...
        time_surface_map_prev, alts_params, filter_mask, polarity_map, counts);
    
        % Set any retention variables
        time_surface_map_prev = time_surface_map_raw;

        % Store the adaptive map score
        alts_activity_score.mean(frameIndex) = mean(adaptive_gains...
            (abs(adaptive_gains)>0));
        alts_activity_score.median(frameIndex) = median(adaptive_gains...
            (abs(adaptive_gains)>0));
        alts_activity_score.std(frameIndex) = std(adaptive_gains...
            (abs(adaptive_gains)>0));
    
    elseif strcmp(accumulatorSelection, 'AGD') == 1

        % ------------------- NUNES GLOBAL ADAPTIVE (AGD) ----------------%
        % ----------------------------------------------------------------%
    
        % Run the AGD algorithm
        [agd_surface, agd_state, ~] = accumulator.adaptiveGlobalDecay(agd_surface,...
            filtered_x, filtered_y, filtered_t, agd_state, agd_params);

        normalized_output_frame = agd_surface;
        agd_activity_store(frameIndex) = agd_state.activity;

    elseif strcmp(accumulatorSelection, 'HOTS') == 1        
    
        % ------------- TIME-SURFACE ACCUMULATION (HOTS) -----------------%
        % ----------------------------------------------------------------%
        
        [ts_t_map, normalized_output_frame] = accumulator.timeSurface(ts_t_map,...
         filtered_x, filtered_y, filtered_t, imgSz, ts_time_constant);

    elseif strcmp(accumulatorSelection, 'SITS') == 1

        % ------------ SPEED INVARIENT TIME-SURFACE (SITS) ---------------%
        % ----------------------------------------------------------------%
        
        [sits_t_map, normalized_output_frame] = ...
            accumulator.speedInvariantTimeSurface(sits_t_map, filtered_x,...
            filtered_y, sits_R);

    elseif strcmp(accumulatorSelection, 'METS') == 1

        % ------------ MOTION-ENCODED TIME-SURFACE (METS) ----------------%
        % ----------------------------------------------------------------%

        [mets_surface, mets_state, normalized_output_frame] = ...
            accumulator.motionEncodedTimeSurface(...
            filtered_x, filtered_y, filtered_t, p_signed, ...
            imgSz, mets_state, mets_params);

    elseif strcmp(accumulatorSelection, 'EVO-ATS') == 1

        % ---------- ZHU ADAPTIVE TIME SURFACE (EVO-ATS) -----------------%
        % ----------------------------------------------------------------%
    
        [zhu_surface, zhu_state, normalized_output_frame, zhu_tau_map] = ...
            accumulator.adaptiveTimeSurfaceZhu(...
            filtered_x, filtered_y, filtered_t, p_signed, ...
            imgSz, zhu_state, zhu_params);

    else

        error('An invalid accumulator has been selected!')

    end

    % Start accumulator timer
    accumulator_stop = toc(accumulator_start);

    if frameIndex == 1
        fprintf('\n--------------------------------------------\n\n');
    end

    % ========================== FRAME METRICS ===========================%
    % --------------------------------------------------------------------%
    
    % Determine n_passed and n_total from whichever filter ran.
    % Each filter stores these differently; unify them here.
    recovered_events = process.recoverRemovedEvents(filter_mask, counts);
    switch filterSelection
        case 'STC'
            metrics_n_passed = stc_n_pass;
            metrics_n_total  = stc_n_tot;
            metrics_background_map = recovered_events;
            metrics_foreground_map = filter_mask;
        case 'BAF'
            metrics_n_passed = baf_n_pass;
            metrics_n_total  = baf_n_tot;
            metrics_background_map = recovered_events;
            metrics_foreground_map = filter_mask;
        case 'EDF'
            metrics_n_passed = edf_n_pass;
            metrics_n_total  = edf_n_tot;
            metrics_background_map = recovered_events;
            metrics_foreground_map = filter_mask;
        case 'STCC'
            metrics_n_passed = stcc_n_pass;
            metrics_n_total  = stcc_n_tot;
            metrics_background_map = recovered_events;
            metrics_foreground_map = filter_mask;
        case 'COH'
            metrics_n_passed = coh_n_passed;
            metrics_n_total  = coh_n_total;

            % Coherence requires a special calculation due to the Gaussian
            % blur applied
            temporary_filter_mask = recovered_events.*counts;
            metrics_background_map = (temporary_filter_mask>0).*1.0;
            metrics_foreground_map = (filter_mask>0).*1.0;
        case 'MCF'
            metrics_n_passed = mcf_n_pass;
            metrics_n_total  = mcf_n_tot;
            metrics_background_map = recovered_events;
            metrics_foreground_map = filter_mask;
        case 'NONE'
            metrics_n_passed = numel(sorted_x);
            metrics_n_total  = numel(sorted_x);
            metrics_background_map = (counts>0).*1.0;
            metrics_foreground_map = (counts>0).*1.0;
        otherwise
            metrics_n_passed = numel(sorted_x);
            metrics_n_total  = numel(sorted_x);
            metrics_background_map = (counts>0).*1.0;
            metrics_foreground_map = (counts>0).*1.0;
    end
    
    % Compute all per-frame metrics
    frame_metrics.SRR(frameIndex) = ...
        metrics.computeSignalRetentionRate(metrics_n_passed, metrics_n_total);
    frame_metrics.ClarkEvansRemoved(frameIndex) = ...
        metrics.computeClarkEvans(metrics_background_map);
    frame_metrics.ClarkEvansRemaining(frameIndex) = ...
        metrics.computeClarkEvans(metrics_foreground_map);
    frame_metrics.ComputeTimeFilter(frameIndex) = filtering_stop;
    frame_metrics.ComputeTimeAccumulator(frameIndex) = accumulator_stop;
    frame_metrics.EventsInFrame(frameIndex) = metrics_n_total;
    frame_metrics.FilteredEvents(frameIndex) = ...
        metrics_n_total-metrics_n_passed;
    frame_metrics.FilteringMEVs(frameIndex) = ...
        metrics_n_total/filtering_stop;
    frame_metrics.AccumulatorMEVs(frameIndex) = ...
        metrics_n_total/accumulator_stop;

    % Update the previous frame reference for next iteration
    prev_output_for_metrics = normalized_output_frame;

    % ------------------------ EXPORTING VIDEO ---------------------------%
    % --------------------------------------------------------------------%

    % Convert to uint8 (0-255 range)
    grayscale_normalized_output_frame = ...
        uint8(normalized_output_frame .* 255);

    if atsOut
        % Capture the frame for the video writer
        set(hImgs{1}, 'CData', grayscale_normalized_output_frame');
        set(hImgs{1}, 'AlphaData', ~isnan(grayscale_normalized_output_frame'));
        set(hAxs{1}, 'Visible','off');
        colormap(gray);
        clim([0 255]);
        writeVideo(videoWriters{1}, grayscale_normalized_output_frame'); 
        
        frameOutputFolder =  string(videoOutPath(1:end-4)) + "-FRAMES";
        % Also write each frame as a PNG to a folder
        % Ensure output folder exists
        if frameIndex == 1
            if ~exist(frameOutputFolder, 'dir')
                mkdir(frameOutputFolder);
            else 
                rmdir(frameOutputFolder, 's');
                mkdir(frameOutputFolder); % Create the output folder
            end
        end
        % Build filename with zero-padded frame index
        fname = fullfile(frameOutputFolder, sprintf('frame_%05d.png', frameIndex));
        imwrite(grayscale_normalized_output_frame', fname);
    end

    % Print progress
    stats.printPercentComplete(frameIndex, frame_total, ...
        frame_metrics.ComputeTimeFilter(frameIndex), ...
        frame_metrics.ComputeTimeAccumulator(frameIndex));

end

% Close the video writer
for videosIdx = 1:length(videoWriters)
    close(videoWriters{videosIdx});
end
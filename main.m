clear;
clc;
close all;

%% Processing pipeline reminder
% To import AEDAT4 -->
% /home/alexandercrain/Repositories/CNN/import/importAEDAT4toHDF5.py

% All datasets --> /home/alexandercrain/Dropbox/Graduate Documents/
% Doctor of Philosophy/Thesis Research/Datasets/SPOT

% Datasets for MATLAB --> /home/alexandercrain/Dropbox/Graduate Documents/
% Doctor of Philosophy/Thesis Research/Datasets/SPOT/HDF5

%% Select algorithms for filtering and accumulation
% Set 'None' for filter selection to skip filtering entirely

% Set filtering
% Options: 'NONE', 'STC', 'BAF', 'EDF', 'STCC', 'COHERENCE'
filterSelection = 'COHERENCE';

% Set accumulator
% Options: 'HOTS', 'SITS', 'METS', 'IEI-ATS', 'AGD', 'EVO-ATS'
accumulatorSelection = 'IEI-ATS';

%% Define processing range
% Define start and end time to process [seconds]
t_start_process = 0; 
t_end_process   = 1000; 

%% Import events for inspection

% Set path to datasets
hdf5Path = ['/home/alexandercrain/Dropbox/Graduate Documents' ...
    '/Doctor of Philosophy/Thesis Research/Datasets/SPOT/HDF5/'];

% Set dataset name
%fileName = 'recording_20260127_145247.hdf5';  % Jack W. (LED Cont)
fileName = 'recording_20251029_131131.hdf5';  % EVOS - NOM - ROT
%fileName = 'recording_20251029_135047.hdf5';  % EVOS - SG - ROT
%fileName = 'recording_20251029_134602.hdf5';  % EVOS - DARK - ROT

% Set output video name
videoOutPath = fullfile('/home/alexandercrain/Videos/Research', ...
    sprintf('Normalized-Output-Fltr-%s-Acmtr-%s-%s.avi', ...
    filterSelection, accumulatorSelection, datestr(now,'yyyymmdd-HHMMSS'))); %#ok<TNOW1,DATST>

% Load the data
tk = double(h5read([hdf5Path fileName], '/timestamp'));
xk = single(h5read([hdf5Path fileName], '/x'));
yk = single(h5read([hdf5Path fileName], '/y'));
pk = single(h5read([hdf5Path fileName], '/polarity'));

% Convert time to seconds
tk = (tk - tk(1))/1e6;

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

%% Initialize all tunable parameters for the algorithms
% GENERAL PARAMETERS
% ------------------
% Set the image size
imgSz                       = [640, 480]; 

% Set the time interval to accumulate over
t_interval                  = 0.33;  % [s]
t_total                     = max(tk);  % [s]
frame_total                 = floor(t_total/t_interval);

% Boolean controls
filter_output_image         = false;

% INTER-EVENT-INTERVAL
% --------------------
% Persistent inter-event-interval map (EMA) parameter
iei_alpha                   = 0.8;     

% Initialize EMA-IEI 
iei_map                     = zeros(imgSz);

% COHERENCE PARAMETERS
% --------------------
% Define coherence parameters - Must be tuned based on t_interval
coh_params.r_s                         = 30/imgSz(1);  % spatial radius [pixels norm]
% coh_params.trace_threshold             = 1.3;
% coh_params.persistence_threshold       = 0.0002;
% coh_params.coherence_threshold         = 0.06;
% coh_params.similarity_threshold        = 0.5;
coh_params.trace_threshold             = 0.0;
coh_params.persistence_threshold       = inf;
coh_params.coherence_threshold         = 0.0;
coh_params.similarity_threshold        = inf;


coh_logs.trace_threshold               = zeros(frame_total, 1);
coh_logs.persistence_threshold         = zeros(frame_total, 1);
coh_logs.similarity_threshold          = zeros(frame_total, 1);

% Initialize mask for filter
filter_mask                 = ones(imgSz);

% BACKGROUND ACTIVITY FILTER PARAMETERS (Delbruck 2008)
% ---------------------------------------------------------------
% T: Support time window [s]. An event passes if any of its 8
%    neighbours fired within the last T seconds. This is the sole
%    tunable parameter.
%
%    - Shorter T → aggressive filtering, only tight spatiotemporal
%      clusters survive (fast edges with high event density)
%    - Longer T  → permissive, retains slower/sparser activity

baf_params.T                = 10e-1;    % [s] support time window

% Initialize the persistent last-event timestamp map at full
% pixel resolution. -Inf ensures no false passes on startup.
baf_lastTimesMap            = -Inf(imgSz);

% Optional: pre-allocate storage for BAF metrics
baf_n_passed_store          = zeros(1, frame_total);
baf_n_total_store           = zeros(1, frame_total);

% SPATIOTEMPORAL CORRELATION FILTER PARAMETERS (Liu et al. 2015)
% ---------------------------------------------------------------
% dT: Correlation time window [s]. An event passes if the time
%     since the last event at its subsampled cell is less than dT.
%     Equivalent to C*(Vrs-Vth)/I1 in the hardware (Eq. in paper).
%     - Shorter dT → more aggressive filtering, only fast activity
%       passes (good for high-speed edges, rejects slow drift)
%     - Longer dT  → more permissive, retains slower activity
stc_params.dT               = 10e-1;    % [s] correlation window
stc_params.subsample_rate   = 1;        % 2×2 spatial subsampling

% Initialize the persistent last-event timestamp map at cell
% resolution. -Inf ensures the first event at every cell always
% fails the ISI check (matching the hardware's Fig. 5 behaviour
% where e1 is always rejected).
stc_block_size = 2^stc_params.subsample_rate;
stc_cellSz     = ceil(imgSz ./ stc_block_size);
stc_lastTimesMap = -Inf(stc_cellSz);

% Optional: pre-allocate storage for STC metrics (event counts
% per frame, for comparison against coherence filtering)
stc_n_passed_store = zeros(1, frame_total);
stc_n_total_store  = zeros(1, frame_total);

% EVENT DENSITY FILTER PARAMETERS (Feng et al. 2020)
% ---------------------------------------------------------------
% L:   Spatial neighbourhood side length (odd integer). The density
%      matrix D is L x L centred on the event pixel. Larger L gives
%      more spatial context but may blur the decision boundary at
%      object edges. The paper uses L = 5.
%
% dt:  Temporal neighbourhood [s]. Only events within [t - dt, t]
%      contribute to the density count. Must be tuned to match the
%      sensor's BA rate and the target's event generation rate.
%      The paper uses dt = 5 ms for a CeleX-IV at 768x640.
%      For the DVXplorer Micro (EVOS dataset), you may need a
%      larger dt (10-50 ms) due to lower event rates from the
%      slowly rotating spacecraft target.
%
% Psi: Density threshold. An event passes coarse filtering only if
%      the total event count in its L x L x dt neighbourhood >= Psi.
%      The paper uses Psi = 3. Lower values are more permissive;
%      higher values reject more aggressively.
%
% Two-stage operation:
%   Stage 1 (Coarse): Event density d = ||D||_1 >= Psi?
%   Stage 2 (Fine):   Hot pixel check — are there other passed
%                     events in the 3x3 neighbourhood?

edf_params.L                = 5;        % Spatial neighbourhood size
edf_params.dt               = 10e-1;    % [s] temporal window
edf_params.Psi              = 3;        % Density threshold

% Initialize the cross-frame event buffer (empty)
edf_eventBuffer.x           = [];
edf_eventBuffer.y           = [];
edf_eventBuffer.t           = [];

% Optional: pre-allocate storage for EDF metrics
edf_n_passed_store          = zeros(1, frame_total);
edf_n_total_store           = zeros(1, frame_total);

% STCC-FILTER PARAMETERS (Li et al. 2024)
% ---------------------------------------------------------------
% N:       Number of neighbouring events in the context window for
%          POS and W(p) computation. The paper uses N = 1000.
%          For the EVOS dataset with ~4000 events/frame at
%          t_interval = 0.33 s, N = 500 is a reasonable starting
%          point to balance accuracy and speed.
%
% sigma_d: Std dev of the spatial distance Gaussian [pixels].
%          Controls how sharply spatial proximity is weighted.
%          Paper default: 8 pixels.
%
% sigma_t: Std dev of the temporal distance Gaussian [seconds].
%          Controls how sharply temporal proximity is weighted.
%          Paper default: 5e-3 (5 ms). For the slower EVOS
%          dynamics, you may need 10–30 ms.
%
% sigma_p: Std dev of the polarity distance Gaussian.
%          Paper default: 1. Since |Δp|² ∈ {0, 1} (for 0/1
%          polarity) or {0, 4} (for ±1 polarity), this controls
%          the ratio of weight given to polarity-matched vs
%          mismatched neighbours.
%
% TH:      Base discrimination threshold for POS.
%          Events pass if POS > TH / W(p).
%          This is the primary tuning knob. Lower TH → more
%          permissive (more events pass). Higher TH → more
%          aggressive filtering. Start with 0.05–0.2 and adjust
%          based on visual inspection of the output.
%
% COMPUTATIONAL NOTE:
%   The STCC-Filter is O(M × N) per frame. With M = 4000 events
%   and N = 500, that's 2M Gaussian evaluations per frame.
%   This is slower than BAF/STC (O(M)) but comparable to the
%   coherence filter's KD-tree operations. If speed is critical,
%   reduce N.

stcc_params.N               = 500;      % Context window size
stcc_params.sigma_d         = 8;        % [pixels]
stcc_params.sigma_t         = 10e-3;    % [s] (wider than paper for EVOS)
stcc_params.sigma_p         = 1;        % Polarity Gaussian std dev
stcc_params.TH              = 0.1;      % Base threshold (tune)

% No persistent state needed (stateless filter)

% Optional: pre-allocate storage for STCC metrics
stcc_n_passed_store         = zeros(1, frame_total);
stcc_n_total_store          = zeros(1, frame_total);

% ADAPTIVE LOCAL TIME-SURFACE PARAMETERS
% --------------------------------------
% Set adaptive local time-surface parameters
alts_params.dt                   = t_interval;
alts_params.filter_size          = 11;
alts_params.filter_sigma         = 9.0;
alts_params.surface_tau_release  = 3.0;
alts_params.div_norm_exp         = 1.0;
alts_params.symmetric_tone_scale = 3.0;
alts_activity_score.mean         = zeros(frame_total, 1);
alts_activity_score.median       = zeros(frame_total, 1);
alts_activity_score.std          = zeros(frame_total, 1);

% Create a cell array to store per-frame data (preallocate for frame_total)
alts_frame_storage      = cell(frame_total,1);

% TIME SURFACE PARAMETERS
% -----------------------
% Initialize map
ts_t_map = -inf(imgSz); 

% Decay constant for the visual 
ts_time_constant = 0.05;  % [seconds]

% SPEED INVARIENT TIME SURFACE PARAMETERS
% ---------------------------------------
% REF: https://arxiv.org/pdf/1903.11332
% ---------------------------------------
% Initialize map
sits_t_map = zeros(imgSz);

% Radius of neighbourhood
sits_R = 3;

% ADAPTIVE GLOBAL DECAY TIME SURFACE PARAMETERS
% ---------------------------------------------
% REF: https://ieeexplore.ieee.org/document/10205486/
% ---------------------------------------------
% Initialize map
agd_surface = zeros(imgSz); 

% State structure to hold history
agd_state.last_t_map = zeros(imgSz);  % Stores timestamp of last event
agd_state.activity = 0;  % Initial Activity level
agd_state.last_update_time = t_start_process; 

% Tuning parameters
% Note: for unfiltered events with large spikes, an aggressive smoothing 
% factor is needed. Otherwise, the scene rests anytime it sees a hot pixel
agd_params.alpha   = 0.001;   % Smoothing factor for activity 
agd_params.K       = 5000000.0;  % Scaling factor (Controls "memory length")
agd_activity_store = zeros(length(frame_total),1);

% MOTION-ENCODED TIME-SURFACE (METS) PARAMETERS
% ----------------------------------------------
% REF: Xu et al., "METS: Motion-Encoded Time-Surface for Event-Based
%      High-Speed Pose Tracking," IJCV, vol. 133, pp. 4401-4419, 2025.
%      DOI: 10.1007/s11263-025-02379-6
% ----------------------------------------------
% State: polarity-separated timestamp maps (Eq. 4)
mets_state.t_last_pos = zeros(imgSz);   % Last event timestamp, positive polarity
mets_state.t_last_neg = zeros(imgSz);   % Last event timestamp, negative polarity
mets_state.t_last_any = zeros(imgSz);   % Last event timestamp, either polarity
mets_state.p_last     = zeros(imgSz);   % Polarity of last event at each pixel

% Parameters (Section 3.3 of the paper — identical to their defaults)
mets_params.R    = 4;    % Observation window half-size [pixels]
mets_params.n    = 3;    % Velocity estimation range [pixels]
mets_params.d    = 5;    % Decay step [pixels]
mets_params.d_th = 8;    % Decay distance threshold [pixels]

% ZHU ADAPTIVE TIME SURFACE (ATS) PARAMETERS
% -------------------------------------------
% REF: Zhu, S. et al., "Event Camera-based Visual Odometry for
%      Dynamic Motion Tracking of a Legged Robot Using Adaptive
%      Time Surface," IEEE/RSJ IROS, 2023, pp. 3475-3482.
%      DOI: 10.1109/IROS55552.2023.10342048
% -------------------------------------------
% State: single timestamp map (no polarity separation)
zhu_state.t_last = zeros(imgSz);

% Parameters (paper does not specify exact values; these are
% reasonable defaults for the EVOS dataset frame intervals)
zhu_params.tau_u       = 0.5;    % Upper bound on decay [s]
zhu_params.tau_l       = 0.01;   % Lower bound on decay [s]
zhu_params.n_neighbors = 8;      % # of most-recent neighbors (paper uses n)
zhu_params.blur_sigma  = 1.0;    % Gaussian blur sigma [pixels]
zhu_params.median_sz   = 3;      % Median filter kernel size [pixels, odd]

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
n_events        = length(tk);
frame_time      = zeros(frame_total, 1);

% Initialize per-pixel timestamp tracking
last_event_timestamp    = zeros(imgSz);
norm_trace_map_prev     = zeros(imgSz);
time_surface_map_prev   = zeros(imgSz);

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

    % Start loop timer
    tic;
    
    % Increment the interval
    t_range_c = (frameIndex - 1) * t_interval;
    t_range_n = (t_range_c+t_interval);

    % Slice the events to a valid range
    [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
    process.sliceToValidRange(t_range_n, xk, yk, tk, pk, imgSz, current_idx);

    % Store the valid event data in the output structure
    output_struct.x_valid{frameIndex} = x_valid;
    output_struct.y_valid{frameIndex} = y_valid;
    output_struct.t_valid{frameIndex} = t_valid;
    output_struct.p_valid{frameIndex} = p_valid;

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

    % Ensure polarity is -1 and 1 (if it's 0 and 1)
    p_signed = double(sorted_p);
    p_signed(p_signed == 0) = -1;

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
    iei_map(new_obs_mask) = (1 - iei_alpha) .* iei_map(new_obs_mask) +...
        iei_alpha .* t_mean_diff(new_obs_mask);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    % =========================== FILTERING ==============================%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

    if frameIndex == 1
        fprintf(['FILTER SELECTION: ' filterSelection '\n']);
    end

    if strcmp(filterSelection, 'STC') == 1

        % -------------------- STC CORRELATION FILTER --------------------%
        % ----------------------------------------------------------------%

        % Liu et al. (2015) spatiotemporal correlation filter.
        % Processes events in temporal order within the frame, checking
        % each event's ISI against the programmable dT window. Events
        % with ISI < dT at their subsampled cell are passed; others
        % are rejected as uncorrelated background activity.
        %
        % The lastTimesMap persists across frames, providing temporal
        % continuity: the first event in a new frame can still be
        % correlated with the last event of the previous frame at
        % the same cell.

        [filter_mask, stc_lastTimesMap, stc_n_pass, stc_n_tot] = ...
            filters.spatiotemporalCorrelation(sorted_x, sorted_y, ...
            sorted_t, imgSz, stc_lastTimesMap, stc_params);

        % Store metrics for post-processing analysis
        stc_n_passed_store(frameIndex) = stc_n_pass;
        stc_n_total_store(frameIndex)  = stc_n_tot;

    elseif strcmp(filterSelection, 'BAF') == 1

        % ------------- BACKGROUND ACTIVITY FILTER (BAF) -----------------%
        % ----------------------------------------------------------------%
        % Delbruck (2008) neighbour-support BA filter.
        % For each event: (1) check if any 8-connected neighbour
        % fired within the last T seconds by reading the timestamp
        % at the current pixel (which was written by a neighbour),
        % (2) write the current timestamp to all 8 neighbour
        % locations for future events to check against.

        [filter_mask, baf_lastTimesMap, baf_n_pass, baf_n_tot] = ...
            filters.backgroundActivityFilter(sorted_x, sorted_y, ...
            sorted_t, imgSz, baf_lastTimesMap, baf_params);

        % Store metrics for post-processing analysis
        baf_n_passed_store(frameIndex) = baf_n_pass;
        baf_n_total_store(frameIndex)  = baf_n_tot;

    elseif strcmp(filterSelection, 'EDF') == 1

        % ------------ EVENT DENSITY FILTER (EDF) ------------------------%
        % ----------------------------------------------------------------%
        % Feng et al. (2020) two-stage event density denoising.
        %
        % Stage 1: Count events in the L x L spatial x dt temporal
        % neighbourhood for each event. Reject if count < Psi.
        %
        % Stage 2: From Stage 1 survivors, reject hot pixels by
        % checking for the presence of other survivors in the 3x3
        % neighbourhood (masking out the centre pixel).

        [filter_mask, edf_eventBuffer, edf_n_pass, edf_n_tot] = ...
            filters.eventDensityFilter(sorted_x, sorted_y, ...
            sorted_t, imgSz, edf_eventBuffer, edf_params);

        % Store metrics for post-processing analysis
        edf_n_passed_store(frameIndex) = edf_n_pass;
        edf_n_total_store(frameIndex)  = edf_n_tot;

    elseif strcmp(filterSelection, 'STCC') == 1

        % ------- STCC-FILTER (Space-Time-Content Correlation) -----------%
        % ----------------------------------------------------------------%
        % Li et al. (2024) Gaussian-weighted POS with polarity-based
        % content correlation self-adjusting threshold.
        %
        % For each event:
        %   1. Compute POS = Σ G(Δd, Δt) over N neighbours (Eq. 10)
        %   2. Compute W(p) = Σ G(Δp) over N neighbours (Eq. 16)
        %   3. Self-adjust threshold: TH' = TH / W(p) (Eq. 17)
        %   4. Pass if POS > TH' (Eq. 18)

        % Compute signed polarity for this filter
        p_signed_for_stcc = sorted_p * 2 - 1;

        [filter_mask, stcc_n_pass, stcc_n_tot] = ...
            filters.stccFilter(sorted_x, sorted_y, ...
            sorted_t, p_signed_for_stcc, imgSz, stcc_params);

        % Store metrics for post-processing analysis
        stcc_n_passed_store(frameIndex) = stcc_n_pass;
        stcc_n_total_store(frameIndex)  = stcc_n_tot;
        

    elseif strcmp(filterSelection, 'COHERENCE') == 1

        % --------------------- COHERENCE FILTER -------------------------%
        % ----------------------------------------------------------------%

        [norm_trace_map, norm_similarity_map, norm_persist_map,...
            filtered_coherence_map, trace_threshold_log, ...
            persistence_threshold_log, similarity_threshold_log] = filters.computeCoherenceMask(sorted_x,...
            sorted_y, sorted_t, imgSz, t_interval, unique_idx, pos, ...
            group_ends, coh_params, frameIndex, norm_trace_map_prev, t_mean_diff);
    
        coh_logs.trace_threshold(frameIndex,1) = trace_threshold_log;
        coh_logs.persistence_threshold(frameIndex,1) = persistence_threshold_log;
        coh_logs.similarity_threshold(frameIndex,1)  = similarity_threshold_log;
        
        % Set any retention variables
        norm_trace_map_prev = norm_trace_map;

        % Create the filter
        filter_mask = filtered_coherence_map;
        filter_mask(isnan(filter_mask)) = 0;
        filter_mask = single(imgaussfilt(single(filter_mask), 5.0, "FilterSize", 9));
        filter_mask(filter_mask<coh_params.coherence_threshold) = 0;

    elseif strcmp(filterSelection, 'NONE') == 1

        % Use a unity mask instead
        filter_mask = ones(imgSz);

    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    % ========================= ACCUMULATION =============================%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

    if frameIndex == 1
        fprintf(['ACCUMULATOR SELECTION: ' accumulatorSelection '\n']);
    end

    if strcmp(accumulatorSelection, 'IEI-ATS') == 1

        % ------------ ADAPTIVE LOCAL TIME-SURFACE UPDATE ----------------%
        % ----------------------------------------------------------------%

        % Accumulate polarity into a 2D grid
        % If multiple events land on one pixel, we sum their polarities 
        polarity_map = accumarray([sorted_x, sorted_y], p_signed, imgSz, @sum, 0);
    
        [normalized_output_frame, time_surface_map_raw, tau_filtered, adaptive_gains] = ...
        accumulator.localAdaptiveTimeSurface(iei_map,...
        time_surface_map_prev, alts_params, filter_mask, polarity_map, counts);
    
        % Set any retention variables
        time_surface_map_prev = time_surface_map_raw;

        % Store the processed frames for later use
        alts_frame_storage{frameIndex} = normalized_output_frame;

        % Store the adaptive map score
        alts_activity_score.mean(frameIndex) = mean(adaptive_gains...
            (abs(adaptive_gains)>0));
        alts_activity_score.median(frameIndex) = median(adaptive_gains...
            (abs(adaptive_gains)>0));
        alts_activity_score.std(frameIndex) = std(adaptive_gains...
            (abs(adaptive_gains)>0));

    
    elseif strcmp(accumulatorSelection, 'AGD') == 1

        % --------------- NUNES GLOBAL ADAPTIVE ACCUMULATION--------------%
        % ----------------------------------------------------------------%
    
        % Run the AGD algorithm
        [agd_surface, agd_state, ~] = accumulator.adaptiveGlobalDecay(agd_surface,...
            sorted_x, sorted_y, sorted_t, agd_state, agd_params);
    
        % Store the surface into the standard normalized frame
        normalized_output_frame = agd_surface;
    
        % Store activity data for later inspection
        agd_activity_store(frameIndex) = agd_state.activity;

    elseif strcmp(accumulatorSelection, 'HOTS') == 1        
    
        % ------------------- TIME-SURFACE ACCUMULATION ------------------%
        % ----------------------------------------------------------------%
        
        % Run the normal time-surface accumulation algorithm
        [ts_t_map, normalized_output_frame] = accumulator.timeSurface(ts_t_map,...
         sorted_x, sorted_y, sorted_t, imgSz, ts_time_constant);

    elseif strcmp(accumulatorSelection, 'SITS') == 1

        % -------- SPEED INVARIENT TIME-SURFACE ACCUMULATION -------------%
        % ----------------------------------------------------------------%
        
        % % Run the speed invarient time-surface accumulation algorithm
        % [sits_t_map, normalized_output_frame] = ...
        %     accumulator.speedInvariantTimeSurface(sits_t_map, sorted_x,...
        %     sorted_y, sits_R);

    elseif strcmp(accumulatorSelection, 'METS') == 1

        % ------------ MOTION-ENCODED TIME-SURFACE (METS) ----------------%
        % ----------------------------------------------------------------%
    
        % % Run Motion Encoded Time Surface algorithm
        % [mets_surface, mets_state, normalized_output_frame] = ...
        %     accumulator.motionEncodedTimeSurface(...
        %     sorted_x, sorted_y, sorted_t, p_signed, ...
        %     imgSz, mets_state, mets_params);

    elseif strcmp(accumulatorSelection, 'EVO-ATS') == 1

        % ------------ ZHU ADAPTIVE TIME SURFACE (ATS) -------------------%
        % ----------------------------------------------------------------%
    
        [zhu_surface, zhu_state, normalized_output_frame, zhu_tau_map] = ...
            accumulator.adaptiveTimeSurfaceZhu(...
            sorted_x, sorted_y, sorted_t, p_signed, ...
            imgSz, zhu_state, zhu_params);

    else

        error('An invalid accumulator has been selected!')

    end

    if frameIndex == 1
        fprintf('\n--------------------------------------------\n\n');
    end

    % ------------------------ EXPORTING VIDEO ---------------------------%
    % --------------------------------------------------------------------%

    % Convert to uint8 (0-255 range)
    grayscale_normalized_output_frame = ...
        uint8(normalized_output_frame .* 255);

    % Apply a Gaussian filter to help smooth out the final image
    if filter_output_image == true 
         grayscale_normalized_output_frame = ...
             imgaussfilt(grayscale_normalized_output_frame, ...
             3.0,"FilterSize",3); 
    end

    % Log processing time
    frame_time(frameIndex) = toc;

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
    stats.printPercentComplete(frameIndex, frame_total, frame_time(frameIndex));

end

% Close the video writer
for videosIdx = 1:length(videoWriters)
    close(videoWriters{videosIdx});
end



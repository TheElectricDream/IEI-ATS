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

%% Define processing range
% Define start and end time to process [seconds]
t_start_process = 60; 
t_end_process   = 120; 

%% Import events for inspection

% Set path to datasets
hdf5_path = ['/home/alexandercrain/Dropbox/Graduate Documents' ...
    '/Doctor of Philosophy/Thesis Research/Datasets/SPOT/HDF5/'];

% Set dataset name
%file_name = 'recording_20260127_145247.hdf5';  % Jack W. (LED Cont)
file_name = 'recording_20251029_131131.hdf5';  % EVOS - NOM - ROT
%file_name = 'recording_20251029_135047.hdf5';  % EVOS - SG - ROT
%file_name = 'recording_20251029_134602.hdf5';  % EVOS - DARK - ROT

% Load the data
tk = double(h5read([hdf5_path file_name], '/timestamp'));
xk = single(h5read([hdf5_path file_name], '/x'));
yk = single(h5read([hdf5_path file_name], '/y'));
pk = single(h5read([hdf5_path file_name], '/polarity'));

% Convert time to seconds
tk = (tk - tk(1))/1e6;

% Convert to single data type to use less memory
tk = single(tk);

% Define whether a Gaussian filter should be applied to the output
filter_image = false;

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
t_interval                  = 0.033;  % [s]
t_total                     = max(tk);  % [s]
frame_total                 = floor(t_total/t_interval);

% COHERENCE PARAMETERS
% --------------------
% Define coherence constants
r_s                         = 30/imgSz(1);  % spatial radius [pixels norm]

% Thresholds
similarity_threshold        = 0.0;
trace_threshold             = 0.0;
persistence_threshold_high  = 0.0001;
persistence_threshold_low   = 0.00001;
coherence_threshold         = 0.06;

% ADAPTIVE LOCAL TIME-SURFACE PARAMETERS
% --------------------------------
% Set adaptive local time-surface parameters
alts_params.surface_tau_min      = 0.05;
alts_params.surface_tau_max      = 0.9;
alts_params.dt                   = t_interval;
alts_params.recency_filter_size  = 9;
alts_params.recency_filter_sigma = 3.0;

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
agd_params.K       = 50000.0;  % Scaling factor (Controls "memory length")
agd_activity_store = zeros(length(frame_total),1);

%% Initialize all figure code for video output
% Indicate which videos should be saved
cohOut = false;
atsOut = true;

% Initialize the videos
[hFigs, hAxs, hImgs, videoWriters] = plot.initializeEventVideos(cohOut,...
    atsOut, imgSz);

%% Initialize data storage and perform data optimizations
% Identify number of events
current_idx     = 1;
n_events        = length(tk);
frame_time      = zeros(frame_total, 1);

% Initialize per-pixel timestamp tracking
last_event_timestamp    = zeros(imgSz);
norm_trace_map_prev     = zeros(imgSz);
time_surface_map_prev   = zeros(imgSz);
filter_mask             = ones(imgSz);

%% Data processing starting point
% Loop through the figures to capture each frame
for frameIndex = 1:frame_total  

    % Start loop timer
    tic;
    
    % Increment the interval
    t_range_c = (frameIndex - 1) * t_interval;
    t_range_n = (t_range_c+t_interval);

    % % Initialize display map for this frame
    % % We will visualize the activity state at the end of the frame
    % nunes_activity_map = G_activity; 

    % Slice the events to a valid range
    [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
    process.sliceToValidRange(t_range_n, xk, yk, tk, pk, imgSz, current_idx);

    % Confirm the presence of valid events in the packet
    % If no events are present, we skip this frame
    if isempty(t_valid)
        
        fprintf('There are no events in this slice, skipping... \n');
        continue;
        
    end   

    % ---------------------- EVENT PREPERATION------------------------%
    % ----------------------------------------------------------------%
    
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
        coherence.computeNeighborhoodStats(sorted_t, unique_idx, pos, ...
        group_ends, imgSz);

    % ---------------------- EVENT COHERENCE -----------------------------%
    % --------------------------------------------------------------------%

    % [t_mean, t_max, t_min, t_std, norm_trace_map, norm_similarity_map, ...
    % norm_persist_map, filtered_coherence_map] = ...
    % coherence.findCoherentEvents(sorted_x, sorted_y, sorted_t ,...
    % imgSz, r_s, t_interval, unique_idx, pos, group_ends,...
    % trace_threshold, similarity_threshold, persistence_threshold_high, ...
    % persistence_threshold_low, frameIndex, norm_trace_map_prev);
    % 
    % filtered_coherence_map(filtered_coherence_map<coherence_threshold) = nan;
    % 
    % % Set any retention variables
    % norm_trace_map_prev = norm_trace_map;
    % 
    % % Extract filter mask
    % filter_mask = (filtered_coherence_map>0.00);   
    % 
    % % Blur mask
    % filter_mask = imgaussfilt(filter_mask.*1, 5.0, "FilterSize", 9);

    % [norm_trace_map, norm_similarity_map, ...
    % norm_persist_map, filtered_coherence_map] = ...
    % coherence.computeCoherenceMask(sorted_x, sorted_y, sorted_t,...
    % imgSz, r_s, t_interval, unique_idx, pos, group_ends, trace_threshold, ...
    % similarity_threshold, persistence_threshold_high, persistence_threshold_low, ...
    % frameIndex, norm_trace_map_prev);
    % 
    % filtered_coherence_map(filtered_coherence_map<coherence_threshold) = nan;
    % 
    % % Set any retention variables
    % norm_trace_map_prev = norm_trace_map;
    % 
    % % Extract filter mask
    % filter_mask = (filtered_coherence_map>0.00);   
    % 
    % % Blur mask
    % filter_mask = imgaussfilt(filter_mask.*1, 5.0, "FilterSize", 9);

    % -------------- ADAPTIVE LOCAL TIME-SURFACE UPDATE ------------------%
    % --------------------------------------------------------------------%
    
    % Accumulate polarity into a 2D grid
    % If multiple events land on one pixel, we sum their polarities (e.g., +1 +1 -1 = +1)
    polarity_map = accumarray([sorted_x, sorted_y], p_signed, imgSz, @sum, 0);

    % Normalize the timestamp 
    % Note: Tried this function out of curiosity
    % https://www.mathworks.com/help/images/ref/entropyfilt.html
    t_entropy_mean = entropyfilt(t_mean);
    log_t_mean = log1p(t_entropy_mean);
    norm_t_mean_diff = log_t_mean./max(log_t_mean(:));

    [normalized_output_frame, time_surface_map, tau_filtered, decayed_surface] = ...
    accumulator.localAdaptiveTimeSurface(norm_t_mean_diff, last_event_timestamp,...
    time_surface_map_prev, alts_params, filter_mask, polarity_map);

    % Set any retention variables
    time_surface_map_prev = time_surface_map;

    % Update the last event timestamp
    last_event_timestamp = max(t_max-t_range_c, eps);

    % Normalize the timestamp 
    log_last_event_timestamp = log1p(last_event_timestamp);
    last_event_timestamp = log_last_event_timestamp./max(log_last_event_timestamp(:));

    % ----------------- NUNES GLOBAL ADAPTIVE ACCUMULATION----------------%
    % --------------------------------------------------------------------%
    
    % % Run the AGD algorithm
    % [agd_surface, agd_state, ~] = accumulator.adaptiveGlobalDecay(agd_surface,...
    %     sorted_x, sorted_y, sorted_t, agd_state, agd_params);
    % 
    % % Store the surface into the standard normalized frame
    % normalized_output_frame = agd_surface;
    % 
    % % Store activity data for later inspection
    % agd_activity_store(frameIndex) = agd_state.activity;

    % --------------------- TIME-SURFACE ACCUMULATION --------------------%
    % --------------------------------------------------------------------%
    
    % % Run the normal time-surface accumulation algorithm
    % [ts_t_map, normalized_output_frame] = accumulator.timeSurface(ts_t_map,...
    %  sorted_x, sorted_y, sorted_t, imgSz, ts_time_constant);

    % ---------- SPEED INVARIENT TIME-SURFACE ACCUMULATION ---------------%
    % --------------------------------------------------------------------%
    
    % % Run the speed invarient time-surface accumulation algorithm
    % [sits_t_map, normalized_output_frame] = ...
    %     accumulator.speedInvariantTimeSurface(sits_t_map, sorted_x,...
    %     sorted_y, sits_R);

    % ------------------------ EXPORTING VIDEO ---------------------------%
    % --------------------------------------------------------------------%

    % Convert to uint8 (0-255 range)
    grayscale_normalized_output_frame = ...
        uint8(normalized_output_frame .* 255);

    % Apply a Gaussian filter to help smooth out the final image
    if filter_image 
         grayscale_normalized_output_frame = ...
             imgaussfilt(grayscale_normalized_output_frame, ...
             3.0,"FilterSize",3); %#ok<UNRCH>
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
        writeVideo(videoWriters{1}, grayscale_normalized_output_frame'); %getframe(hFigs{1}));
        
        frameOutputFolder = '/home/alexandercrain/Videos/output_frames';
        % Also write each frame as a PNG to a folder
        % Ensure output folder exists
        if ~exist(frameOutputFolder, 'dir')
            mkdir(frameOutputFolder);
        end
        % Build filename with zero-padded frame index
        fname = fullfile(frameOutputFolder, sprintf('frame_%05d.png', frameIndex));
        % Write the PNG. imwrite expects HxWx1 or HxWx3; transpose to match display
        imwrite(grayscale_normalized_output_frame', fname);
    end

    % if cohOut
    %     % Capture the frame
    %     set(hImgs{2}, 'CData', coherence_map');
    %     set(hImgs{2}, 'AlphaData', ~isnan(coherence_map'));
    %     set(hAxs{2}, 'Visible','off');
    %     writeVideo(videoWriters{2}, getframe(hFigs{2})); 
    % end

    % Print progress
    stats.printPercentComplete(frameIndex, frame_total, frame_time(frameIndex));

end

% Close the video writer
for videosIdx = 1:length(videoWriters)
    close(videoWriters{videosIdx});
end
%% Clear out all variables
clear;
clc;
close all

%% Define general parameters

% Define start / end time to process
t_start_process = 0;        % [s]
t_end_process   = 1000;     % [s]

% Set the interval size for the event accumulation and filtering
t_interval      = 0.3;     % [s]

% Set the image size
img_size        = [640, 480]; % [pixels]

% Define if a buffered dataset is used (helpful for low RAM computers)
use_buffer      = false;  % Stores event data in a buffer

% Set a flag to save data
save_data       = true;  % Saves basic data

% Set a flag to generate figures
gen_figures     = true;  % Warning: slows down code

% Set flags to determine what kind of videos should be saved
coh_out         = false;  % If true, outputs coherence filter video 
lats_out        = true;  % If true, outputs LATS accumulator video

% Set flags to determine which filters should be run
compare_filters  = false;  % Set to true to run a second filter for comparisons
filter_selection = 'NONE';  % 'NONE', 'STC', 'BAF', 'EDF', 'STCC', 'MCF'

% Set flags to determine which accumulators should be run
compare_accumulators  = false;
accumulator_selection = 'IEI-ATS';  % 'HOTS', 'SITS', 'METS', 'LATS', 'AGD', 'EVO-ATS', 'IEI-ATS'

% Set dataset name
file_name = 'recording_20251029_131131.hdf5';  % EVOS - NOM - ROT
%file_name = 'recording_20251029_135047.hdf5';  % EVOS - SG - ROT
%file_name = 'recording_20251029_134602.hdf5';  % EVOS - DARK - ROT
%file_name = 'recording_20251029_153226.hdf5';  % EVOS - NOM - CC+ROT

% Set the plotting_frame and plotting_type 
plotting_frame               = round(205.*(0.33/t_interval));  % 205 (NOM), 179 (SG), 181 (DARK)
plotting_type                = '-NOM';  % -NOM, -SG, -DARK

%% Import events for inspection

% Import the events
[video_out_path, xk, yk, tk, pk, t_total, buf] =...
    import.importHDF5EventData(file_name, filter_selection, accumulator_selection, ...
    t_start_process, t_end_process, use_buffer);

% Calculate the total number of frames based on the total time and the time
% interval
frame_total                 = floor(t_total / t_interval);

%% Initialize video output

[h_figs, h_axs, h_imgs, video_writers] = plot.initializeEventVideos(coh_out,...
    lats_out, img_size, video_out_path);

%% Initialize comparison filters

% Initialize the Background Activity Filter by Delbruck
baf = filters.initializeDelbruckBAF(img_size, frame_total);

% Initialize the Spatiotemporal Correlation Filter by Liu
stc = filters.initializeLiuSTC(img_size, frame_total);

% Initialize the Motion Consistency Filter by Wang
mcf = filters.initializeWangMCF(img_size, frame_total);

% Initialize the Event Density Filter by Feng
edf = filters.initializeFengEDF(img_size, frame_total);

% Initialize the Space-Time-Content Correlation Filter by Li
stcc = filters.initializeLiSTCC(img_size, frame_total);

%% Initialize comparison accumulators

% Initialize the Time Surface Accumulator by Lagorce
hots = accumulator.initializeLagorceHOTS(img_size);

% Initialize the Speed Invarient Time Surface Accumulator by Manderscheid
sits = accumulator.initializeManderscheidSITS(img_size);

% Initialize the Adaptive Global Accumulator by Nunes
agd = accumulator.initializeNunesAGD(img_size, frame_total);

% Initialize the Motion Encoded Time Surface by Xu
mets = accumulator.initializeXuMETS(img_size);

% Initialize the Adaptive Time Surface by Zhu
zhu = accumulator.initializeZhuATS(img_size);

%% Initialize filter thesis research filter & accumulator parameters

% Initialize the Coherence Filter paramters
coh_params.r_s              = 60/img_size(1);  % Density kernal diameter
coh_params.s_bnd            = 0.6;             % Regularity bound 
coh_params.hpa_decay        = 0.9;             % Decay per frame for HPA calculation
coh_params.hpa_bnd          = 3;               % Number of warm-up frames before statistcs are valid

% Initialize the Adaptive Local Time Surface parameters
alts_params.dt                   = t_interval;
alts_params.buf_size             = 10;
alts_params.surface_tau_release  = 0.1;
alts_params.div_norm_exp         = 0.01;
alts_params.symmetric_tone_scale = 0.35;
alts_params.kernel_size          = 5;
alts_params.sigma_base           = 0.5;
alts_params.sigma_outlier        = 3;

%% Pre-allocate matrices and vectors for processing speed

filter_mask             = ones(img_size);
global_hot_mask         = zeros(img_size);
last_event_timestamp    = zeros(img_size);
norm_trace_map_prev     = zeros(img_size);
time_surface_map_prev   = zeros(img_size);
hot_pixel_accumulator   = zeros(img_size);
alts_activity_score.mean         = zeros(frame_total, 1);
alts_activity_score.median       = zeros(frame_total, 1);
alts_activity_score.std          = zeros(frame_total, 1);
pixel_history = NaN(prod(img_size), alts_params.buf_size);

if use_buffer == false
    n_events        = length(tk);
    frame_time      = zeros(frame_total, 1);
end

% Pre-allocate metrics storage .
frame_metrics.SRR                    = zeros(1, frame_total);
frame_metrics.ClarkEvansRemoved      = zeros(1, frame_total);
frame_metrics.ClarkEvansRemaining    = zeros(1, frame_total);
frame_metrics.ComputeTimeFilter      = zeros(1, frame_total);
frame_metrics.ComputeTimeAccumulator = zeros(1, frame_total);
frame_metrics.EventsInFrame          = zeros(1, frame_total);
frame_metrics.FilteredEvents         = zeros(1, frame_total);
frame_metrics.FilteringMEVs          = zeros(1, frame_total);
frame_metrics.AccumulatorMEVs        = zeros(1, frame_total);
frame_metrics.HotPixelCount          = zeros(1, frame_total);
frame_metrics.HotPixelThresh         = zeros(1, frame_total);

% Pre-allocate cells to store the surviving events of each chunk
global_filtered_x = cell(frame_total, 1);
global_filtered_y = cell(frame_total, 1);
global_filtered_t = cell(frame_total, 1);
global_filtered_p = cell(frame_total, 1);

%% Data processing starting point
% Loop through the figures to capture each frame

% Initialize the current index
current_idx     = 1;

fprintf('\n')
fprintf('------------------ OPTIONS -----------------\n')
fprintf('--------------------------------------------\n')
fprintf('\n')

% Start the loop, process over all frames
for frame_index = 1:frame_total 
    
    % Increment the interval
    t_range_c = (frame_index - 1) * t_interval;
    t_range_n = (t_range_c+t_interval);

    % Slice the events to a valid range
    if use_buffer == false
        [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
        process.sliceToValidRange(t_range_n, xk, yk, tk, pk,...
        current_idx);
    else
        [x_valid, y_valid, t_valid, p_valid] = buf.nextWindow(t_range_n, ...
            img_size);
    end

    % Confirm the presence of valid events in the packet
    % If no events are present, we skip this frame
    if isempty(t_valid)
        
        fprintf('There are no events in this slice, skipping... \n');
        continue;
        
    end   
    
    % Sort events using linear indexing (for computational efficiency)
    [x_sorted, y_sorted, t_sorted, p_sorted, pos, unique_idx, linear_idx,...
        group_ends, ~] = process.sortEventsByLinearIndex(x_valid, y_valid, ...
        t_valid, p_valid, img_size);

    % Reset the frames for the current loop and calculate the total event
    % count per (x,y) pixel 
    counts              = zeros(img_size);
    counts(unique_idx)  = group_ends - pos + 1;
    counts_raw          = counts;

    % Calculate the window statistics
    [t_mean, ~, t_max, t_min, t_mean_diff, t_std_diff, stats_time] = ...
        stats.computeNeighborhoodStats(t_sorted - min(t_sorted), ...
        unique_idx, pos, group_ends, img_size);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    % =========================== FILTERING ==============================%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    
    % Start filter timer
    filtering_start = tic;

    % Print the current chosen filter
    if frame_index == 1
        fprintf(['FILTER SELECTION: ' filter_selection '\n']);
    end
    
    % If comparing to an existing filter, run said filtering algorithm
    if compare_filters
        [filter_mask_comparison, stc, baf, edf, stcc, mcf] = ...
            filters.applyChosenEventFilter(filter_selection, stc, baf, edf,...
            stcc, mcf, x_sorted, y_sorted, t_sorted, p_sorted, counts); %#ok<UNRCH>
    end
    
    % Run the Coherence Filtering algorithm
    [norm_trace_map, norm_trace_map_nofilt, norm_regularity_map, norm_persist_map,...
        norm_persist_map_raw, filtered_coherence_map, hot_pixel_accumulator, aperiodic_mask,...
        local_hot_mask, global_hot_mask, secondary_cleaning, filtered_counts_mask] = filters.computeCoherenceMask(x_sorted,...
        y_sorted, t_sorted, img_size, t_interval, coh_params, frame_index, ...
        norm_trace_map_prev, t_std_diff,...
        t_mean_diff, counts, hot_pixel_accumulator,...
        gen_figures, global_hot_mask);

    % Set any retention variables
    norm_trace_map_prev = norm_trace_map;
    filtered_coherence_map_raw = filtered_coherence_map;
    
    % Create the filter mask
    filter_mask = filtered_coherence_map;
    filter_mask(isnan(filter_mask)) = 0;

    % Apply the Gaussian blur to remove last few points
    filter_mask = single(imgaussfilt(single(filter_mask), 5.0, "FilterSize", 9));
    
    % Calculate the Rosin threshold for the filter mask
    [ros_threshold, threshold_diagnostic] = stats.rosinThreshold(filter_mask);
    
    % Use the threshold on the mask & log the result
    filtered_coherence_map_nogauss = filtered_coherence_map.*filter_mask;
    filter_mask(filter_mask < ros_threshold) = 0;
    filtered_coherence_map = filtered_coherence_map.*filter_mask;

    % Compute event-level pass/total counts for the COH filter.
    % The other filters return these directly; COH needs them
    % derived from the pixel-level mask and event locations.
    coh_linear_idx  = sub2ind(img_size, x_sorted(:), y_sorted(:));
    coh_n_total     = numel(x_sorted);
    coh_n_passed    = sum(filter_mask(coh_linear_idx) > 0);

    % Grab the indices of the filtered mask
    filter_mask_idx = find(filter_mask>0);
   
    % Remove the filtered results from the statistics
    counts          = counts.*(filter_mask>0);
    t_mean          = t_mean.*(filter_mask>0); %#ok<NASGU>
    t_mean_diff     = t_mean_diff.*(filter_mask>0); %#ok<NASGU>

    % If desired, generate plots for a journal article
    if frame_index == plotting_frame && gen_figures
        
        journal.generateAllJournalPlots(t_mean_diff, norm_trace_map, norm_regularity_map, ...
            norm_persist_map, filtered_coherence_map, x_sorted, y_sorted, t_sorted, coh_params, ...
            local_hot_mask, norm_trace_map_nofilt, filtered_coherence_map_raw, filtered_coherence_map_nogauss, ...
            filter_mask, threshold_diagnostic, ros_threshold, plotting_type) %#ok<UNRCH>
    end

    % Stop filter timer
    filtering_stop = toc(filtering_start);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    % ========================= ACCUMULATION =============================%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    
    % Start accumulator timer
    accumulator_start = tic;

    if frame_index == 1
        fprintf(['ACCUMULATOR SELECTION: ' accumulator_selection '\n']);
    end

    filtered_x = x_sorted;
    filtered_y = y_sorted;
    filtered_t = t_sorted;
    filtered_p = p_sorted;
    
    % Strip the defective pixels globally
    hot_pixel_idx = find(local_hot_mask);
    
    if ~isempty(hot_pixel_idx)
        % Strip from the coordinate arrays
        sorted_lin_idx = sub2ind(img_size, filtered_x, filtered_y);
        remove_mask = ismember(sorted_lin_idx, hot_pixel_idx);
        
        filtered_x(remove_mask) = [];
        filtered_y(remove_mask) = [];
        filtered_t(remove_mask) = [];
        filtered_p(remove_mask) = [];
        
    end

    if frame_index == plotting_frame && gen_figures
        journal.showScatterPlotOfEventVector(filtered_x, filtered_y,...
            filtered_t, ['Hot-Pixel-Only-Filtered-Event-Data-Nominal-Rot'...
            plotting_type '.pdf'], false) %#ok<UNRCH>
    end

    % Any events which fall within the filter mask should not be included
    % in the accumulation
    event_linear_idx = sub2ind(img_size, filtered_x, filtered_y);
    remove_mask = filter_mask(event_linear_idx) == 0;

    % Apply to all event vectors simultaneously
    filtered_x(remove_mask) = [];
    filtered_y(remove_mask) = [];
    filtered_t(remove_mask) = [];
    filtered_p(remove_mask) = [];

    % Ensure polarity is -1 and 1 (if it's 0 and 1)
    p_signed = double(filtered_p);
    p_signed(p_signed == 0) = -1;

    % Drop the current chunk's surviving events into the global storage
    global_filtered_x{frame_index} = filtered_x;
    global_filtered_y{frame_index} = filtered_y;
    global_filtered_t{frame_index} = filtered_t;
    global_filtered_p{frame_index} = filtered_p;
    
    if frame_index == plotting_frame && gen_figures
        journal.showScatterPlotOfEventVector(filtered_x, filtered_y, ...
            filtered_t, ['Fully-Filtered-Event-Data-Nominal-Rot' ...
            plotting_type '.pdf'], false); %#ok<UNRCH>
    end

    if compare_accumulators
        [agd, hots, sits, mets, zhu] = ...
            accumulators.applyChosenEventAccumulator(accumulator_selection,...
            agd, hots, sits, mets, zhu) %#ok<UNRCH>
    end

    % ------------ ADAPTIVE LOCAL TIME-SURFACE UPDATE ----------------%
    % ----------------------------------------------------------------%
      
    [t_mean, t_std, t_mean_diff, t_std_diff, ~, pixel_history] = ...
        stats.backfillIEIStatisticsReference(t_sorted, unique_idx, pos, ...
        group_ends, img_size, counts, alts_params.buf_size, pixel_history,...
        alts_params.kernel_size);

    % Overwrite the stats
    t_mean = t_mean.*(filter_mask>0);
    t_mean_diff = t_mean_diff.*(filter_mask>0);

    if frame_index == plotting_frame && gen_figures
        journal.showScatterPlotOfRuleMaps3D(t_mean_diff,...
            ['Fully-Filtered-Backfilled-Mean-IEI' plotting_type '.pdf'],...
            true); %#ok<UNRCH>
    end

    if frame_index == 86
        fprintf('Stop')
    end

    % Accumulate polarity into a 2D grid
    % If multiple events land on one pixel, we sum their polarities 
    polarity_map = accumarray([filtered_x, filtered_y],...
        p_signed, img_size, @sum, 0);

    % Accumulate the frames
    lats_outputs = accumulator.localAdaptiveTimeSurface(t_mean_diff,...
        time_surface_map_prev, alts_params, filter_mask, polarity_map, counts);
    
    % Extract variables from the LATS outputs
    normalized_output_frame = lats_outputs.normalized_output_frame;
    time_surface_map_raw = lats_outputs.time_surface_map_raw;
    tau_filtered = lats_outputs.tau_filtered;
    adaptive_gains = lats_outputs.adaptive_gains;

    % Set any retention variables
    time_surface_map_prev = time_surface_map_raw;

    % Store the adaptive map score
    alts_activity_score.mean(frame_index) = mean(adaptive_gains...
        (abs(adaptive_gains)>0));
    alts_activity_score.median(frame_index) = median(adaptive_gains...
        (abs(adaptive_gains)>0));
    alts_activity_score.std(frame_index) = std(adaptive_gains...
        (abs(adaptive_gains)>0));

    % Start accumulator timer
    accumulator_stop = toc(accumulator_start);

    if frame_index == 1
        fprintf('\n--------------------------------------------\n\n');
    end

    % ========================== FRAME METRICS ===========================%
    % --------------------------------------------------------------------%
    
    % Determine n_passed and n_total from whichever filter ran.
    % Each filter stores these differently; unify them here.
    recovered_events = process.recoverRemovedEvents(filter_mask, counts_raw);
    switch filter_selection
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
            temporary_filter_mask = recovered_events.*counts_raw;
            metrics_background_map = (temporary_filter_mask>0).*1.0;
            metrics_foreground_map = (filter_mask>0).*1.0;
        case 'MCF'
            metrics_n_passed = mcf_n_pass;
            metrics_n_total  = mcf_n_tot;
            metrics_background_map = recovered_events;
            metrics_foreground_map = filter_mask;
        case 'NONE'
            metrics_n_passed = numel(x_sorted);
            metrics_n_total  = numel(x_sorted);
            metrics_background_map = (counts>0).*1.0;
            metrics_foreground_map = (counts>0).*1.0;
        otherwise
            metrics_n_passed = numel(x_sorted);
            metrics_n_total  = numel(x_sorted);
            metrics_background_map = (counts>0).*1.0;
            metrics_foreground_map = (counts>0).*1.0;
    end

    % Store all important data for future plotting
    if save_data
        frame_metrics.HotPixelCount(frame_index) ...
            = sum(local_hot_mask(:));
        frame_metrics.FilterThreshold(frame_index)...
            = ros_threshold;
        frame_metrics.ElbowDiagnostics{frame_index}...
            = threshold_diagnostic;
        frame_metrics.SRR(frame_index) = ...
            metrics.computeSignalRetentionRate(metrics_n_passed, metrics_n_total);
        frame_metrics.ClarkEvansRemoved(frame_index) = ...
            metrics.computeClarkEvans(metrics_background_map);
        frame_metrics.ClarkEvansRemaining(frame_index) = ...
            metrics.computeClarkEvans(metrics_foreground_map);
        frame_metrics.ComputeTimeFilter(frame_index) = filtering_stop;
        frame_metrics.ComputeTimeAccumulator(frame_index) = accumulator_stop;
        frame_metrics.EventsInFrame(frame_index) = metrics_n_total;
        frame_metrics.FilteredEvents(frame_index) = ...
            metrics_n_total-metrics_n_passed;
        frame_metrics.FilteringMEVs(frame_index) = ...
            metrics_n_total/filtering_stop;
        frame_metrics.AccumulatorMEVs(frame_index) = ...
            metrics_n_total/accumulator_stop;
    end

    % ------------------------ EXPORTING VIDEO ---------------------------%
    % --------------------------------------------------------------------%

    % Convert to uint8 (0-255 range)
    grayscale_normalized_output_frame = ...
        uint8(normalized_output_frame .* 255);

    if lats_out
        % Capture the frame for the video writer
        set(h_imgs{1}, 'CData', grayscale_normalized_output_frame');
        set(h_imgs{1}, 'AlphaData', ~isnan(grayscale_normalized_output_frame'));
        set(h_axs{1}, 'Visible','off');
        colormap(gray);
        clim([0 255]);
        writeVideo(video_writers{1}, grayscale_normalized_output_frame'); 
        
        frameOutputFolder =  string(video_out_path(1:end-4)) + "-FRAMES";
        % Also write each frame as a PNG to a folder
        % Ensure output folder exists
        if frame_index == 1
            if ~exist(frameOutputFolder, 'dir')
                mkdir(frameOutputFolder);
            else 
                rmdir(frameOutputFolder, 's');
                mkdir(frameOutputFolder); % Create the output folder
            end
        end
        % Build filename with zero-padded frame index
        fname = fullfile(frameOutputFolder, sprintf('frame_%05d.png', frame_index));
        imwrite(grayscale_normalized_output_frame', fname);
    end

    % Print progress
    stats.printPercentComplete(frame_index, frame_total, ...
        frame_metrics.ComputeTimeFilter(frame_index), ...
        frame_metrics.ComputeTimeAccumulator(frame_index));

    % Create a distinct filename for the chunk and save it
    dataOutputFolder =  string(video_out_path(1:end-4)) + "-METRICS";
    if frame_index == 1
        if ~exist(dataOutputFolder, 'dir')
            mkdir(dataOutputFolder);
        else 
            rmdir(dataOutputFolder, 's');
            mkdir(dataOutputFolder); % Create the output folder
        end
    end

    % Check if the chunk is full, OR if it's the very last frame of the dataset
    if (frame_index == frame_total) && save_data

        filename = fullfile(dataOutputFolder, 'metrics.mat');
        save(filename, 'frame_metrics');
        fprintf('\n[IO] Saved metrics to %s\n', filename);

    end

end

if gen_figures
    % Exporting figures
    journal.showALTSActivityScoreStatistics(alts_activity_score, ...
        ['ALTS-Adaptive-Gain-Mean-Rot-Nom' plotting_type],false) %#ok<UNRCH>
    theta_series = frame_metrics.FilterThreshold;   
    srr_series   = frame_metrics.SRR;  
    srr_series(srr_series == 0) = [];
    journal.showRosinThresholdStability(theta_series, srr_series, ...
        ['Rosin-Sequence-Threshold-Stability' plotting_type '.pdf'], false);
end

% Close the video writer
for videosIdx = 1:length(video_writers)
    close(video_writers{videosIdx});
end

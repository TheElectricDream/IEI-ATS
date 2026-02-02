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

% Load the data
tk = double(h5read([hdf5_path file_name], '/timestamp'));
xk = single(h5read([hdf5_path file_name], '/x'));
yk = single(h5read([hdf5_path file_name], '/y'));
pk = single(h5read([hdf5_path file_name], '/polarity'));

% Convert time to seconds
tk = (tk - tk(1))/1e6;

% Convert to single data type to use less memory
tk = single(tk);

% Define whether polarity should be ignored
include_polarity = false;

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
% Set the image size
imgSz = [640, 480]; 

% Define coherence constants
r_s = 30/imgSz(1);  % spatial radius [pixels norm]

% Set persistance variables
gamma_trace = 0.5;

% Thresholds
similarity_threshold = 0;
trace_threshold = 0;
persistance_threshold = 0;

% Set coherence weights
coherence_s = 1;
coherence_p = 1;
coherence_t = 1;

% Set the coherence threshold
coherence_threshold = 0.01;

% Set time-surface parameters
surface_k_tau = 1.0;
surface_tau_min = 0.1;
surface_tau_max = 1.2;

% Set recency time constant
recency_T = 0.5;

% Set the time interval to accumulate over
t_interval = 0.01;  % [s]
t_total = max(tk);  % [s]
frame_total = floor(t_total/t_interval);

%% Initialize all figure code for video output
% Indicate which videos should be saved
cohOut = false;
atsOut = true;

% Initialize the videos
[hFigs, hAxs, hImgs, videoWriters] = plotting.initializeEventVideos(cohOut, atsOut, imgSz);

%% Initialize data storage and perform data optimization before algorithm
% Identify number of events
current_idx = 1;
n_events = length(tk);
frame_time = zeros(frame_total, 1);

% Initialize per-pixel timestamp tracking
last_event_timestamp = zeros(imgSz);
adaptive_surface_prev = zeros(imgSz);

% Initialize persistent variables for Temporal Coherence
prev_cloud_xy = [];
prev_cloud_vals = [];

% Keep a running log of the max value and always scale by that number to
% generate consistent frames
max_time_surface_global = eps;

%% Data processing starting point
% Loop through the figures to capture each frame
for frameIndex = 1:frame_total  

    % Start loop timer
    tic;
    
    % Increment the interval
    t_range_c = (frameIndex - 1) * t_interval;
    t_range_n = (t_range_c+t_interval);

    % Slice the events to a valid range
    [current_idx, x_valid, y_valid, t_valid] = ...
    process.sliceToValidRange(t_range_n, xk, yk, tk, imgSz, current_idx);

    if ~isempty(t_valid)
        % Convert 2D subscripts (x,y) to 1D linear indices
        % Imagine you are a post-man with a disorganized stack of letters. 
        % Instead of dealing with letters for
        % 3rd Avenue, 5th street, the sub2ind function assigns and "ID" for
        % each house. So going forward, (3,5) might just be "House #1".
        % Programatically this just means that (1,1) is "1", (1,2) is "2".
        % So you will have a list at the end which is of size x*y. In this
        % case with a frame of size 480 by 640, the MAXIMUM size of the list
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

        % Reset the frames for the current loop
        t_diff = nan(imgSz);
        t_max = nan(imgSz);
        t_min = nan(imgSz);
        t_std = nan(imgSz);
        sum_exp_dist_map = zeros(imgSz);
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

        % ---------------------- EVENT COHERENCE -------------------------%
        % ----------------------------------------------------------------%

        % Now we want to start implementing some coherence rules from one
        % window to the next. What does that mean. It means that we start
        % look at frameIndex and pick an event within that current list. We
        % then calculate some parameters using the closest events. 

        % First we find spatial neighbours in the same plane 
        [neighbours_db, distances_db] = coherence.findSpatialNeighbours(...
            sorted_x, sorted_y, sorted_t, r_s, imgSz, t_interval);

        % Calculate the sum of the cells
        sum_exp_event = cellfun(@sum, distances_db);
        
        % Finally, we can extract each chunk sequentially and calculate the
        % maximum and minimum of that "chunk".
        for k = 1:length(unique_idx)
            val_chunk_t = sorted_t(pos(k):group_ends(k));
            val_chunk_exp = sum_exp_event(pos(k):group_ends(k));
            idx = unique_idx(k);
            t_diff(idx) = mean(diff(val_chunk_t));
            t_max(idx) = max(val_chunk_t);
            t_min(idx) = min(val_chunk_t);
            t_std(idx) = std(val_chunk_t);
            sum_exp_dist_map(idx) = max(val_chunk_exp);
        end
        
        % Normalize the trace map by calculating the LOG first, then
        % normalizing it
        log_trace_map = log1p(sum_exp_dist_map'); 
        norm_trace_map = log_trace_map' ./ max(log_trace_map(:));
       
        % Calculate the actual mean (range / (N-1)). This is the average of
        % the current map, but if we want to get the moving average accross
        % % the frames we will need to do a bit more work
        % t_max(isnan(t_max))=0;
        % t_min(isnan(t_min))=0;
        % tk_diff_map = (t_max - t_min) ./ (counts - 1);
        tk_diff_map = t_diff;
        tk_diff_map(isinf(tk_diff_map)) = NaN;
        tk_diff_map(isnan(tk_diff_map)) = 0.0;

        % Reset statistics maps
        cv_map = zeros(imgSz);
        cv_values = zeros(imgSz);
        regularity_map = zeros(imgSz);
        
        % Calculate the coefficient of variation map
        mask = tk_diff_map ~= 0;
        cv_values(mask) = t_std(mask) ./ tk_diff_map(mask);
        cv_map = max(cv_map, cv_values);

        % Also calculate the inverse as this is more representative of
        % similarity
        regularity_map(mask) = 1 ./ (cv_values(mask) + 1); 
        regularity_map(regularity_map==1)=0; 

    else
        tk_diff_map = nan(imgSz);

    end   

    % Discretize the event array into voxels so that convolutions can
    % be applied
    opts_voxel.dx = 1.0;
    opts_voxel.dy = 1.0;
    opts_voxel.dt = t_interval;
    [Vocc, xEdges, yEdges, tEdges] = voxelization.discretizeEventsToVoxels(sorted_x, ...
        sorted_y, sorted_t, opts_voxel);

    % 3D convolution on the voxel volume
    h = ones(3,3,3,'single')/8;                 
    Vsm = convn(single(Vocc), h, 'same');        

    % Remove disconnected regions
    cleanVsm = bwareaopen(Vsm, 2);

    % Create a clean map
    cleanMap = (Vsm(:,:,1) + Vsm(:,:,2));

    % Keep above 0.5
    cleanMap(cleanMap<0.5)=0;
    cleanMap(cleanMap>=0.5)=1;
    
    % Convert the map to a mesh for additional processing
    [normalized_cv_point_cloud] = process.generateMeshFromFrame(regularity_map');

    % Calculate point to point similarity in the CV map
    [similarity_score] = coherence.findSimilarities(sorted_x,...
        sorted_y, sorted_t./t_interval, imgSz);

    % Initialize the empty map
    similarity_map = nan(imgSz(2), imgSz(1)); 

    % Assign values directly using linear indexing
    linear_idx_cv_map = sub2ind(size(similarity_map), sorted_y,...
        sorted_x);
    similarity_map(linear_idx_cv_map) = similarity_score;

    % Log normalize the similarity map
    log_similarity_map = log1p(similarity_map-min(similarity_map(:))); 
    norm_similarity_map = log_similarity_map' ./ max(log_similarity_map(:));

    % Reset the background to zero for visualization purposes
    norm_similarity_map(isnan(norm_similarity_map)) = 0;

    % Calculate the persistance map
    [rowsExists, colsExists] = find(norm_trace_map>0);
    logPersistance = zeros(length(rowsExists),1);
    persist_map = zeros(size(norm_trace_map));
    
    if frameIndex == 1
    else
        % for numValidIdx = 1:length(colsExists)
        %     [closestRows, closestCols, sortedIdx,...
        %         t_row, t_col, t_val, ...
        %         b_rows, b_cols, b_vals, minDist] = coherence.findPersistence(norm_trace_map, ...
        %         norm_trace_map_prev, imgSz, rowsExists(numValidIdx), colsExists(numValidIdx));
        %     % plotting.visualizePersistance(t_row, t_col, t_val, ...
        %     %             b_rows, b_cols, b_vals, sortedIdx);
        %     logPersistance(numValidIdx) = minDist;
        % end
        % 
        % % Map each minDist to its (row, col) in the new map
        % persist_map(sub2ind(size(persist_map), ...
        % rowsExists, colsExists)) = logPersistance;
        rowsExists = rowsExists(:);
        colsExists = colsExists(:);
        
        logPersistance = coherence.findPersistenceVectorized( ...
            norm_trace_map, norm_trace_map_prev, imgSz, rowsExists, colsExists);

        persist_map(sub2ind(size(persist_map), rowsExists, colsExists)) = logPersistance;
    end

    norm_trace_map_prev = norm_trace_map;

    % Invert persistance
    inverted_persistance = 1./(persist_map + 1);
    inverted_persistance(inverted_persistance==1)=0; 

    % Log normalize the persistance map
    log_persist_map = log1p(inverted_persistance); 
    norm_persist_map = log_persist_map ./ max(log_persist_map(:));    

    % % Calculate the coherence map
    % coherence_map = ((coherence_s .* norm_similarity_map) + ... 
    %             (coherence_t .* norm_trace_map));

    % Use the coherence map to filter out the values
    filtered_similarity_map = zeros(size(norm_similarity_map));
    filtered_trace_map = zeros(size(norm_trace_map));
    filtered_persistance_map = zeros(size(norm_persist_map));

    similarity_mask = (norm_similarity_map >= similarity_threshold);
    trace_mask = (norm_trace_map >= trace_threshold);
    persistance_mask = (norm_persist_map >= persist_map);

    filtered_similarity_map(similarity_mask) = norm_similarity_map(similarity_mask);
    filtered_trace_map(trace_mask) = norm_trace_map(trace_mask);
    filtered_persistance_map(persistance_mask) = norm_persist_map(persistance_mask);

    filtered_coherence_map = ((coherence_s .* filtered_similarity_map) + ... 
                (coherence_t .* filtered_trace_map) +...
                coherence_p .* filtered_persistance_map);

    filtered_coherence_map(filtered_persistance_map>0.999 & filtered_trace_map < 0.1)=0;
    filtered_coherence_map(filtered_coherence_map<=1.65)=0;

    % ----------------- ADAPTIVE TIME SURFACE UPDATE ---------------------%
    % --------------------------------------------------------------------%

    % Calculate the candidate decay time
    candidate_tau_n = surface_k_tau.*max(tk_diff_map, eps).*cleanMap;

    % Bound the candidate map
    candidate_tau_n(candidate_tau_n>surface_tau_max) = surface_tau_max;
    candidate_tau_n(candidate_tau_n<surface_tau_min) = surface_tau_min;

    % Calculate recency 
    recency_function = zeros(size(candidate_tau_n));
    recency_function(isnan(recency_function)) = 0;
    recency_function(last_event_timestamp~=0) = ...
        exp(-last_event_timestamp(last_event_timestamp~=0) ./ recency_T);
    recency_weighted_tau = candidate_tau_n.*recency_function +...
        surface_tau_min.*(1 - recency_function);
    
    % Accumulate the time surface
    if frameIndex == 1
        time_surface_map = filtered_coherence_map;
        decayed_surface = filtered_coherence_map;
        tau_current = recency_weighted_tau;
        norm_log_decayed_surface = filtered_coherence_map;
    else
        beta_tau = 0.1;

        tau_current = (1 - beta_tau) .* tau_map_prev + beta_tau .* candidate_tau_n;

        decayed_surface = tau_current .*...
              exp(-t_interval./tau_current);

        time_surface_map = tau_current;
    end
    time_surface_map_prev = time_surface_map;
    tau_map_prev = tau_current;

    % Update the last event timestamp
    last_event_timestamp = t_max.*cleanMap;

    % Normalize the time surface map
    log_time_surface_map = log1p(time_surface_map);
    norm_log_time_surface_map = log_time_surface_map ./ max(log_time_surface_map(:));
    min_log_time_surface_map = min(norm_log_time_surface_map(:));
    norm_log_time_surface_map = norm_log_time_surface_map-min_log_time_surface_map;
    
    % Convert to uint8 (0-255 range)
    mapdatauint8 = uint8(norm_log_time_surface_map .* 255);
    %mapdatauint8_gauss = imgaussfilt(mapdatauint8, 3.0,"FilterSize",3);

    % Apply a Gaussian filter to help smooth out the final image
    % if filter_image 
    %     mapdatauint8_gauss = imgaussfilt(mapdatauint8, 3.0,"FilterSize",3);
    % end

    % ------------------------ EXPORTING VIDEO ---------------------------%
    % --------------------------------------------------------------------%

    % Log processing time
    frame_time(frameIndex) = toc;

    if atsOut
        % Capture the frame
        set(hImgs{1}, 'CData', mapdatauint8');
        set(hImgs{1}, 'AlphaData', ~isnan(mapdatauint8'));
        set(hAxs{1}, 'Visible','off');
        colormap(gray);
        clim([0 255]);
        writeVideo(videoWriters{1}, mapdatauint8'); %getframe(hFigs{1}));
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
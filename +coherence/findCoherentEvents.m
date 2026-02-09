function [t_mean, t_max, t_min, t_std, norm_trace_map, norm_similarity_map, ...
    norm_persist_map, filtered_coherence_map] = ...
    findCoherentEvents(sorted_x, sorted_y, sorted_t ,...
    imgSz, r_s, t_interval, unique_idx, pos, group_ends, trace_threshold, ...
    similarity_threshold, persistence_threshold_high, persistence_threshold_low, ...
    frameIndex, norm_trace_map_prev)

    % Reset the frames for the current loop
    t_mean = zeros(imgSz);
    t_max = zeros(imgSz);
    t_min = zeros(imgSz);
    t_std = zeros(imgSz);
    sum_exp_dist_map = zeros(imgSz);

   if ~isempty(sorted_t)

        % Now we want to start implementing some coherence rules from one
        % window to the next. What does that mean. It means that we start
        % look at frameIndex and pick an event within that current list. We
        % then calculate some parameters using the closest events. 

        % First we find spatial neighbours in the same plane 
        [~, distances_db] = coherence.findSpatialNeighbours(...
            sorted_x, sorted_y, sorted_t, r_s, imgSz, t_interval);

        % Calculate the sum of the cells
        sum_exp_event = cellfun(@sum, distances_db);

        % Finally, we can extract each chunk sequentially and calculate the
        % maximum and minimum of that "chunk".
        for k = 1:length(unique_idx)
            val_chunk_t = sorted_t(pos(k):group_ends(k));
            val_chunk_exp = sum_exp_event(pos(k):group_ends(k));
            idx = unique_idx(k);
            t_mean(idx) = mean(diff(val_chunk_t));
            t_max(idx) = max(val_chunk_t);
            t_min(idx) = min(val_chunk_t);
            t_std(idx) = std(diff(val_chunk_t));
            sum_exp_dist_map(idx) = max(val_chunk_exp);
        end

        % Remove events which are not consistant enough in space
        trace_mask = (sum_exp_dist_map <= trace_threshold);
        sum_exp_dist_map(trace_mask) = 0;

        % Normalize the trace map by calculating the LOG first, then
        % normalizing it
        log_trace_map = log1p(sum_exp_dist_map'); 
        norm_trace_map = log_trace_map' ./ max(log_trace_map(:));

        % % Reset statistics maps
        % cv_map = zeros(imgSz);
        % cv_values = zeros(imgSz);
        % regularity_map = zeros(imgSz);
        % 
        % % Calculate the coefficient of variation map
        % mask = t_mean ~= 0;
        % cv_values(mask) = t_std(mask) ./ t_mean(mask);
        % cv_map = max(cv_map, cv_values);
        % 
        % % Also calculate the inverse as this is more representative of
        % % similarity
        % regularity_map(mask) = 1 ./ (cv_values(mask) + 1); 
        % regularity_map(regularity_map==1)=0; 

    else

        t_mean = nan(imgSz);

    end   

    % Calculate point to point similarity in the CV map
    [similarity_score] = coherence.findSimilarities(sorted_x,...
        sorted_y, sorted_t./t_interval, imgSz);

    % Initialize the empty map
    similarity_map = nan(imgSz(2), imgSz(1)); 

    % Assign values directly using linear indexing
    linear_idx_cv_map = sub2ind(size(similarity_map), sorted_y,...
        sorted_x);
    similarity_map(linear_idx_cv_map) = similarity_score;

    % Remove events which do not meet the similarity criteria
    similarity_mask = (similarity_map <= similarity_threshold);
    similarity_map(similarity_mask) = 0;

    % Log normalize the similarity map
    log_similarity_map = log1p(similarity_map); 
    norm_similarity_map = log_similarity_map' ./ max(log_similarity_map(:));

    % Reset the background to zero for visualization purposes
    norm_similarity_map(isnan(norm_similarity_map)) = 0;

    if frameIndex == 1
        persist_map = norm_trace_map;
    else
        % Reset maps
        persist_map = zeros(size(norm_trace_map));

        % Assess persistence across frames
        [~, ~, minDists, validIdx] = ...
            coherence.findPersistenceVectorized(norm_trace_map,... 
            norm_trace_map_prev, imgSz);

        % Assign results directly to the map
        if ~isempty(validIdx)
            persist_map(validIdx) = minDists;
        end
    end

    if frameIndex ~= 1
        persist_map(persist_map<persistence_threshold_low)=nan;
        persist_map(persist_map>persistence_threshold_high)=nan;
    end
    
    % % Invert persistence
    % inverted_persistence = 1./(persist_map + 1);
    % inverted_persistence(inverted_persistence==1)=0; 

    % Log normalize the persistence map
    log_persist_map = log1p(persist_map); 
    norm_persist_map = log_persist_map ./ max(log_persist_map(:));    

    % % Calculate the coherence map
    % coherence_map = ((coherence_s .* norm_similarity_map) + ... 
    %             (coherence_t .* norm_trace_map));

    % filtered_coherence_map = ((coherence_s .* norm_similarity_map .* 0) .* ... 
    %             (coherence_t .* norm_trace_map) .*...
    %             coherence_p .* norm_persist_map);

    opts.dx = 1;
    opts.dy = 1;
    opts.dt = t_interval;
    [voxelOccupancy, ~, ~,...
    ~] = voxelization.discretizeEventsToVoxels(sorted_x, sorted_y, sorted_t, opts);
    voxelCleanMap = bwareaopen(voxelOccupancy, 10);


    cleanMap = voxelCleanMap(:,:,1)+voxelCleanMap(:,:,2);

    filtered_coherence_map = (norm_trace_map .* norm_persist_map).*cleanMap;

end
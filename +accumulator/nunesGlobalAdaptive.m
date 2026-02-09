function [normalized_output_frame, state] = nunesGlobalAdaptive(xk, yk, pk, imgSz, state, params)
%ADAPTIVEGLOBALDECAY Implements the Nunes et al. (2023) adaptive decay algorithm.
%
%   This function accumulates events into a 2D surface using a decay factor
%   that depends on the global event activity (Global Event Count).
%
%   Inputs:
%       xk, yk  - Vectors of x and y coordinates (1-based).
%       pk      - Vector of polarities (typically -1 and 1, or 0 and 1).
%       imgSz   - Size of the image [height, width].
%       state   - Struct containing the persistent state of the algorithm:
%                   .S       : The current accumulated surface (H x W).
%                   .N_last  : The global event count at the last update for each pixel (H x W).
%                   .N_global: The running total of global events processed.
%       params  - Struct containing algorithm parameters:
%                   .tau_N   : The decay constant in "number of events" domain.
%                              (e.g., 5000 - 50000 depending on resolution).
%
%   Outputs:
%       normalized_output_frame - The visualization frame (0 to 1).
%       state                   - Updated state struct.
    
    % Ensure polarity values are providing contrast
    pk(pk == 0) = -1;

    % Initialize State if empty
    if isempty(state)
        state.S = zeros(imgSz);
        state.N_last = zeros(imgSz); % Initialize with 0 count
        state.N_global = 0;
    end
    
    % Ensure inputs are column vectors
    xk = double(xk(:));
    yk = double(yk(:));
    pk = double(pk(:));
    
    num_events = length(xk);
    if num_events == 0
        % If no events, just decay the existing surface to current global count
        % (Global count doesn't change, so typically no decay happens, 
        % or we could decay if we considered time. But in this algorithm 
        % decay is driven by NEW events. So we return current state.)
        normalized_output_frame = normalizeSurface(state.S);
        return;
    end
    
    % Algorithm Parameters
    tau_N = params.tau_N; % Decay constant in "Event Number" domain
    
    % Vectorized Processing
    % To process efficiently in MATLAB, we group events by pixel index.
    % However, strict sequential order matters for the global counter.
    % For a batch implementation, we can calculate the final value for each
    % pixel analytically based on the events that fell on it in this batch.
    
    % Global Event Indices for the current batch
    % The global counter increments by 1 for every event in the stream.
    global_indices = state.N_global + (1:num_events)';
    
    % Linear Pixel Indices
    linear_indices = sub2ind(imgSz, xk, yk);
   
    % We need to accumulate the effect of events for each unique pixel.
    % Since accumulation is recursive: S_new = S_old * decay + polarity,
    % we must process events at the same pixel in temporal order.
    
    % Current Surface and Last Count maps
    S = state.S;
    N_last = state.N_last;
    
    % Identify all events happening at each pixel
    % Sort by pixel index, but maintain temporal order (stable sort)
    [sorted_pix_idx, sort_order] = sort(linear_indices);
    sorted_global_indices = global_indices(sort_order);
    sorted_pk = pk(sort_order);
    
    % Find boundaries where pixel index changes
    % diff(sorted_pix_idx) != 0
    boundaries = [0; find(diff(sorted_pix_idx) ~= 0); length(sorted_pix_idx)];
    
    % Iterate over pixels (implicitly)
    % We will update S and N_last for all touched pixels.
    % Since we can't easily vectorise the variable-length recurrence per pixel
    % without a Mex file, we use the analytical sum for the batch.
    
    % Get previous state for all affected pixels
    unique_pix_indices = sorted_pix_idx(boundaries(2:end));
    S_prev = S(unique_pix_indices);
    N_last_prev = N_last(unique_pix_indices);
    
    % Since MATLAB loops are slow, we iterate through the 'boundaries'
    % This loop length is == number of active pixels (e.g. 5% of image), 
    % which is much smaller than number of events.
    
    S_new_values = zeros(size(unique_pix_indices));
    N_last_new_values = zeros(size(unique_pix_indices));
    
    for i = 1:length(unique_pix_indices)
        start_ptr = boundaries(i) + 1;
        end_ptr = boundaries(i+1);
        
        % Events for this specific pixel
        k_seq = sorted_global_indices(start_ptr:end_ptr);
        p_seq = sorted_pk(start_ptr:end_ptr);
        
        % Previous stats
        n_prev = N_last_prev(i);
        val = S_prev(i);
        k_end = k_seq(end);
        
        % Term 1: Decay old value to current time (k_end)
        decay_old = exp(-(k_end - n_prev) / tau_N);
        term1 = val * decay_old;
        
        % Term 2: Sum of new events decayed to current time (k_end)
        % weights = exp(-(k_end - k_seq) / tau_N)
        weights = exp(-(k_end - k_seq) / tau_N);
        term2 = sum(p_seq .* weights);
        
        S_new_values(i) = term1 + term2;
        N_last_new_values(i) = k_end;
    end
    
    % Update State
    S(unique_pix_indices) = S_new_values;
    N_last(unique_pix_indices) = N_last_new_values;
    state.S = S;
    state.N_last = N_last;
    state.N_global = state.N_global + num_events;
    
    % Prepare Output Frame
    % To visualize, we must decay ALL pixels to the current global time.
    % Otherwise, pixels not updated in this batch look "stuck" in the past.
    
    current_global_time = state.N_global;
    
    % Decay factor for every pixel based on how long since its last update
    % Decay = exp(-(Current_Global - Last_Update_Map) / tau)
    % Note: We do NOT update state.S here, only the visualization.
    % The state.S must preserve the value AT its last event for accurate integration later.
    
    decay_map = exp(-(current_global_time - state.N_last) / tau_N);
    output_surface = state.S .* decay_map;
    
    % Normalize for display (Robust min-max or sigmoid)
    normalized_output_frame = normalizeSurface(output_surface);
    
end

function norm_S = normalizeSurface(S)

    % Clip values to a reasonable contrast integration range    
    % Robust Auto-scale (Ignore outliers)
    mean_val = mean(S(:));
    std_val = std(S(:));
    min_v = mean_val - 3*std_val;
    max_v = mean_val + 3*std_val;
    S_clipped = max(min(S, max_v), min_v);
    norm_S = (S_clipped - min_v) / (max_v - min_v);

end
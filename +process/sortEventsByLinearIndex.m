function [x_sorted, y_sorted, t_sorted, p_sorted, pos, unique_idx, linear_idx,...
    group_ends, fcn_time] = sortEventsByLinearIndex(x_valid, y_valid, t_valid,...
    p_valid, img_size)
    % SORTEVENTSBYLINEARINDEX Sorts events by their (x,y) coordinates.
    % 
    % [x_sorted, y_sorted, t_sorted, p_sorted, pos, unique_idx,...
    %     group_ends] = sortEventsByLinearIndex(x_valid, y_valid, t_valid,...
    %     p_valid, img_size)
    % 
    % Inputs:
    %   x_valid   - [M x 1] Unsorted x-coordinates for all events in window
    %   y_valid   - [M x 1] Unsorted t-coordinates for all events in window
    %   t_valid   - [M x 1] Unsorted timestamps for all events in window
    %   p_valid   - [M x 1] Unsorted polarity for all events in window
    %   img_size  - [1 x 2] Image dimensions [nRows, nCols]
    % 
    % Outputs:
    %   x_sorted   - [M x 1] Sorted x-coordinates for all events in window
    %   y_sorted   - [M x 1] Sorted t-coordinates for all events in window
    %   t_sorted   - [M x 1] Sorted timestamps for all events in window
    %   p_sorted   - [M x 1] Sorted polarity for all events in window
    %   pos        - [K x 1] Start index of each group in sorted_t
    %   unique_idx - [K x 1] Linear pixel index of each group
    %   group_ends - [K x 1] End index of each group in sorted_t

    % Start internal timer
    tic;

    % Convert 2D subscripts (x,y) to 1D linear indices
    linear_idx = sub2ind(img_size, x_valid, y_valid);
    
    % Sort primarily by space (linear_idx), and secondarily by time (t_valid)
    [~, sort_order] = sortrows([linear_idx, t_valid]);
    
    sorted_idx = linear_idx(sort_order);
    t_sorted = t_valid(sort_order);
    x_sorted = x_valid(sort_order);
    y_sorted = y_valid(sort_order);
    p_sorted = p_valid(sort_order);
    
    % Extract starting position of each group in the sorted list
    [unique_idx, pos, ~] = unique(sorted_idx);
    
    % Extract ending positions for each group in the sorted list
    group_ends = [pos(2:end)-1; length(sorted_idx)];

    % Get total function time
    fcn_time = toc; 
    
end
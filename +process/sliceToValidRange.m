function [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
    sliceToValidRange(t_range_n, xk, yk, tk, pk, current_idx)
    % SLICETOVALIDRANGE Extract the next time window of events and 
    % filter to valid image bounds.
    %
    %     [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
    %         sliceToValidRange(t_range_n, xk, yk, tk, pk, img_size, ...
    %         current_idx)
    %
    %     Inputs:
    %       t_range_n   - Scalar upper time bound of the window,
    %                     INCLUSIVE, in seconds.
    %       xk, yk      - [M x 1] Full event coordinate vectors
    %                     (x = row, y = col).
    %       tk          - [M x 1] Full timestamp vector in seconds,
    %                     sorted ascending.
    %       pk          - [M x 1] Full polarity vector.
    %       img_size    - [1 x 2] Image dimensions [nRows, nCols].
    %       current_idx - Scalar start index into the full vectors; acts
    %                     as the implicit lower bound of the window.
    %
    %     Outputs:
    %       current_idx - Updated index: the first event with
    %                     tk > t_range_n (length(tk)+1 if exhausted).
    %                     Pass back in on the next call.
    %       x_valid, y_valid, t_valid, p_valid
    %                   - Event vectors within the window and inside the
    %                     spatial bounds [1, nRows] x [1, nCols]. Empty
    %                     ([]) if no events fall in the window.

    % Find the window bounds
    start_node = current_idx;
    while current_idx <= length(tk) && tk(current_idx) <= t_range_n
        current_idx = current_idx + 1;
    end
    end_node = current_idx - 1;
    
    % Extract the desired slice of data
    if end_node >= start_node

        x_valid = xk(start_node:end_node);
        y_valid = yk(start_node:end_node);
        t_valid = tk(start_node:end_node);
        p_valid = pk(start_node:end_node);

    else

        x_valid = [];
        y_valid = [];
        t_valid = [];
        p_valid = [];

    end
end
function [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
    sliceToValidRange(t_range_n, xk, yk, tk, pk, imgSz, current_idx)
% SLICETOVALIDRANGE Extract the next time window of events and filter to valid image bounds.
%
%     [current_idx, x_valid, y_valid, t_valid, p_valid] = ...
%         sliceToValidRange(t_range_n, xk, yk, tk, pk, imgSz, ...
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
%       imgSz       - [1 x 2] Image dimensions [nRows, nCols].
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
%
%     Notes:
%       Stateful sequential slicer: there is no explicit lower time
%       bound; the window starts wherever the previous call left
%       current_idx. The linear scan visits each event exactly once
%       across a full sequence of calls (amortized O(M) total),
%       relying on tk being sorted ascending.
%
%       Out-of-bounds events are consumed, not deferred: current_idx
%       advances past them, so they are permanently discarded and
%       will not appear in any later window.
%
%     See also MAIN.

    % ----------------------------------------------------------------
    % 1. Linear scan to find window bounds
    % ----------------------------------------------------------------
    start_node = current_idx;
    while current_idx <= length(tk) && tk(current_idx) <= t_range_n
        current_idx = current_idx + 1;
    end
    end_node = current_idx - 1;
    
    % ----------------------------------------------------------------
    % 2. Slice and filter
    % ----------------------------------------------------------------
    if end_node >= start_node
        x_slice = xk(start_node:end_node);
        y_slice = yk(start_node:end_node);
        t_slice = tk(start_node:end_node);
        p_slice = pk(start_node:end_node);
        % Spatial bounds filter
        valid_mask = x_slice >= 1 & x_slice <= imgSz(1) & ...
            y_slice >= 1 & y_slice <= imgSz(2);
        x_valid = x_slice(valid_mask);
        y_valid = y_slice(valid_mask);
        t_valid = t_slice(valid_mask);
        p_valid = p_slice(valid_mask);
    else
        x_valid = [];
        y_valid = [];
        t_valid = [];
        p_valid = [];
    end
end
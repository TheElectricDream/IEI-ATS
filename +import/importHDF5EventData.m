function [video_out_path, xk, yk, tk, pk, t_total, buf] =...
    importHDF5EventData(file_name, filter_selection, accumulator_selection, ...
    t_start_process, t_end_process, use_buffer)

    % Set path to datasets
    hdf5_path = ['/home/alexandercrain/Dropbox/Graduate Documents' ...
        '/Doctor of Philosophy/Thesis Research/Datasets/SPOT/HDF5/'];
    
    % Set output video name
    video_out_path = fullfile('/home/alexandercrain/Videos/Research', ...
        sprintf('Normalized-Output-Fltr-%s-Acmtr-%s-%s.avi', ...
        filter_selection, accumulator_selection, char(datetime('now','Format','yyyyMMdd-HHmmss'))));
    
    if use_buffer == false
    
        % Load the data (actual event data)
        tk = double(h5read([hdf5_path file_name], '/timestamp'));
        xk = single(h5read([hdf5_path file_name], '/x'));
        yk = single(h5read([hdf5_path file_name], '/y'));
        pk = single(h5read([hdf5_path file_name], '/polarity'));
    
        % Convert time to seconds
        tk = (tk - tk(1))./1e6;
    
        % Convert to single data type to use less memory
        tk = single(tk);
    
        % Find indices within the valid range
        valid_idx = tk >= t_start_process & tk <= t_end_process;
    
        % Filter the data vectors, note the +1 because the raw data is
        % zero-indexed, and MATLAB does not handle those
        tk = tk(valid_idx);
        xk = xk(valid_idx)+1;
        yk = yk(valid_idx)+1;
        pk = pk(valid_idx);
    
        % Shift time to start at 0 for the new window
        % This ensures your frame loop starts correctly at frame 1
        tk = tk - t_start_process; 
    
        % Clear unused variables for memory
        clearvars valid_idx;
    
        % Set the time interval to accumulate over
        t_total                     = max(tk);  % [s]

        % Define an empty buffer
        buf = [];
    
    else
    
        % Create buffered event reader
        buf = import.EventBuffer(fullfile(hdf5_path, file_name), ...
            t_start_process, t_end_process);
    
        t_total                     = buf.t_total;

        % Define empty event arrays
        xk = []; yk = []; tk = []; pk = [];
    
    end

end
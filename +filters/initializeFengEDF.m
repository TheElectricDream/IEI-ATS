function [edf] = initializeFengEDF(img_size, frame_total)

    % EVENT DENSITY FILTER PARAMETERS
    % REF: https://doi.org/10.3390/app10062024

    % Tuning Guide:
    %   The Event Density Filter classifies each incoming event as signal or
    %   noise by counting the number of events that have occurred within a
    %   spatiotemporal neighbourhood. If the total count (the event density)
    %   falls below a threshold, the event is rejected as background
    %   activity. A second pass removes hot-pixel noise by checking whether
    %   the surviving events in a small 3x3 neighbourhood have any spatial
    %   support from surrounding pixels.
    %
    %   L (spatial neighbourhood size): Defines the side length, in pixels,
    %   of the square region centred on the candidate event within which
    %   supporting events are counted. Increasing L draws support from a
    %   wider area, improving robustness in sparse scenes but potentially
    %   allowing spatially distant noise events to support one another.
    %   Decreasing L restricts support to the immediate vicinity, suiting
    %   dense scenes with fine spatial structure.
    %
    %   dt (temporal neighbourhood size): Defines how far back in time, in
    %   seconds, the filter looks for supporting events. Increasing dt
    %   extends the temporal window, which is appropriate for slow-moving
    %   objects or low event rates where supporting events are spread over
    %   longer intervals. Decreasing it tightens the window, suiting fast
    %   motion where only recent events are likely to be correlated with
    %   the candidate.
    %
    %   Psi (event density threshold): The minimum number of events that
    %   must be present in the spatiotemporal neighbourhood for the
    %   candidate event to be retained. Raising Psi demands more supporting
    %   evidence before accepting an event, producing more aggressive noise
    %   removal at the risk of discarding valid signal events at the edges
    %   of moving objects. Lowering it retains more events at the risk of
    %   admitting background activity noise.
    
    edf.edf_params.L                = 5;        
    edf.edf_params.dt               = 10e-1;    
    edf.edf_params.Psi              = 3;      
    edf.edf_eventBuffer.x           = [];
    edf.edf_eventBuffer.y           = [];
    edf.edf_eventBuffer.t           = [];
    edf.edf_n_passed_store          = zeros(1, frame_total);
    edf.edf_n_total_store           = zeros(1, frame_total);

end
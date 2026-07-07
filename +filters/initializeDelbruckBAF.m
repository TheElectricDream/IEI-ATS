function [baf] = initializeDelbruckBAF(img_size, frame_total)

    % BACKGROUND ACTIVITY FILTER PARAMETERS (Delbruck 2008)
    % REF: https://www.zora.uzh.ch/handle/20.500.14742/41289
    
    % Tuning Guide:
    %   The Background Activity Filter removes uncorrelated noise events
    %   caused by thermal noise and junction leakage currents in the
    %   sensor hardware. For each incoming event, its timestamp is written
    %   into the timestamp memory of all eight immediately neighbouring
    %   pixels. The event is then accepted only if the timestamp already
    %   stored at its own pixel location indicates that a nearby event
    %   occurred within the support time T. Because genuine activity from
    %   moving objects produces spatiotemporally correlated events in
    %   adjacent pixels, these events mutually support one another and
    %   pass through. Isolated noise events, which lack recent
    %   neighbouring activity, are rejected.
    %
    %   T (support time): The maximum allowable time, in seconds, between
    %   the current event and the most recent neighbouring event for the
    %   current event to be accepted as correlated activity. Increasing T
    %   extends the support window, which is appropriate for slow-moving
    %   objects or low event rates where the time between correlated
    %   events at adjacent pixels may be long. Decreasing it tightens the
    %   window, improving rejection of background activity noise at the
    %   risk of discarding valid events from slow or sparse motion. The
    %   optimal value depends on the expected speed of objects in the
    %   scene relative to the sensor resolution: faster motion produces
    %   shorter inter-event intervals at neighbouring pixels, permitting
    %   a smaller T.
    
    baf.baf_params.T                = 10e-1;   
    baf.baf_lastTimesMap            = -Inf(img_size);
    baf.baf_n_passed_store          = zeros(1, frame_total);
    baf.baf_n_total_store           = zeros(1, frame_total);

end
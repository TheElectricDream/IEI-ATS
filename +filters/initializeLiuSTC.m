function [stc] = initializeLiuSTC(img_size, frame_total)

    % SPATIOTEMPORAL CORRELATION FILTER PARAMETERS 
    % REF: http://ieeexplore.ieee.org/document/7168735/

    % Tuning Guide:
    %   The Spatiotemporal Correlation Filter removes background activity
    %   (BA) noise by exploiting the fact that real events arising from
    %   moving features are temporally correlated with nearby events,
    %   whereas BA noise events occur randomly and in isolation. The
    %   sensor pixel array is spatially subsampled so that blocks of
    %   pixels map onto a single filter cell. When an event arrives at a
    %   cell, it opens a time window during which any subsequent event
    %   mapping to the same cell is accepted as correlated activity. If
    %   no prior event has occurred in the cell within the time window,
    %   the incoming event is rejected as uncorrelated noise.
    %
    %   dT (temporal correlation window): The duration, in seconds, for
    %   which a filter cell remains "open" after receiving an event.
    %   Any subsequent event that maps to the same cell within this
    %   window is classified as correlated and allowed to pass.
    %   Increasing dT extends the window, which is appropriate for
    %   slow-moving objects or low event rates where the inter-event
    %   interval within a cell may be long. Decreasing it tightens the
    %   window, improving noise rejection at the risk of discarding
    %   valid events from slow or sparse activity. The original hardware
    %   design supports a range of 5 ns to 1 s, with typical values of
    %   1-10 ms.
    %
    %   subsample_rate (spatial neighbourhood exponent): Controls the
    %   size of the pixel block that maps onto a single filter cell.
    %   A value of S causes each block of 2^S x 2^S sensor pixels to
    %   share one filter cell (e.g. S=0 gives 1x1, S=1 gives 2x2, S=2
    %   gives 4x4). Increasing this value widens the spatial
    %   neighbourhood, so that events from a larger region of the sensor
    %   provide temporal support for one another. This improves
    %   robustness for sparse or slow-moving scenes but reduces spatial
    %   selectivity, potentially allowing spatially distant noise events
    %   to support each other. Decreasing it restricts support to fewer
    %   pixels, preserving spatial precision at the cost of requiring
    %   denser local activity for events to pass.
    
    stc.stc_params.dT               = 10e-3;    
    stc.stc_params.subsample_rate   = 2;       
    stc.stc_block_size              = 2^stc.stc_params.subsample_rate;
    stc.stc_cellSz                  = ceil(img_size ./ stc.stc_block_size);
    stc.stc_lastTimesMap            = -Inf(stc.stc_cellSz);
    stc.stc_n_passed_store          = zeros(1, frame_total);
    stc.stc_n_total_store           = zeros(1, frame_total);

end
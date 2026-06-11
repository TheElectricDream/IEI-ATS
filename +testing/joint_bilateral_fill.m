function Map_Filled = joint_bilateral_fill(Map, SpatialSigma, RangeSigma)
    % Map: 640x480 matrix with NaNs
    if nargin < 2, SpatialSigma = 7; end
    if nargin < 3, RangeSigma = 20; end % Adjusted for the raw [0, 180] data range
    
    Hole_Mask = Map==0;
    
    % Step 1: Create a structurally sound guidance map using regionfill.
    % This cleanly bridges the hole geometrically without blurring the edge.
    Guidance = regionfill(Map, Hole_Mask); 
    
    % Step 2: Use imbilatfilt (highly compatible across older toolbox versions).
    % We pass DegreeOfSmoothing (derived from RangeSigma) and SpatialSigma.
    % To ensure sharp boundaries are preserved across the 0-180 step, 
    % the DegreeOfSmoothing should match the expected variance scale.
    Map_Filtered = imbilatfilt(Guidance, RangeSigma^2, SpatialSigma);
                    
    % Step 3: Keep the original valid data pristine, only use the 
    % edge-preserving filtered data to patch the original holes.
    Map_Filled = Map;
    Map_Filled(Hole_Mask) = Map_Filtered(Hole_Mask);
end
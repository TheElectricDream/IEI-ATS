function Map_Filled = fast_bilateral_inpainting(Map, SpatialSigma, RangeSigma)
% Map: 640x480 matrix with NaNs
if nargin < 2, SpatialSigma = 5; end
if nargin < 3, RangeSigma = 25; end % Threshold for the 180-unit step

% 1. Create binary mask of valid data (1 = valid, 0 = hole)
Mask = double(Map~=0);

% 2. Initialize NaNs to 0 so matrix multiplication works cleanly
Map_Zeroed = Map;
Map_Zeroed(Map==0) = 0;

% 3. Create a fast spatial Gaussian kernel
KernelSize = ceil(SpatialSigma * 3) * 2 + 1;
K_spatial = fspecial('gaussian', [KernelSize, KernelSize], SpatialSigma);

% 4. Get an initial fast structural estimate to compute range distances
% (Using a fast box filter over the mask to get local consensus)
Local_Mean = imfilter(Map_Zeroed, K_spatial) ./ (imfilter(Mask, K_spatial) + eps);

% 5. Compute the Range Weight Matrix dynamically
% This penalizes pixels that deviate too much from the local structure
Range_Weights = exp(-(Map_Zeroed - Local_Mean).^2 / (2 * RangeSigma^2));

% 6. Apply the combined spatial and range weights simultaneously
% Effectively: Filled = Convolve(Data * Mask * Range) / Convolve(Mask * Range)
Weighted_Data = Map_Zeroed .* Mask .* Range_Weights;
Weighted_Mask = Mask .* Range_Weights;

Data_Buffered = imfilter(Weighted_Data, K_spatial, 'replicate');
Mask_Buffered = imfilter(Weighted_Mask, K_spatial, 'replicate');

% Divide by the weight map to normalize (plus eps to avoid division by zero)
Map_Filled = Data_Buffered ./ (Mask_Buffered + eps);

% 7. Preserve pristine original valid pixels
Map_Filled(Mask == 1) = Map(Mask == 1);
end
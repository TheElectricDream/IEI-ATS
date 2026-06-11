function Map_Filled = anisotropic_diffusion_fill(Map, Iterations, Kappa)
    % Map: 640x480 matrix with holes represented as NaNs
    % Iterations: Number of diffusion steps (e.g., 50-100)
    % Kappa: Conduction threshold (controls sharpness retention)
    
    if nargin < 2, Iterations = 80; end
    if nargin < 3, Kappa = 15; end 
    
    Map_Filled = Map;
    Hole_Mask = Map==0;
    
    % FIX: Use regionfill to smoothly interpolate the initial guess 
    % across the 2D matrix boundaries so gradients don't break on step 1.
    Map_Filled = regionfill(Map, Hole_Mask);
    
    % Integration constant (dt <= 1/4 for stability in 2D)
    delta_t = 0.25; 
    
    for i = 1:Iterations
        % Compute forward spatial differences
        [Gx, Gy] = gradient(Map_Filled);
        Grad_Mag = sqrt(Gx.^2 + Gy.^2);
        
        % Perona-Malik conduction coefficient 
        c = exp(-(Grad_Mag ./ Kappa).^2);
        
        % Compute divergence via discrete Laplacian variations
        [N, S, E, W] = get_neighbor_diffs(Map_Filled);
        
        % Update map *only* inside the original hole locations
        Flux = (c.*N + c.*S + c.*E + c.*W);
        Map_Filled(Hole_Mask) = Map_Filled(Hole_Mask) + delta_t * Flux(Hole_Mask);
    end
end

function [N, S, E, W] = get_neighbor_diffs(I)
    % Shift differences for numerical divergence
    N = imfilter(I, [0 1 0; 0 -1 0; 0 0 0], 'replicate');
    S = imfilter(I, [0 0 0; 0 -1 0; 0 1 0], 'replicate');
    E = imfilter(I, [0 0 0; 0 -1 1; 0 0 0], 'replicate');
    W = imfilter(I, [0 0 0; 1 -1 0; 0 0 0], 'replicate');
end
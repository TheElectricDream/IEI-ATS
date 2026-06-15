function [grad_mag, nx, ny, nz] = calculateLocalFrameGradient(frame)

    % Calculate the gradient information
    [dz_dy, dz_dx] = gradient(frame);  
    
    % Compute gradient magnitude and normal vectors for the surface z = Z(x,y).
    grad_mag = sqrt(dz_dx.^2 + dz_dy.^2);
    
    % Surface normals (unnormalized): n = [-dz_dx; -dz_dy; 1]
    nx = -dz_dx;
    ny = -dz_dy;
    nz = ones(size(filtered_coherence_map));
    n_norm = sqrt(nx.^2 + ny.^2 + nz.^2);
    nx = nx ./ n_norm;
    ny = ny ./ n_norm;
    nz = nz ./ n_norm;

end
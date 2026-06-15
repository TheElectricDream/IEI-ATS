function [similarity_score, cv_map, regularity_map] = findRegularity(sorted_x, sorted_y, std_map, mean_map, imgSz)

    % --- 1. observation mask ----------------------------------------
    obs_mask = mean_map > 0;

    % --- 2. per-pixel CV --------------------------------------------
    cv_map = std_map ./ max(mean_map, eps);
    cv_map(~obs_mask) = Inf;            % unobserved -> score 0 downstream

    % --- 3. CV-based regularity -------------------------------------
    regularity_map = 1 ./ (1 + cv_map);        % -> 0 where unobserved

    % --- 4. per-event lookup ----------------------------------------
    linear_idx = sub2ind(imgSz, sorted_x(:), sorted_y(:));
    similarity_score = regularity_map(linear_idx);

end

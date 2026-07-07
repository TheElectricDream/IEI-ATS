function [map, map_raw, outlier_mask] = rejectPolarityOutliers(map, map_raw, n_sigma)
% REJECTPOLARITYOUTLIERS Robust, one-sided, per-polarity outlier replacement using MAD.
%
%     [map, map_raw, outlier_mask] = rejectPolarityOutliers(map, map_raw, n_sigma)
%
%     Outputs:
%       map          - Input with outliers replaced by the median of its population.
%       map_raw      - Companion array updated identically.
%       outlier_mask - Logical array of the same size, true where outliers were found.

    % Initialize the universal tracking mask
    outlier_mask = false(size(map));
    
    for pol = [1, -1]
    
        if pol == 1
            mask = map > 0;
        else
            mask = map < 0;
        end
    
        vals = map(mask);
    
        if isempty(vals), continue; end
    
        med = median(vals);
        mad_val = median(abs(vals - med));
        
        robust_sig = 1.4826 * mad_val;
        if robust_sig == 0
            robust_sig = std(vals);
        end
    
        if pol == 1
            outliers = map > (med + n_sigma * robust_sig);
        else
            outliers = map < (med - n_sigma * robust_sig);
        end
        
        % Record these specific polarity outliers into the universal mask
        outlier_mask(outliers) = true;
        
        % Replace outliers with the robust median
        map(outliers)     = med;
        map_raw(outliers) = median(map_raw(mask));
    end
end
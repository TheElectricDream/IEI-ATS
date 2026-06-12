function [map, map_raw] = rejectPolarityOutliers(map, map_raw, n_sigma)
% REJECTPOLARITYOUTLIERS One-sided, per-polarity outlier replacement with the polarity median.
%
%     [map, map_raw] = rejectPolarityOutliers(map, map_raw, n_sigma)
%
%     Inputs:
%       map     - Signed array; positive and negative values are
%                 treated as separate polarity populations (zeros
%                 belong to neither).
%       map_raw - Companion array of the same size, clipped with the
%                 outlier mask derived from map so that display and
%                 feedback paths stay consistent.
%       n_sigma - Positive scalar; number of standard deviations
%                 defining the outlier cut.
%
%     Outputs:
%       map     - Input with outliers replaced by the median of its
%                 own polarity population.
%       map_raw - Input with the same outlier pixels replaced by the
%                 median of map_raw over that polarity's mask.
%
%     Notes:
%       The cut is ONE-SIDED, away from zero: positive values are
%       flagged only if above mu + n_sigma*sigma of the positive
%       population; negative values only if below mu - n_sigma*sigma
%       of the negative population. Values near zero are never
%       flagged.
%
%       Statistics (mean, std, median) are computed in a single pass
%       over each polarity population including the outliers
%       themselves; this is not iterative sigma clipping.
%
%     See also MEDIAN, STD, FILLOUTLIERS.
    
    for pol = [1, -1]
    
        if pol == 1
            mask = map > 0;
        else
            mask = map < 0;
        end
    
        vals = map(mask);
    
        if isempty(vals), continue; end
    
        mu  = mean(vals);
        sig = std(vals);
        med = median(vals);
    
        if pol == 1
            outliers = map > (mu + n_sigma * sig);
        else
            outliers = map < (mu - n_sigma * sig);
        end

        map(outliers)     = med;
        map_raw(outliers) = median(map_raw(mask));
    end
end
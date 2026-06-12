function [removed_pixels] = recoverRemovedEvents(filter_mask, counts)
% RECOVERREMOVEDEVENTS Binary map of active pixels rejected by a filter mask.
%
%     removed_pixels = recoverRemovedEvents(filter_mask, counts)
%
%     Inputs:
%       filter_mask - [H x W] Filter survival mask; nonzero means the
%                     pixel was kept, zero means rejected. Assumed
%                     nonnegative.
%       counts      - [H x W] Per-pixel event counts for the window.
%
%     Outputs:
%       removed_pixels - [H x W] Double map in {0, 1}; 1 where a
%                        pixel received events (counts > 0) but was
%                        rejected by the mask (filter_mask == 0).
%
%     Notes:
%       Equivalent to double(filter_mask == 0 & counts > 0) for a
%       nonnegative mask. Negative filter_mask entries are silently
%       treated as kept (they survive the inversion with a negative
%       value, zeroing the comparison).
%
%     See also COHERENCE.COMPUTECOHERENCEMASK.

    temp = filter_mask;
    temp(filter_mask>0)=0;
    temp(filter_mask==0)=1;
    counts_temp = counts.*temp;
    removed_pixels = (counts_temp>0).*1.0;
end
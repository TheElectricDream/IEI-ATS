function printPercentComplete(index, total, duration_filter, duration_accumulator)
% PRINTPERCENTCOMPLETE Print processing progress and per-stage timing to the console.
%
%     printPercentComplete(index, total, duration_filter, ...
%         duration_accumulator)
%
%     Inputs:
%       index                - Scalar current iteration index.
%       total                - Scalar total number of iterations.
%       duration_filter      - Elapsed time of the current filtering
%                              stage, in seconds.
%       duration_accumulator - Elapsed time of the current
%                              accumulation stage, in seconds.
%
%     Notes:
%       Prints a three-line block: percentage complete (one decimal
%       place) followed by the two stage timings (millisecond
%       precision). No outputs.
%
%     See also FPRINTF.

    percent_complete = index / total * 100;
    fprintf('\nProcessing... %.1f%% complete. Stats: \n \t Filtering: (%.3f seconds)\n  \t Accumulation: (%.3f seconds)\n', ...
        percent_complete, duration_filter, duration_accumulator);
end
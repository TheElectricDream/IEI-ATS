function [filter_mask, stc, baf, edf, stcc, mcf] = ...
    applyChosenEventFilter(filter_selection, stc, baf, edf, stcc, mcf,...
    x_sorted, y_sorted, t_sorted, p_sorted, counts)
% APPLYCHOSENFILTER runs the currently selected filter for comparison with
% the filtering approach.

    if strcmp(filter_selection, 'STC') == 1

        % -------------------- STC CORRELATION FILTER --------------------%
        % ----------------------------------------------------------------%

        [filter_mask, stc.stc_lastTimesMap, stc_n_pass, stc_n_tot] = ...
            filters.spatiotemporalCorrelation(x_sorted, y_sorted, ...
            t_sorted, img_size, stc.stc_lastTimesMap, stc.stc_params);

        % Store metrics for post-processing analysis
        stc.stc_n_passed_store(frame_index) = stc_n_pass;
        stc.stc_n_total_store(frame_index)  = stc_n_tot;

    elseif strcmp(filter_selection, 'BAF') == 1

        % ------------- BACKGROUND ACTIVITY FILTER (BAF) -----------------%
        % ----------------------------------------------------------------%

        [filter_mask, baf.baf_lastTimesMap, baf_n_pass, baf_n_tot] = ...
            filters.runDelbruckBAF(x_sorted, y_sorted, ...
            t_sorted, img_size, baf.baf_lastTimesMap, baf.baf_params);

        % Store metrics for post-processing analysis
        baf.baf_n_passed_store(frame_index) = baf_n_pass;
        baf.baf_n_total_store(frame_index)  = baf_n_tot;

    elseif strcmp(filter_selection, 'EDF') == 1

        % ------------ EVENT DENSITY FILTER (EDF) ------------------------%
        % ----------------------------------------------------------------%

        [filter_mask, edf.edf_eventBuffer, edf_n_pass, edf_n_tot] = ...
            filters.eventDensityFilter(x_sorted, y_sorted, ...
            t_sorted, img_size, edf.edf_eventBuffer, edf.edf_params);

        % Store metrics for post-processing analysis
        edf.edf_n_passed_store(frame_index) = edf_n_pass;
        edf.edf_n_total_store(frame_index)  = edf_n_tot;

    elseif strcmp(filter_selection, 'STCC') == 1

        % ------- STCC-FILTER (Space-Time-Content Correlation) -----------%
        % ----------------------------------------------------------------%

        % Compute signed polarity for this filter
        p_signed_for_stcc = p_sorted * 2 - 1;

        [filter_mask, stcc_n_pass, stcc_n_tot, stcc_threshold] = ...
            filters.stccFilter(x_sorted, y_sorted, ...
            t_sorted, p_signed_for_stcc, img_size, stcc.stcc_params);

        % Store metrics for post-processing analysis
        stcc.stcc_n_passed_store(frame_index) = stcc_n_pass;
        stcc.stcc_n_total_store(frame_index)  = stcc_n_tot;
        stcc.stcc_threshold(frame_index) = stcc_threshold;

    elseif strcmp(filter_selection, 'MCF') == 1
 
        % ----------- MOTION CONSISTENCY FILTER (MCF) ----------------%
        % ------------------------------------------------------------%

        [filter_mask, mcf.mcf_eventBuffer, mcf_n_pass, mcf_n_tot] = ...
            filters.motionConsistencyFilter(x_sorted, y_sorted, ...
            t_sorted, img_size, mcf.mcf_eventBuffer, mcf.mcf_params);
 
        % Store metrics for post-processing analysis
        mcf.mcf_n_passed_store(frame_index) = mcf_n_pass;
        mcf.mcf_n_total_store(frame_index)  = mcf_n_tot;
        
    elseif strcmp(filter_selection, 'NONE') == 1

        % Use a unity mask instead
        filter_mask = (counts>0).*1.0;

    end

end
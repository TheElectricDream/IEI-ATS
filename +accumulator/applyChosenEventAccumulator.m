function [agd, hots, sits, mets, zhu] = ...
    applyChosenEventAccumulator(accumulator_selection, agd, hots, sits, mets, zhu)

    if strcmp(accumulator_selection, 'AGD') == 1

        % ------------------- NUNES GLOBAL ADAPTIVE (AGD) ----------------%
        % ----------------------------------------------------------------%
    
        % Run the AGD algorithm
        [agd.agd_surface, agd.agd_state, ~] = accumulator.adaptiveGlobalDecay(agd.agd_surface,...
            filtered_x, filtered_y, filtered_t, agd.agd_state, agd.agd_params);

        agd.normalized_output_frame = agd.agd_surface;
        agd.agd_activity_store(frame_index) = agd.agd_state.activity;

    elseif strcmp(accumulator_selection, 'HOTS') == 1        
    
        % ------------- TIME-SURFACE ACCUMULATION (HOTS) -----------------%
        % ----------------------------------------------------------------%
        
        [hots.ts_t_map, hots.normalized_output_frame] = accumulator.timeSurface(hots.ts_t_map,...
         filtered_x, filtered_y, filtered_t, img_size, hots.ts_time_constant);

    elseif strcmp(accumulator_selection, 'SITS') == 1

        % ------------ SPEED INVARIENT TIME-SURFACE (SITS) ---------------%
        % ----------------------------------------------------------------%
        
        [sits.sits_t_map, sits.normalized_output_frame] = ...
            accumulator.speedInvariantTimeSurface(sits.sits_t_map, filtered_x,...
            filtered_y, sits.sits_R);

    elseif strcmp(accumulator_selection, 'METS') == 1

        % ------------ MOTION-ENCODED TIME-SURFACE (METS) ----------------%
        % ----------------------------------------------------------------%

        [mets.mets_surface, mets.mets_state, mets.normalized_output_frame] = ...
            accumulator.motionEncodedTimeSurface(...
            filtered_x, filtered_y, filtered_t, p_signed, ...
            img_size, mets.mets_state, mets.mets_params);

    elseif strcmp(accumulator_selection, 'EVO-ATS') == 1

        % ---------- ZHU ADAPTIVE TIME SURFACE (EVO-ATS) -----------------%
        % ----------------------------------------------------------------%
    
        [zhu.zhu_surface, zhu.zhu_state, zhu.normalized_output_frame, zhu.zhu_tau_map] = ...
            accumulator.adaptiveTimeSurfaceZhu(...
            filtered_x, filtered_y, filtered_t, p_signed, ...
            img_size, zhu.zhu_state, zhu.zhu_params);

    else

        error('An invalid accumulator has been selected!')

    end

end
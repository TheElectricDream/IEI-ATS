function [stcc] = initializeLiSTCC(img_size, frame_total)

    % STCC-FILTER PARAMETERS
    % REF: https://doi.org/10.1016/j.image.2024.117136

    % Tuning Guide:
    %   The Space–Time-Content Correlation computes a probability-of-signal for
    %   each incoming event by evaluating its spatiotemporal correlation
    %   with N neighbouring events. Events whose POS falls below the
    %   discrimination threshold TH are classified as noise and rejected.
    %
    %   N (context window size): The number of recent events used to
    %   evaluate each candidate. Larger values provide more statistical
    %   support and improve robustness, but increase computation time.
    %   Smaller values are faster but may produce unreliable POS estimates.
    %
    %   sigma_d (spatial scale): Controls how rapidly the influence of
    %   neighbouring events decays with spatial distance, in pixels.
    %   Increasing this value extends the spatial support region, which is
    %   appropriate for sparse scenes or broad edge structures. Decreasing
    %   it restricts support to the immediate vicinity, suiting dense or
    %   fine-featured scenes.
    %
    %   sigma_t (temporal scale): Controls how rapidly the influence of
    %   neighbouring events decays with temporal distance, in seconds.
    %   Increasing this value extends the temporal support window, which is
    %   appropriate for slow-moving objects or low event rates. Decreasing
    %   it tightens the window, suiting fast motion where temporally
    %   distant events are unlikely to be correlated.
    %
    %   sigma_p (polarity scale): Controls how strongly polarity agreement
    %   is weighted. Decreasing this value enforces stricter same-polarity
    %   agreement between neighbouring events. Increasing it relaxes this
    %   constraint, permitting greater mixed-polarity support.
    %
    %   TH (base discrimination threshold): The minimum POS value an event
    %   must exceed to be retained. Raising TH produces more aggressive
    %   filtering at the risk of discarding valid signal events. Lowering
    %   it retains more events at the risk of admitting noise. Note that
    %   the content-correlation mechanism adjusts TH at run-time, so this
    %   value serves as the initial operating point for that adaptation.
    
    stcc.stcc_params.N               = 100;      
    stcc.stcc_params.sigma_d         = 6;       
    stcc.stcc_params.sigma_t         = 10e-1;   
    stcc.stcc_params.sigma_p         = 5;        
    stcc.stcc_params.TH              = 0.1;     
    stcc.stcc_n_passed_store         = zeros(1, frame_total);
    stcc.stcc_n_total_store          = zeros(1, frame_total);
    stcc.stcc_threshold              = zeros(1, frame_total);

end
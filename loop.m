%% SETUP
% Clear out all variables
clear;    
clc;     
close all; 

% Define list of filter names to iterate over
filtersNames = {'NONE', 'STC', 'BAF', 'EDF', 'STCC', 'MCF', 'COH'};

% Define the list of accumulators to iterate over
accumulatorNames = {'HOTS', 'SITS', 'METS', 'IEI-ATS', 'AGD', 'EVO-ATS'};

% Define the list of experiments to loop over
%   1 - NOM-ROT
%   2 - SG-ROT
%   3 - DARK-ROT
% experimentsNames = {'recording_20251029_131131.hdf5',...
%                     'recording_20251029_135047.hdf5',...
%                     'recording_20251029_134602.hdf5'}; 

experimentsNames = {'recording_20251029_131131.hdf5'}; 

%% LOOPING
% Preallocate cell array to store metrics for each filter run
loopLogFilters = cell(1, numel(filtersNames));
loopLogAccumulators = cell(1, numel(accumulatorNames));
loopLogExperiments = cell(1, numel(experimentsNames));

% Loop through all filters and accumulators for all experiments
for expIdx = 1:numel(experimentsNames)

    fileName = experimentsNames{expIdx};

    % % Loop through each filter name, set selection and run main.m
    % for loopIdx = 1:numel(filtersNames)
    %     filterSelection = filtersNames{loopIdx}; 
    %     accumulatorSelection = 'IEI-ATS';
    %     run("main.m");                           
    %     loopLogFilters{loopIdx} = frame_metrics;       
    % end
    
    % Loop through each filter name, set selection and run main.m
    for loopIdx = 1:numel(accumulatorNames)
        filterSelection = 'COH';
        accumulatorSelection = accumulatorNames{loopIdx}; 
        run("main.m");                           
        loopLogAccumulators{loopIdx} = frame_metrics;   
    end

    loopLogExperiments{expIdx} = {loopLogFilters, loopLogAccumulators};
end
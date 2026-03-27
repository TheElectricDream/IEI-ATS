%% SETUP
% Clear out all variables
clear;    
clc;     
close all; 

% Define list of filter names to iterate over
filtersNames = {'NONE', 'STC', 'BAF', 'EDF', 'STCC', 'MCF', 'COH'};

% Define the list of accumulators to iterate over
accumulatorNames = {'HOTS', 'SITS', 'METS', 'IEI-ATS', 'AGD', 'EVO-ATS'};

%% LOOPING
% Preallocate cell array to store metrics for each filter run
loopLogFilters = cell(1, numel(filtersNames));
loopLogAccumulators = cell(1, numel(accumulatorNames));

% Loop through each filter name, set selection and run main.m
for loopIdx = 1:numel(filtersNames)
    filterSelection = filtersNames{loopIdx}; % Current filter to use in main.m
    accumulatorSelection = 'IEI-ATS';
    run("main.m");                           % Execute main script which should use filterSelection
    loopLogFilters{loopIdx} = frame_metrics;        % Save resulting metrics from main.m into log
end

% Loop through each filter name, set selection and run main.m
for loopIdx = 1:numel(accumulatorNames)
    filterSelection = 'EDF';
    accumulatorSelection = accumulatorNames{loopIdx}; % Current filter to use in main.m
    run("main.m");                           % Execute main script which should use filterSelection
    loopLogAccumulators{loopIdx} = frame_metrics;        % Save resulting metrics from main.m into log
end
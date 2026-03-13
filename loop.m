clear;    
clc;     
close all; 

% Define list of filter names to iterate over
filtersNames = {'NONE', 'STC', 'BAF', 'EDF', 'STCC', 'MCF', 'COH'};

% Preallocate cell array to store metrics for each filter run
loopLog = cell(1, numel(filtersNames));

% Loop through each filter name, set selection and run main.m
for loopIdx = 1:numel(filtersNames)
    filterSelection = filtersNames{loopIdx}; % Current filter to use in main.m
    run("main.m");                           % Execute main script which should use filterSelection
    loopLog{loopIdx} = frame_metrics;        % Save resulting metrics from main.m into log
end
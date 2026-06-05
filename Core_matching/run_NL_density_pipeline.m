close all; clc; clear;

projectDir = 'C:\Users\evsalg001\Documents\MATLAB\Density parametrization\Nansen Legacy\Core_matching';
rawDir = 'C:\Users\evsalg001\Downloads\Nansen Legacy raw data';

codeDir = fullfile(projectDir, 'code');
intermediateDir = fullfile(projectDir, 'data_intermediate');
exportDir = fullfile(projectDir, 'export');

if ~exist(intermediateDir, 'dir')
    mkdir(intermediateDir);
end

if ~exist(exportDir, 'dir')
    mkdir(exportDir);
end

addpath(codeDir)

matchFile = NL_1_core_matching(rawDir, intermediateDir);

valuesFile = NL_2_import_of_matched_cores(matchFile, rawDir, intermediateDir);

[sectionsFile, coreAverageFile] = NL_3_core_processing(valuesFile, exportDir);

fprintf('\nPipeline complete.\n');
fprintf('Sections: %s\n', sectionsFile);
fprintf('Core average: %s\n', coreAverageFile);
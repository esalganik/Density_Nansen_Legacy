function matchFile = NL_1_core_matching(rootDir, outputDir)

densityFiles = findFiles(rootDir, '*density*.csv', "density");
temperatureFiles = findFiles(rootDir, '*temperature*.csv', "temperature");
salinityFiles = findFiles(rootDir, '*salinity*.csv', "salinity");

densityTable = filesToTable(densityFiles, "density");
temperatureTable = filesToTable(temperatureFiles, "temperature");
salinityTable = filesToTable(salinityFiles, "salinity");

tripleMatchTable = densityTable(:, {'MatchKey','DensityFile','DensityFolder','DensityPath'});

tripleMatchTable.TemperatureFile = strings(height(tripleMatchTable), 1);
tripleMatchTable.TemperatureFolder = strings(height(tripleMatchTable), 1);
tripleMatchTable.TemperaturePath = strings(height(tripleMatchTable), 1);

tripleMatchTable.SalinityFile = strings(height(tripleMatchTable), 1);
tripleMatchTable.SalinityFolder = strings(height(tripleMatchTable), 1);
tripleMatchTable.SalinityPath = strings(height(tripleMatchTable), 1);

for i = 1:height(tripleMatchTable)

    key = tripleMatchTable.MatchKey(i);

    tempIdx = find(temperatureTable.MatchKey == key, 1, 'first');
    salIdx = find(salinityTable.MatchKey == key, 1, 'first');

    if ~isempty(tempIdx)
        tripleMatchTable.TemperatureFile(i) = temperatureTable.TemperatureFile(tempIdx);
        tripleMatchTable.TemperatureFolder(i) = temperatureTable.TemperatureFolder(tempIdx);
        tripleMatchTable.TemperaturePath(i) = temperatureTable.TemperaturePath(tempIdx);
    end

    if ~isempty(salIdx)
        tripleMatchTable.SalinityFile(i) = salinityTable.SalinityFile(salIdx);
        tripleMatchTable.SalinityFolder(i) = salinityTable.SalinityFolder(salIdx);
        tripleMatchTable.SalinityPath(i) = salinityTable.SalinityPath(salIdx);
    end
end

tripleMatchTable.HasTemperature = tripleMatchTable.TemperatureFile ~= "";
tripleMatchTable.HasSalinity = tripleMatchTable.SalinityFile ~= "";
tripleMatchTable.HasTriple = tripleMatchTable.HasTemperature & tripleMatchTable.HasSalinity;

densityWithTriple = tripleMatchTable(tripleMatchTable.HasTriple, :);
densityWithoutTriple = tripleMatchTable(~tripleMatchTable.HasTriple, :);

% fprintf('\n========================================\n');
% fprintf('Density - Temperature - Salinity filename matching\n');
% fprintf('========================================\n');
% fprintf('Density files: %d\n', height(densityTable));
% fprintf('Temperature files: %d\n', height(temperatureTable));
% fprintf('Salinity files: %d\n', height(salinityTable));
% fprintf('Density files with full triple: %d\n', height(densityWithTriple));
% fprintf('Density files without full triple: %d\n\n', height(densityWithoutTriple));
% 
% disp(tripleMatchTable)
% 
% fprintf('\nDensity files WITH full triple:\n');
% disp(densityWithTriple(:, {'MatchKey','DensityFile','TemperatureFile','SalinityFile'}))
% 
% fprintf('\nDensity files WITHOUT full triple:\n');
% disp(densityWithoutTriple(:, {'MatchKey','DensityFile','TemperatureFile','SalinityFile'}))

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

matchFile = fullfile(outputDir, ...
    'density_temperature_salinity_file_matches.csv');

writetable(tripleMatchTable, matchFile);

fprintf('Saved matching table to:\n%s\n', matchFile);

%% Helpers
function files = findFiles(rootDir, pattern, variableName)

files = dir(fullfile(rootDir, '**', pattern));
files = files(~[files.isdir]);

fileNames = {files.name}';

isSummary = contains(fileNames, 'summary', 'IgnoreCase', true);
isComputedQC = contains(fileNames, 'computed_qc', 'IgnoreCase', true);

if variableName == "temperature"
    isVariableQC = contains(fileNames, 'temperature_qc', 'IgnoreCase', true);
else
    isVariableQC = false(size(fileNames));
end

isGenerated = contains(fileNames, 'density_temperature_salinity_file_matches', ...
    'IgnoreCase', true) | ...
    contains(fileNames, 'processed_density', ...
    'IgnoreCase', true);

files = files(~isSummary & ~isComputedQC & ~isVariableQC & ~isGenerated);

end

function T = filesToTable(files, variableName)

fileNames = string({files.name}');
folders = string({files.folder}');
paths = fullfile(folders, fileNames);
matchKeys = makeMatchKey(fileNames);

switch variableName
    case "density"
        T = table(matchKeys, fileNames, folders, paths, ...
            'VariableNames', {'MatchKey','DensityFile','DensityFolder','DensityPath'});
    case "temperature"
        T = table(matchKeys, fileNames, folders, paths, ...
            'VariableNames', {'MatchKey','TemperatureFile','TemperatureFolder','TemperaturePath'});
    case "salinity"
        T = table(matchKeys, fileNames, folders, paths, ...
            'VariableNames', {'MatchKey','SalinityFile','SalinityFolder','SalinityPath'});
end

end

function key = makeMatchKey(fileNames)

key = string(fileNames);

key = regexprep(key, '\.csv$', '', 'ignorecase');

key = regexprep(key, '(?i)_density$', '');
key = regexprep(key, '(?i)_temperature$', '');
key = regexprep(key, '(?i)_salinity$', '');

key = regexprep(key, '(?i)_?computed_qc', '');

key = regexprep(key, '(?i)_main$', '');
key = regexprep(key, '(?i)_optics$', '');
key = regexprep(key, '(?i)_ROV$', '');
key = regexprep(key, '(?i)_site2$', '');
key = regexprep(key, '(?i)_ramses$', '');
key = regexprep(key, '(?i)_thinice$', '');
key = regexprep(key, '(?i)_gem2calcsite$', '');

key = regexprep(key, '\s*_\s*', '_');
key = strtrim(key);

end

end
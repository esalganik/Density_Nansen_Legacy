function valuesOutFile = NL_2_import_of_matched_cores(matchFile, dataRootDir, outputDir)

matchTable = readtable(matchFile, 'TextType', 'string');
matchTable = matchTable(matchTable.HasTriple == true, :);

matchTable.DensityPath = findFilePath(dataRootDir, matchTable.DensityFile);
matchTable.SalinityPath = findFilePath(dataRootDir, matchTable.SalinityFile);
matchTable.TemperaturePath = findFilePath(dataRootDir, matchTable.TemperatureFile);

allCoreValuesTable = table();
coreAverageTable = table();

for i = 1:height(matchTable)

    densityData = readDensityCore(matchTable.DensityPath(i));
    salinityData = readSalinityCore(matchTable.SalinityPath(i));
    temperatureData = readTemperatureCore(matchTable.TemperaturePath(i));

    if isempty(densityData) || isempty(salinityData) || isempty(temperatureData)
        fprintf('Skipped incomplete triple: %s\n', matchTable.DensityFile(i));
        continue
    end

    densityDepth = mean([densityData.MinDepth densityData.MaxDepth], 2);
    salinityDepth = mean([salinityData.MinDepth salinityData.MaxDepth], 2);
    temperatureDepth = temperatureData.Depth;

    maxDensityDepth = max(densityDepth, [], 'omitnan');
    maxSalinityDepth = max(salinityDepth, [], 'omitnan');
    maxTemperatureDepth = max(temperatureDepth, [], 'omitnan');

    densityIceThickness = densityData.IceThickness(1);
    salinityIceThickness = salinityData.IceThickness(1);
    temperatureIceThickness = temperatureData.IceThickness(1);
    densityFreeboard = densityData.Freeboard(1);

    if isnan(densityIceThickness)
        densityIceThickness = maxDensityDepth;
    end

    if isnan(salinityIceThickness)
        salinityIceThickness = maxSalinityDepth;
    end

    if isnan(temperatureIceThickness)
        temperatureIceThickness = maxTemperatureDepth;
    end

    densityNormDepth = densityDepth ./ maxDensityDepth;
    salinityNormDepth = salinityDepth ./ maxSalinityDepth;
    temperatureNormDepth = temperatureDepth ./ maxTemperatureDepth;

    salinityInterp = interp1( ...
        salinityNormDepth, ...
        salinityData.Salinity, ...
        densityNormDepth, ...
        'linear', ...
        'extrap');

    temperatureInterp = interp1( ...
        temperatureNormDepth, ...
        temperatureData.Temperature, ...
        densityNormDepth, ...
        'linear', ...
        'extrap');

    oneCoreTable = densityData;
    oneCoreTable.IceThickness(:) = densityIceThickness;
    oneCoreTable.DensityDepth = densityDepth;
    oneCoreTable.InterpolatedSalinity = salinityInterp;
    oneCoreTable.InterpolatedTemperature = temperatureInterp;
    oneCoreTable.SalinityFile = repmat(matchTable.SalinityFile(i), height(oneCoreTable), 1);
    oneCoreTable.TemperatureFile = repmat(matchTable.TemperatureFile(i), height(oneCoreTable), 1);
    oneCoreTable.SalinityIceThickness = repmat(salinityIceThickness, height(oneCoreTable), 1);
    oneCoreTable.TemperatureIceThickness = repmat(temperatureIceThickness, height(oneCoreTable), 1);

    allCoreValuesTable = [allCoreValuesTable; oneCoreTable];

    oneAverage = table( ...
        matchTable.MatchKey(i), ...
        matchTable.DensityFile(i), ...
        matchTable.SalinityFile(i), ...
        matchTable.TemperatureFile(i), ...
        mean(densityData.Density, 'omitnan'), ...
        mean(densityData.LabTemperature, 'omitnan'), ...
        mean(densityData.DensitySalinity, 'omitnan'), ...
        mean(salinityInterp, 'omitnan'), ...
        mean(temperatureInterp, 'omitnan'), ...
        densityIceThickness, ...
        salinityIceThickness, ...
        temperatureIceThickness, ...
        densityFreeboard, ...
        maxDensityDepth, ...
        maxSalinityDepth, ...
        maxTemperatureDepth, ...
        'VariableNames', { ...
            'MatchKey', ...
            'DensityFile', ...
            'SalinityFile', ...
            'TemperatureFile', ...
            'MeanDensity', ...
            'MeanLabTemperature', ...
            'MeanDensitySalinity', ...
            'MeanInterpolatedSalinity', ...
            'MeanInterpolatedTemperature', ...
            'DensityIceThickness', ...
            'SalinityIceThickness', ...
            'TemperatureIceThickness', ...
            'DensityFreeboard', ...
            'MaxDensityDepth', ...
            'MaxSalinityDepth', ...
            'MaxTemperatureDepth'} ...
    );

    coreAverageTable = [coreAverageTable; oneAverage];

end

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

valuesOutFile = fullfile(outputDir, ...
    'processed_density_with_salinity_temperature_values.csv');

averagesOutFile = fullfile(outputDir, ...
    'processed_density_with_salinity_temperature_averages.csv');

writetable(allCoreValuesTable, valuesOutFile);
writetable(coreAverageTable, averagesOutFile);

fprintf('\nSaved:\n');
fprintf('%s\n', valuesOutFile);
fprintf('%s\n', averagesOutFile);

%% Helpers
function T = readDensityCore(filePath)

raw = readcell(filePath, 'DatetimeType', 'text');

fileName = string(getFileName(filePath));
folder = string(fileparts(filePath));

time = parseTime(raw{3,2});
lon = toNumber(raw{4,2});
lat = toNumber(raw{5,2});
freeboard = toNumber(raw{17,2});
snowThickness = toNumber(raw{18,2});
iceThickness = toNumber(raw{19,2});

headerRow = find(strcmpi(string(raw(:,1)), "min_depth"), 1, 'first');

if isempty(headerRow)
    T = table();
    return
end

dataStartRow = headerRow + 1;

minDepth = cellToNumber(raw(dataStartRow:end, 1));
maxDepth = cellToNumber(raw(dataStartRow:end, 2));

labTemperature = cellToNumberKeepNaN(raw(dataStartRow:end, 3));
densitySalinity = cellToNumberKeepNaN(raw(dataStartRow:end, 4));

density = cellToNumber(raw(dataStartRow:end, 11));

n = min([numel(minDepth), numel(maxDepth), numel(density)]);

minDepth = minDepth(1:n);
maxDepth = maxDepth(1:n);
density = density(1:n);

labTemperature = padOrTrim(labTemperature, n);
densitySalinity = padOrTrim(densitySalinity, n);

T = table( ...
    repmat(fileName, n, 1), ...
    repmat(folder, n, 1), ...
    repmat(time, n, 1), ...
    repmat(lon, n, 1), ...
    repmat(lat, n, 1), ...
    repmat(snowThickness, n, 1), ...
    repmat(freeboard, n, 1), ...
    repmat(iceThickness, n, 1), ...
    minDepth(:), ...
    maxDepth(:), ...
    labTemperature(:), ...
    densitySalinity(:), ...
    density(:), ...
    'VariableNames', { ...
        'DensityFile', ...
        'DensityFolder', ...
        'Time', ...
        'Lon', ...
        'Lat', ...
        'SnowThickness', ...
        'Freeboard', ...
        'IceThickness', ...
        'MinDepth', ...
        'MaxDepth', ...
        'LabTemperature', ...
        'DensitySalinity', ...
        'Density'} ...
);

end

function T = readSalinityCore(filePath)

raw = readcell(filePath, 'DatetimeType', 'text');

iceThickness = toNumber(raw{19,2});

headerRow = find(strcmpi(string(raw(:,1)), "min_depth"), 1, 'first');

if isempty(headerRow)
    T = table();
    return
end

header = string(raw(headerRow,:));
salinityCol = find(strcmpi(strtrim(header), "salinity"), 1, 'first');

if isempty(salinityCol)
    salinityCol = find(contains(lower(strtrim(header)), "salinity"), 1, 'first');
end

if isempty(salinityCol)
    T = table();
    return
end

dataStartRow = headerRow + 1;

minDepth = cellToNumber(raw(dataStartRow:end, 1));
maxDepth = cellToNumber(raw(dataStartRow:end, 2));
salinity = cellToNumber(raw(dataStartRow:end, salinityCol));

n = min([numel(minDepth), numel(maxDepth), numel(salinity)]);

T = table( ...
    minDepth(1:n), ...
    maxDepth(1:n), ...
    salinity(1:n), ...
    repmat(iceThickness, n, 1), ...
    'VariableNames', ...
    {'MinDepth','MaxDepth','Salinity','IceThickness'} ...
);

end

function T = readTemperatureCore(filePath)

raw = readcell(filePath, 'DatetimeType', 'text');

iceThickness = toNumber(raw{19,2});

nRows = size(raw, 1);
headerRow = [];
tempCol = [];

for r = 1:nRows
    thisRow = string(raw(r,:));
    tempColHere = find(strcmpi(strtrim(thisRow), "temperature"), 1, 'first');

    if ~isempty(tempColHere)
        headerRow = r;
        tempCol = tempColHere;
        break
    end
end

if isempty(headerRow)
    T = table();
    return
end

dataStartRow = headerRow + 1;

depth = cellToNumber(raw(dataStartRow:end, 1));
temperature = cellToNumber(raw(dataStartRow:end, tempCol));

n = min([numel(depth), numel(temperature)]);

T = table( ...
    depth(1:n), ...
    temperature(1:n), ...
    repmat(iceThickness, n, 1), ...
    'VariableNames', ...
    {'Depth','Temperature','IceThickness'} ...
);
end

function x = toNumber(value)

if isnumeric(value)
    x = value;
else
    txt = strtrim(string(value));

    if txt == "" || strcmpi(txt, "NA") || strcmpi(txt, "NaN") || strcmpi(txt, "missing")
        x = NaN;
    else
        x = str2double(txt);
    end
end

end

function values = cellToNumber(cellsIn)

values = nan(size(cellsIn));

for j = 1:numel(cellsIn)
    values(j) = toNumber(cellsIn{j});
end

values = values(~isnan(values));
values = values(:);

end

function values = cellToNumberKeepNaN(cellsIn)

values = nan(size(cellsIn));

for j = 1:numel(cellsIn)
    values(j) = toNumber(cellsIn{j});
end

values = values(:);

end

function values = padOrTrim(values, n)

values = values(:);

if numel(values) >= n
    values = values(1:n);
else
    values(end+1:n, 1) = NaN;
end

end

function t = parseTime(value)

if isdatetime(value)
    t = value;
    return
end

txt = char(strtrim(string(value)));

try
    t = datetime(txt, 'InputFormat', 'dd/MM/yyyy', 'Format', 'dd/MM/yyyy');
    return
catch
end

try
    t = datetime(txt, 'InputFormat', 'dd-MM-yyyy', 'Format', 'dd/MM/yyyy');
    return
catch
end

try
    t = datetime(txt, 'InputFormat', 'dd-MM-yyyy''T''HH:mm:ss.SS');
    return
catch
end

try
    t = datetime(txt, 'InputFormat', 'dd-MM-yyyy''T''HH:mm:ss');
    return
catch
end

t = NaT;

end

function fileName = getFileName(filePath)

[~, name, ext] = fileparts(filePath);
fileName = name + ext;

end

function paths = findFilePath(rootDir, fileNames)

paths = strings(size(fileNames));

for i = 1:numel(fileNames)

    if fileNames(i) == ""
        paths(i) = "";
        continue
    end

    found = dir(fullfile(rootDir, '**', char(fileNames(i))));

    if isempty(found)
        warning('Could not find file: %s', fileNames(i));
        paths(i) = "";
    else
        paths(i) = string(fullfile(found(1).folder, found(1).name));
    end

end

end

end
function loopless_fva(directoryPath, fractionOfOptimum)

% Get a list of all files in the directory with .xml or .sbml extension
files = dir(fullfile(directoryPath, '*.xml'));

% Loop through each file
for file = files'
    inputPath = fullfile(directoryPath, file.name);
    
    if endsWith(file.name, '_loopless.xml') || contains(file.name, "P3") || contains(file.name, "P5") || contains(file.name, "fastcore")  ...
        || contains(file.name, "N3") || contains(file.name, "N5")
        continue;
    end

    
    [modelDir, modelName, ext] = fileparts(inputPath);
    disp(inputPath)
    model = readCbModel(inputPath);
    %index = find(strcmp(model.rxns, 'EX_C00009__dra'));
    %model.lb(index) = -0.06;
    disp(optimizeCbModel(model));

    [minFlux, maxFlux] = fluxVariability(model, fractionOfOptimum, 'max', [model.rxns], 1, 'fastSNP');

    numReactions = length(model.rxns);

    for i = 1:numReactions
        model.lb(i) = minFlux(i); % Update the lower bound of the i-th reaction
        model.ub(i) = maxFlux(i); % Update the upper bound of the i-th reaction
    end

    disp(optimizeCbModel(model));

    outputFileName = [modelName, '_loopless', ext];
    outputPath = fullfile(modelDir, outputFileName);

    writeCbModel(model, 'format', 'sbml', 'fileName', outputPath);
end

end

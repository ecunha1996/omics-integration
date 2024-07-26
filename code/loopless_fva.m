function loopless_fva(inputPath, franctionOfOptimum)

[modelDir, modelName, ext] = fileparts(inputPath);

model = readCbModel(inputPath);

model = changeRxnBounds(model, 'EX_C00205__dra', -20, 'l');
model = changeRxnBounds(model, 'EX_C00011__dra', -2, 'l');

disp(optimizeCbModel(model));

[minFlux, maxFlux] = fluxVariability(model, franctionOfOptimum, 'max', [model.rxns], 1, 'fastSNP');

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
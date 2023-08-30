function [trainedClassifier30BilayeredNeuralNet, validationAccuracy] = trainClassifier30BilayeredNeuralNet(trainingData)
% [trainedClassifier30BilayeredNeuralNet, validationAccuracy] = trainClassifier30BilayeredNeuralNet(trainingData)
% Returns a trained Classifier30BilayeredNeuralNet and its accuracy. This code recreates the
% classification model trained in Classification Learner app. Use the
% generated code to automate training the same model with new data, or to
% learn how to programmatically train models.
%
%  Input:
%      trainingData: A table containing the same predictor and response
%       columns as those imported into the app.
%
%
%  Output:
%      trainedClassifier30BilayeredNeuralNet: A struct containing the trained Classifier30BilayeredNeuralNet. The
%       struct contains various fields with information about the trained
%       Classifier30BilayeredNeuralNet.
%
%      trainedClassifier30BilayeredNeuralNet.predictFcn: A function to make predictions on new
%       data.
%
%      validationAccuracy: A double representing the validation accuracy as
%       a percentage. In the app, the Models pane displays the validation
%       accuracy for each model.
%
% Use the code to train the model with new data. To retrain your
% Classifier30BilayeredNeuralNet, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a Classifier30BilayeredNeuralNet trained with the original data set
% T, enter:
%   [trainedClassifier30BilayeredNeuralNet, validationAccuracy] = trainClassifier30BilayeredNeuralNet(T)
%
% To make predictions with the returned 'trainedClassifier30BilayeredNeuralNet' on new data T2,
% use
%   [yfit,scores] = trainedClassifier30BilayeredNeuralNet.predictFcn(T2)
%
% T2 must be a table containing at least the same predictor columns as used
% during training. For details, enter:
%   trainedClassifier30BilayeredNeuralNet.HowToPredict

% Auto-generated by MATLAB on 29-Aug-2023 22:25:19


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = {'Circularity15', 'Circularity16', 'Circularity20', 'Circularity27', 'AR7', 'AR8', 'Roundness25', 'Roundness27', 'Roundness30', 'Solidity11', 'Solidity13', 'Solidity14', 'Eccentricity5', 'Eccentricity12', 'Eccentricity13', 'Eccentricity16', 'ESF25', 'ESF27', 'ESF30', 'ElongationFW12', 'ElongationFW13', 'ElongationFW17', 'ElongationFmF15', 'ElongationFmF24', 'ElongationFmF26', 'ElongationFmF29', 'Convexity7', 'Convexity17', 'Convexity18', 'Convexity19', 'Curl7', 'Curl14', 'Curl20', 'Curl28', 'Curl30', 'CurlW5', 'CurlW15', 'CurlW16', 'CurlW20', 'AreaNorm15', 'AreaNorm20', 'AreaNorm23', 'AreaNorm24', 'AreaNorm30', 'PerimeterNorm16', 'PerimeterNorm18', 'PerimeterNorm20', 'PerimeterNorm26', 'MajorNorm12', 'MajorNorm17', 'MajorNorm18', 'MajorNorm19', 'MinorNorm10', 'MinorNorm15', 'MinorNorm18', 'MinorNorm21', 'HeightNorm11', 'HeightNorm12', 'HeightNorm16', 'HeightNorm19', 'HeightNorm20', 'HeightNorm21', 'HeightNorm22', 'WidthNorm11', 'WidthNorm18', 'WidthNorm19', 'WidthNorm20', 'WidthNorm21', 'WidthNorm22', 'FeretNorm12', 'FeretNorm16', 'FeretNorm17', 'FeretNorm18', 'MinFeretNorm11', 'MinFeretNorm15', 'MinFeretNorm18', 'MinFeretNorm20'};
predictors = inputTable(:, predictorNames);
response = inputTable.Group;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];
classNames = categorical({'AA'; 'ABeta'; 'AS'; 'SCD'});

% Train a Classifier30BilayeredNeuralNet
% This code specifies all the Classifier30BilayeredNeuralNet options and trains the Classifier30BilayeredNeuralNet.
classificationNeuralNetwork = fitcnet(...
    predictors, ...
    response, ...
    'LayerSizes', [10 10], ...
    'Activations', 'relu', ...
    'Lambda', 0, ...
    'IterationLimit', 1000, ...
    'Standardize', true, ...
    'ClassNames', classNames);

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
neuralNetworkPredictFcn = @(x) predict(classificationNeuralNetwork, x);
trainedClassifier30BilayeredNeuralNet.predictFcn = @(x) neuralNetworkPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier30BilayeredNeuralNet.RequiredVariables = {'AR7', 'AR8', 'AreaNorm15', 'AreaNorm20', 'AreaNorm23', 'AreaNorm24', 'AreaNorm30', 'Circularity15', 'Circularity16', 'Circularity20', 'Circularity27', 'Convexity17', 'Convexity18', 'Convexity19', 'Convexity7', 'Curl14', 'Curl20', 'Curl28', 'Curl30', 'Curl7', 'CurlW15', 'CurlW16', 'CurlW20', 'CurlW5', 'ESF25', 'ESF27', 'ESF30', 'Eccentricity12', 'Eccentricity13', 'Eccentricity16', 'Eccentricity5', 'ElongationFW12', 'ElongationFW13', 'ElongationFW17', 'ElongationFmF15', 'ElongationFmF24', 'ElongationFmF26', 'ElongationFmF29', 'FeretNorm12', 'FeretNorm16', 'FeretNorm17', 'FeretNorm18', 'HeightNorm11', 'HeightNorm12', 'HeightNorm16', 'HeightNorm19', 'HeightNorm20', 'HeightNorm21', 'HeightNorm22', 'MajorNorm12', 'MajorNorm17', 'MajorNorm18', 'MajorNorm19', 'MinFeretNorm11', 'MinFeretNorm15', 'MinFeretNorm18', 'MinFeretNorm20', 'MinorNorm10', 'MinorNorm15', 'MinorNorm18', 'MinorNorm21', 'PerimeterNorm16', 'PerimeterNorm18', 'PerimeterNorm20', 'PerimeterNorm26', 'Roundness25', 'Roundness27', 'Roundness30', 'Solidity11', 'Solidity13', 'Solidity14', 'WidthNorm11', 'WidthNorm18', 'WidthNorm19', 'WidthNorm20', 'WidthNorm21', 'WidthNorm22'};
trainedClassifier30BilayeredNeuralNet.ClassificationNeuralNetwork = classificationNeuralNetwork;
trainedClassifier30BilayeredNeuralNet.About = 'This struct is a trained model exported from Classification Learner R2023a.';
trainedClassifier30BilayeredNeuralNet.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  [yfit,scores] = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = {'Circularity15', 'Circularity16', 'Circularity20', 'Circularity27', 'AR7', 'AR8', 'Roundness25', 'Roundness27', 'Roundness30', 'Solidity11', 'Solidity13', 'Solidity14', 'Eccentricity5', 'Eccentricity12', 'Eccentricity13', 'Eccentricity16', 'ESF25', 'ESF27', 'ESF30', 'ElongationFW12', 'ElongationFW13', 'ElongationFW17', 'ElongationFmF15', 'ElongationFmF24', 'ElongationFmF26', 'ElongationFmF29', 'Convexity7', 'Convexity17', 'Convexity18', 'Convexity19', 'Curl7', 'Curl14', 'Curl20', 'Curl28', 'Curl30', 'CurlW5', 'CurlW15', 'CurlW16', 'CurlW20', 'AreaNorm15', 'AreaNorm20', 'AreaNorm23', 'AreaNorm24', 'AreaNorm30', 'PerimeterNorm16', 'PerimeterNorm18', 'PerimeterNorm20', 'PerimeterNorm26', 'MajorNorm12', 'MajorNorm17', 'MajorNorm18', 'MajorNorm19', 'MinorNorm10', 'MinorNorm15', 'MinorNorm18', 'MinorNorm21', 'HeightNorm11', 'HeightNorm12', 'HeightNorm16', 'HeightNorm19', 'HeightNorm20', 'HeightNorm21', 'HeightNorm22', 'WidthNorm11', 'WidthNorm18', 'WidthNorm19', 'WidthNorm20', 'WidthNorm21', 'WidthNorm22', 'FeretNorm12', 'FeretNorm16', 'FeretNorm17', 'FeretNorm18', 'MinFeretNorm11', 'MinFeretNorm15', 'MinFeretNorm18', 'MinFeretNorm20'};
predictors = inputTable(:, predictorNames);
response = inputTable.Group;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];
classNames = categorical({'AA'; 'ABeta'; 'AS'; 'SCD'});

% Perform cross-validation
partitionedModel = crossval(trainedClassifier30BilayeredNeuralNet.ClassificationNeuralNetwork, 'KFold', 10);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

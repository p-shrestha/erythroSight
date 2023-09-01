% summarizeMetrics.m
% Program to read results for classification and generate summaries
% ps UBC 2023

%% Rank classifiers based on AUROC or testing accuracy
tic
classifierNamesAll = ["1FineTree", "2MediumTree", "3CoarseTree", "6EffLogReg", "7EffLinSVM", "9KernelNaiveBayes", "10LinSVM", "11QuadSVM", "12CubicSVM", "13FineGaussSVM", "14MedGaussSVM", "15CoarseGaussSVM", "16FineKNN", "17MedKNN", "18CoarseKNN", "19CosineKNN", "20CubicKNN", "21WeightedKNN", "22BoostedTree", "23BaggedTree", "24SubspaceDis", "25SubspaceKNN", "26RUSBoostedTree", "27NarrowNeuralNet", "28MedNeuralNet", "29WideNeuralNet", "30BilayeredNeuralNet", "31TrilayeredNeuralNet", "32SVMKernel", "33LogisticRegKernel"];


accTestMeanList = zeros(length(classifierNamesAll),1);
aurocMeanList = zeros(length(classifierNamesAll),1);
classifiersList = strings([length(classifierNamesAll),1]);

iCvalid = 0; % Number of classifiers that have results
for iC = 1:1:length(classifierNamesAll)
    if isfolder(classifierNamesAll(iC))
        iCvalid = iCvalid + 1;

        load(strcat(classifierNamesAll(iC),'/evalMetrics.mat'));
        accTestMeanList(iCvalid) = accTestMean;
        aurocMeanList(iCvalid) = aurocMean;
        classifiersList(iCvalid) = classifierNamesAll{iC};


    else
        fprintf('Classifier %s does not have results \n', classifierNamesAll(iC));
    end
end

% Remove zeros and empty strings
accTestMeanList(iCvalid+1:end,:) = [];
aurocMeanList(iCvalid+1:end,:) = [];
classifiersList(iCvalid+1:end,:) = [];

% Sort list in descending order
[accTestSortedList,idxAcc] = sort(accTestMeanList, 'descend');
[aurocSortedList,idxAUROC] = sort(aurocMeanList, 'descend');
classifierListAccTest = classifiersList(idxAcc);
classifierListAUROC = classifiersList(idxAUROC);

% Create a table with the ranked classifiers based on AUROC and accuracy
rankedTable = table(classifierListAUROC,aurocSortedList,classifierListAccTest, accTestSortedList, 'VariableNames',{'Classifiers (ranked by AUROC)', 'Mean AUROC', 'Classifiers (ranked by accuracy)', 'Mean testing accuracy'});
%rankedTable.Properties.Column = {'Classifiers (ranked by AUROC)', 'Mean AUROC', 'Classifiers (ranked by accuracy)', 'Mean testing accuracy'};

folderName = '0Summary';

% Create a folder 
if not(isfolder(folderName))
    mkdir(folderName);
end

% Save tables
fprintf("Saving tables as .csv files .... \n");
writetable(rankedTable, strcat(folderName,'/rankedTable.csv'));

toc

%% Create a table with a summary of overall evaluation metrics and save it

tic
fprintf("Creating overall metric table \n");

group = {'AA', 'ABeta','AS','SCD'};

% Save a table with only overall parameters

metricNames1 = ["AUROC", "F1-Score","Accuracy (OvA)", "Sensitivity", "Specificity", "Positive Predictive Value", "Negative Predictive Value"];
metricNamesAll = ["AUROC (macro-averaged)", "F1-Score (macro-averaged)","Accuracy OvA (macro-averaged)", "Sensitivity (macro-averaged)", "Specificity (macro-averaged)", "Positive Predictive Value (macro-averaged)", "Negative Predictive Value (macro-averaged)", "Validation Accuracy", "Testing Accuracy"];

nVar = length(metricNamesAll);
overallMetricsTable = array2table(strings(length(classifierListAUROC),nVar));

% iC is the classifier name
for iC = 1:1:length(classifierListAUROC)
    load(strcat(classifierListAUROC(iC),'/evalMetrics.mat'));

    % ii is the index of variable from metricNames1: AUROC, F1-Score, Acc(OVA), Sens, Spec, PPV, NPV
    for ii = 1:1:length(metricNames1)
        overallMetricsTable.(ii)(iC) = strcat(num2str(evalMetricsSummary.(17)(ii), 3), ' (', num2str(evalMetricsSummary.(19)(ii), 3), '-', num2str(evalMetricsSummary.(20)(ii), 3), ')');

    end
    % Save validation accuracy
    iRowSource = 1;
    overallMetricsTable.(ii+1)(iC) = strcat(num2str(overallMetricsSummary.(1)(iRowSource), 3), ' (', num2str(overallMetricsSummary.(3)(iRowSource), 3), '-', num2str(overallMetricsSummary.(4)(iRowSource), 3), ')');
    % Save testing accuracy
    iRowSource = 2;
    overallMetricsTable.(ii+2)(iC) = strcat(num2str(overallMetricsSummary.(1)(iRowSource), 3), ' (', num2str(overallMetricsSummary.(3)(iRowSource), 3), '-', num2str(overallMetricsSummary.(4)(iRowSource), 3), ')');

end

overallMetricsTable.Properties.RowNames = classifierListAUROC;
overallMetricsTable.Properties.VariableNames = metricNamesAll;

% Save tables
fprintf("Saving tables as .csv files .... \n");
writetable(overallMetricsTable, strcat(folderName,'/overallMetricsTable.csv'),'WriteRowNames',true);
toc

%% Create a table with a summary of individual evaluation metrics and save it

tic
fprintf("Creating evaluation metric table (individual and macro-averaged) \n");
%group = {'AA', 'ABeta','AS','SCD'};

% Save a table with evaluation metrics (groupwise and macroaveraged) 

metricNames1 = ["AUROC", "F1-Score","Accuracy (OvA)", "Sensitivity", "Specificity", "Positive Predictive Value", "Negative Predictive Value"];
metricNamesAll = ["AUROC", "F1-Score","Accuracy (OvA)", "Sensitivity", "Specificity", "Positive Predictive Value", "Negative Predictive Value", "Validation Accuracy", "Testing Accuracy"];

nVar = length(metricNamesAll);
evalMetricsTable = array2table(strings(length(classifierListAUROC)*6,nVar));

% Six rows per classifier for: 
% 1: Classifier name as row name and empty rows
% 2: Macro-averaged (or overall)
% 3: AA
% 4: ABeta
% 5: AS
% 6: SCD
rowNamesRepeat = ["Macro-averaged (or overall)"; "AA"; "ABeta"; "AS"; "SCD"];
nRowRepeat = length(rowNamesRepeat);
nGroup = length(group);
rowNames = strings(length(classifierListAUROC)*6,1);


% iC is the classifier name
for iC = 1:1:length(classifierListAUROC)

    % Enter the name of the classifier and the rowNamesRepeat to rowNames
    rowNames((nRowRepeat+1)*(iC-1)+1,1) = classifierListAUROC(iC);
    %rowNames(((nRowRepeat+1)*(iC-1)+2):((nRowRepeat+1)*(iC-1)+6),1) = rowNamesRepeat;
    % For having different rowNames use the following
    for rr = 1:1:length(rowNamesRepeat)
        rowNames((nRowRepeat+1)*(iC-1)+rr+1,1) = strcat(rowNamesRepeat(rr),' (', classifierListAUROC(iC),')');
    end

    load(strcat(classifierListAUROC(iC),'/evalMetrics.mat'));

    % ii is the index of variable from metricNames1: AUROC, F1-Score, Acc(OVA), Sens, Spec, PPV, NPV
    for ii = 1:1:length(metricNames1)
        % Add overall or macroaveraged values
        evalMetricsTable.(ii)((nRowRepeat+1)*(iC-1)+2) = strcat(num2str(evalMetricsSummary.(17)(ii), 3), ' (', num2str(evalMetricsSummary.(19)(ii), 3), '-', num2str(evalMetricsSummary.(20)(ii), 3), ')');

        % Add for rest of the individual groups
        for iG = 1:1:nGroup
            evalMetricsTable.(ii)((nRowRepeat+1)*(iC-1)+2+iG) = strcat(num2str(evalMetricsSummary.((nGroup*(iG-1))+1)(ii), 3), ' (', num2str(evalMetricsSummary.((nGroup*(iG-1))+3)(ii), 3), '-', num2str(evalMetricsSummary.((nGroup*(iG-1))+4)(ii), 3), ')');

        end

    end
    % Save validation accuracy
    iRowSource = 1;
    evalMetricsTable.(ii+1)((nRowRepeat+1)*(iC-1)+2) = strcat(num2str(overallMetricsSummary.(1)(iRowSource), 3), ' (', num2str(overallMetricsSummary.(3)(iRowSource), 3), '-', num2str(overallMetricsSummary.(4)(iRowSource), 3), ')');
    % Save testing accuracy
    iRowSource = 2;
    evalMetricsTable.(ii+2)((nRowRepeat+1)*(iC-1)+2) = strcat(num2str(overallMetricsSummary.(1)(iRowSource), 3), ' (', num2str(overallMetricsSummary.(3)(iRowSource), 3), '-', num2str(overallMetricsSummary.(4)(iRowSource), 3), ')');

end

evalMetricsTable.Properties.RowNames = rowNames;
evalMetricsTable.Properties.VariableNames = metricNamesAll;

% Save tables
fprintf("Saving tables as .csv files .... \n");
writetable(evalMetricsTable, strcat(folderName,'/evalMetricsTable.csv'),'WriteRowNames',true);

%% Save variables
%rankedTable, classifierListAccTest, accTestSortedList,
%classifierListAUROC, aurocSortedList, overallMetricsTable, evalMetricsTable

fprintf("Saving the summary variables/tables summaryVariables.mat .... \n");
save(strcat(folderName,'/summaryVariables.mat'), "rankedTable", "classifierListAccTest", "accTestSortedList", "classifierListAUROC", "aurocSortedList", "overallMetricsTable", "evalMetricsTable", '-v7.3');


toc

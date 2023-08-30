% classificationMultipleSplits_MisclassCost10_CS.m
% Testing classification of data into 4 groups using different models
% Added misclassification cost of 5 to AA-ABeta and AA-AS
% ps UBC 2023

% Can be run in linux or the cluster with the following command (example): 
% matlab -nodisplay -r "classifierName = \"12CubicSVM\";classificationMultipleSplits_MisclassCost10_CS"

% ------------------------------------------------------------------------
%                     IMPORT/LOAD DATA
% ------------------------------------------------------------------------


%% 1) Import histData for processing of morphological parameters
% Note: Here, the histDataCS (or coverslip-wise) histcount data is selected
% (For finer analysis, image-wise data or groups of image-wise data can be
% selected, for coarser analysis, donor-wise data can be selected

% Load all histData - [Note: Need to modify path based on where the data is
% stored in local machine or cluster]
tic 
fprintf('Reading data - generalInfo.mat and histData.mat \n');

group = {'AA', 'ABeta','AS','SCD'};
nGroup = length(group);

for i = 1:1:length(group)
    % Run similar command for all groups to read generalInfo.mat and
    % histData.mat
    % Note: Adjust file path if required
    % Example: AA.info = load('../../../1DCharacterization/t02hrs/cnCombined/30b40mp/AA/generalInfo.mat');
    % Example: AA.hD = load('../../../1DCharacterization/t02hrs/cnCombined/30b40mp/AA/histData.mat');
    eval(strcat(group{i},".info = load('../../../../Plots/1DCharacterization/t02hrs/cnCombined/30b40mp/",group{i},"/generalInfo.mat');"));
    eval(strcat(group{i},".hD = load('../../../../Plots/1DCharacterization/t02hrs/cnCombined/30b40mp/",group{i},"/histData.mat');"));

end

toc

classifierNamesAll = ["1FineTree", "2MediumTree", "3CoarseTree", "6EffLogReg", "7EffLinSVM", "9KernelNaiveBayes", "10LinSVM", "11QuadSVM", "12CubicSVM", "13FineGaussSVM", "14MedGaussSVM", "15CoarseGaussSVM", "16FineKNN", "17MedKNN", "18CoarseKNN", "19CosineKNN", "20CubicKNN", "21WeightedKNN", "22BoostedTree", "23BaggedTree", "24SubspaceDis", "25SubspaceKNN", "26RUSBoostedTree", "27NarrowNeuralNet", "28MedNeuralNet", "29WideNeuralNet", "30BilayeredNeuralNet", "31TrilayeredNeuralNet", "32SVMKernel", "33LogisticRegKernel"];

% Note: The classfierName can be supplied as an input in the cluster
%classifierName = "12CubicSVM";
%classifierName = "11QuadSVM";% "24SubspaceDis"; %"29WideNeuralNet"; %
%classifierName = classifierNamesAll(9); %9 = 12CubicSVM

nIterations = 1000;

accVal = zeros(nIterations,1);  % Validation accuracy for each training iteration
accTest = zeros(nIterations,1); % Testing accuracy for each testing iteration
labelsGT = cell(1,nIterations); % Ground truth labels for each testing iteration
labelsPred = cell(1,nIterations);   % Predicted labels for each testing iteration
scores = cell(1,nIterations); % Class scores outputted from the model testing
nTrainPreBal = zeros(nIterations, nGroup);  % Number of samples for training before balancing classes for each group
nTrainPostBal = zeros(nIterations, nGroup);  % Number of samples for training after balancing classes for each group
nTestPreBal = zeros(nIterations, nGroup);  % Number of samples for testing before balancing classes for each group
nTestPostBal = zeros(nIterations, nGroup);  % Number of samples for testing after balancing classes for each group

cMat = zeros(nGroup,nGroup,nIterations); % cMat is a 3D matrix storing all confusion matrices for all iterations

tp = zeros(nIterations, nGroup); % True positive for each class and testing iteration
fp = zeros(nIterations, nGroup); % False positive for each class and testing iteration
tn = zeros(nIterations, nGroup); % True negative for each class and testing iteration
fn = zeros(nIterations, nGroup); % False negative for each class and testing iteration
accOVAClass = zeros(nIterations, nGroup); % Accuracy (one-vs-all) for each class and testing iteration
sensClass = zeros(nIterations, nGroup); % Sensitivity for each class and testing iteration
specClass = zeros(nIterations, nGroup); % Specificity for each class and testing iteration
ppvClass = zeros(nIterations, nGroup); % Positive predictive value for each class and testing iteration
npvClass = zeros(nIterations, nGroup); % Negative predictive value for each class and testing iteration
f1scoreClass = zeros(nIterations, nGroup); % F1-score for each class and testing iteration

aurocClass = zeros(nIterations, nGroup); % Area under receiver operating characteristic curve (AUROC) for each class and testing iteration
aurocMacroAvg = zeros(nIterations, 1); % Macro-averaged accuracy (one-vs-all) testing iteration

accOVAMacroAvg = zeros(nIterations, 1); % Macro-averaged accuracy (one-vs-all) testing iteration
sensMacroAvg = zeros(nIterations, 1); % Macro-averaged sensitivity (one-vs-all) testing iteration
specMacroAvg = zeros(nIterations, 1); % Macro-averaged specificity (one-vs-all) testing iteration
ppvMacroAvg = zeros(nIterations, 1); % Macro-averaged PPV (one-vs-all) testing iteration
npvMacroAvg = zeros(nIterations, 1); % Macro-averaged NPV (one-vs-all) testing iteration
f1scoreMacroAvg = zeros(nIterations, 1); % Macro-averaged F1-score (one-vs-all) testing iteration

%{
% Parameters needed for calculating ROC curve for each iteration
AUCn = zeros(nIterations, length(group));
xROCn = cell(1,nIterations);
yROCn = cell(1,nIterations);
tROCn = cell(1,nIterations);
%}

% ------------------------------------------------------------------------
%       DATA PREPARATION (RANDOM SPLITS, BALANCING, ETC.)
% ------------------------------------------------------------------------

for Ni = 1:1:nIterations
    fprintf('%d of %d Iteration; Model = %s \n', Ni, nIterations, classifierName);
    %% 2) Randomly separate data into training and testing sets
    % Training:testing data split approximately 80:20
    % It is important that the split is participant or donor wise, such that
    % all the data (e.g. coverslips) from one participant only appears in
    % either trainig or testing data

    % The numbers of all groups are "Group (Total) = Training:Testing"
    % The following example is for a trainRatio of 0.8 (train:test = 80:20)
    % AA (30) = 24:6
    % ABeta (23) = 18:5
    % AS (45) = 36:9
    % SCD (40) = 32:8

    tic

    % Training dataset ratio
    trainRatio = 0.8;
    fprintf('\t Randomly separating data, with training:testing = %.f:%.f \n',trainRatio*100, (1-trainRatio)*100);

    % Calculate the number of participants to be separated for testing
    nTrAA = round(length(AA.hD.histDataCS30b40mp)*trainRatio);
    nTrABeta = round(length(ABeta.hD.histDataCS30b40mp)*trainRatio);
    nTrAS = round(length(AS.hD.histDataCS30b40mp)*trainRatio);
    nTrSCD = round(length(SCD.hD.histDataCS30b40mp)*trainRatio);

    % Generate random indices for all groups (full size)
    AArndI = randperm(length(AA.hD.histDataCS30b40mp));
    ABetarndI = randperm(length(ABeta.hD.histDataCS30b40mp));
    ASrndI = randperm(length(AS.hD.histDataCS30b40mp));
    SCDrndI = randperm(length(SCD.hD.histDataCS30b40mp));

    % Store values in training set
    trainAA = {AA.hD.histDataCS30b40mp{AArndI(1:nTrAA)}};
    trainABeta = {ABeta.hD.histDataCS30b40mp{ABetarndI(1:nTrABeta)}};
    trainAS = {AS.hD.histDataCS30b40mp{ASrndI(1:nTrAS)}};
    trainSCD = {SCD.hD.histDataCS30b40mp{SCDrndI(1:nTrSCD)}};

    % Store values in testing set
    testAA = {AA.hD.histDataCS30b40mp{AArndI(nTrAA+1:end)}};
    testABeta = {ABeta.hD.histDataCS30b40mp{ABetarndI(nTrABeta+1:end)}};
    testAS = {AS.hD.histDataCS30b40mp{ASrndI(nTrAS+1:end)}};
    testSCD = {SCD.hD.histDataCS30b40mp{SCDrndI(nTrSCD+1:end)}};

    %toc


    %% 3) Store training data into table
    % Load histData for each morphological parameter as a feature (in the form
    % of an array) and the group information related to that cover-slip or
    % donor. If each morphological parameter has 'nbins' bins, then each
    % parameter occupies 'nbins' rows, e.g. for Area: Area1, Area2, ....
    %tic

    morphParamName = AA.info.morphParamName;
    morphParamNameC = cell(1, length(morphParamName));
    %morphParamNameC = {'Area','Perimeter','Angle','Major','Minor','Height','Width','Mean','StdDev','Median','Skewness','Kurtosis','Feret','FeretAngle','MinFeret','FeretX','FeretY','Circularity','AR','Roundness','ConvexHullArea','Solidity','Eccentricity','ESF','ElongationFW','ElongationFmF','ConvexPerimeter','Convexity','FibreLength','FibreWidth','Curl','CurlW'};
    for mpn = 1:1:length(morphParamName)
        morphParamNameC{mpn} = morphParamName{mpn};
    end
    columnNames = [morphParamNameC, 'Group'];

    % List of morphological parameters to keep and list of parameters to delete
    keep19 = [18,19,20,22,23,24,25,26,28,31,32,33,34,35,36,37,38,39,40];
    delete21 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,21,27,29,30];

    tableRow = 1;

    %fprintf('Storing training data into table \n');
    fprintf('\t Storing training data into table ....');
    % Fill in the first row with the information from AA (first coverslip only)
    % The morphological parameters and the groupName are added as columns
    trainDataCN_4groups_30b19ndmp = cell2table([trainAA{1}{1}, char(AA.info.groupName)]);

    % Add the remaining rows
    % Loop through all donors or UIDs in the group AA
    %fprintf('\t Storing data from AA \n');
    fprintf('AA...');
    for uidI = 1:1:length(trainAA)
        % Loop through all coverslips of the donor
        for csI = 1:1:length(trainAA{1,uidI})
            if tableRow > 1
                trainDataCN_4groups_30b19ndmp = [trainDataCN_4groups_30b19ndmp; cell2table([trainAA{1,uidI}{1,csI}, char(AA.info.groupName)])];
            end
            tableRow = tableRow + 1;
        end
    end

    % Loop through all donors or UIDs in the group ABeta
    %fprintf('\t Storing data from ABeta \n');
    fprintf('ABeta...');
    for uidI = 1:1:length(trainABeta)
        % Loop through all coverslips of the donor
        for csI = 1:1:length(trainABeta{1,uidI})
            trainDataCN_4groups_30b19ndmp = [trainDataCN_4groups_30b19ndmp; cell2table([trainABeta{1,uidI}{1,csI}, char(ABeta.info.groupName)])];
            tableRow = tableRow + 1;
        end
    end

    % Loop through all donors or UIDs in the group AS
    %fprintf('\t Storing data from AS \n');
    fprintf('AS....');
    for uidI = 1:1:length(trainAS)
        % Loop through all coverslips of the donor
        for csI = 1:1:length(trainAS{1,uidI})
            trainDataCN_4groups_30b19ndmp = [trainDataCN_4groups_30b19ndmp; cell2table([trainAS{1,uidI}{1,csI}, char(AS.info.groupName)])];
            tableRow = tableRow + 1;
        end
    end

    % Loop through all donors or UIDs in the group SS
    %fprintf('\t Storing data from SCD (SBeta and SS) \n');
    fprintf('SCD (SBeta and SS) \n');
    for uidI = 1:1:length(trainSCD)
        % Loop through all coverslips of the donor
        for csI = 1:1:length(trainSCD{1,uidI})
            trainDataCN_4groups_30b19ndmp = [trainDataCN_4groups_30b19ndmp; cell2table([trainSCD{1,uidI}{1,csI}, char(SCD.info.groupName)])];
            tableRow = tableRow + 1;
        end
    end

    % After all additions, make tableRow equal to number of rows
    tableRow = tableRow -1;

    % Replace variable names with morphParamNameC
    trainDataCN_4groups_30b19ndmp.Properties.VariableNames = columnNames;
    % Make Group variable categorical
    trainDataCN_4groups_30b19ndmp.Group = categorical(trainDataCN_4groups_30b19ndmp.Group);
    % Add description
    trainDataCN_4groups_30b19ndmp.Properties.Description = 'Table with histDataCS30b40mp for coverslips from training set and 11 non-dimensional morphological parameters';

    % Delete the parameters that are dimensional or intensity based
    for ii = 1:1:length(delete21)
        eval(strcat('trainDataCN_4groups_30b19ndmp.',morphParamNameC{delete21(ii)},'=[];'));
    end
    %toc

    %% 4) Store testing data into table
    % Load histData for each morphological parameter as a feature (in the form
    % of an array) and the group information related to that cover-slip or
    % donor. If each morphological parameter has 'nbins' bins, then each
    % parameter occupies 'nbins' rows, e.g. for Area: Area1, Area2, ....
    %tic

    %{
morphParamName = AA.info.morphParamName;
morphParamNameC = cell(1, length(morphParamName));
%morphParamNameC = {'Area','Perimeter','Angle','Major','Minor','Height','Width','Mean','StdDev','Median','Skewness','Kurtosis','Feret','FeretAngle','MinFeret','FeretX','FeretY','Circularity','AR','Roundness','ConvexHullArea','Solidity','Eccentricity','ESF','ElongationFW','ElongationFmF','ConvexPerimeter','Convexity','FibreLength','FibreWidth','Curl','CurlW'};
for mpn = 1:1:length(morphParamName)
    morphParamNameC{mpn} = morphParamName{mpn};
end
columnNames = [morphParamNameC, 'Group'];
    %}

    tableRowTest = 1;

    fprintf('\t Storing testing data into table ....');
    % Fill in the first row with the information from AA (first coverslip only)
    % The morphological parameters and the groupName are added as columns
    testDataCN_4groups_30b19ndmp = cell2table([testAA{1}{1}, char(AA.info.groupName)]);

    % Add the remaining rows
    % Loop through all donors or UIDs in the group AA
    %fprintf('\t Storing data from AA \n');
    fprintf('AA...');
    for uidI = 1:1:length(testAA)
        % Loop through all coverslips of the donor
        for csI = 1:1:length(testAA{1,uidI})
            if tableRowTest > 1
                testDataCN_4groups_30b19ndmp = [testDataCN_4groups_30b19ndmp; cell2table([testAA{1,uidI}{1,csI}, char(AA.info.groupName)])];
            end
            tableRowTest = tableRowTest + 1;
        end
    end

    % Loop through all donors or UIDs in the group ABeta
    %fprintf('\t Storing data from ABeta \n');
    fprintf('ABeta...');
    for uidI = 1:1:length(testABeta)
        % Loop through all coverslips of the donor
        for csI = 1:1:length(testABeta{1,uidI})
            testDataCN_4groups_30b19ndmp = [testDataCN_4groups_30b19ndmp; cell2table([testABeta{1,uidI}{1,csI}, char(ABeta.info.groupName)])];
            tableRowTest = tableRowTest + 1;
        end
    end

    % Loop through all donors or UIDs in the group AS
    %fprintf('\t Storing data from AS \n');
    fprintf('AS...');
    for uidI = 1:1:length(testAS)
        % Loop through all coverslips of the donor
        for csI = 1:1:length(testAS{1,uidI})
            testDataCN_4groups_30b19ndmp = [testDataCN_4groups_30b19ndmp; cell2table([testAS{1,uidI}{1,csI}, char(AS.info.groupName)])];
            tableRowTest = tableRowTest + 1;
        end
    end

    % Loop through all donors or UIDs in the group SS
    %fprintf('\t Storing data from SCD (SBeta and SS) \n');
    fprintf('SCD (SBeta and SS) \n');
    for uidI = 1:1:length(testSCD)
        % Loop through all coverslips of the donor
        for csI = 1:1:length(testSCD{1,uidI})
            testDataCN_4groups_30b19ndmp = [testDataCN_4groups_30b19ndmp; cell2table([testSCD{1,uidI}{1,csI}, char(SCD.info.groupName)])];
            tableRowTest = tableRowTest + 1;
        end
    end

    % After all additions, make tableRow equal to number of rows
    tableRowTest = tableRowTest -1;

    % Replace variable names with morphParamNameC
    testDataCN_4groups_30b19ndmp.Properties.VariableNames = columnNames;
    % Make Group variable categorical
    testDataCN_4groups_30b19ndmp.Group = categorical(testDataCN_4groups_30b19ndmp.Group);
    % Add description
    testDataCN_4groups_30b19ndmp.Properties.Description = 'Table with histDataCS30b40mp for coverslips from testing set and 11 non-dimensional morphological parameters';

    % Delete the parameters that are dimensional or intensity based
    for ii = 1:1:length(delete21)
        eval(strcat('testDataCN_4groups_30b19ndmp.',morphParamNameC{delete21(ii)},'=[];'));
    end
    %toc

    %% 5) Save tables with training and testing data
    %{
tic
roundN = 'round3';

fprintf("Saving table trainDataCN_4groups_30b19ndmp.... \n");
save(strcat('randomSplits4Groups/',roundN,'/trainDataCN_4groups_30b19ndmp.mat'), 'trainDataCN_4groups_30b19ndmp');

fprintf("Saving table testDataCN_4groups_30b19ndmp.... \n");
save(strcat('randomSplits4Groups/',roundN,'/testDataCN_4groups_30b19ndmp.mat'), 'testDataCN_4groups_30b19ndmp');

toc
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          BALANCE CLASSES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% 5) Find number of all classes

    trainG = cellstr(trainDataCN_4groups_30b19ndmp.Group);
    testG = cellstr(testDataCN_4groups_30b19ndmp.Group);

    AATrainInd = find(contains(trainG,'AA')); 
    AATestInd = find(contains(testG,'AA')); 
    nTrainPreBal(Ni, 1) = length(AATrainInd);
    nTestPreBal(Ni, 1) = length(AATestInd);

    ABetaTrainInd = find(contains(trainG,'ABeta'));
    ABetaTestInd = find(contains(testG,'ABeta'));
    nTrainPreBal(Ni, 2) = length(ABetaTrainInd);
    nTestPreBal(Ni, 2) = length(ABetaTestInd);

    ASTrainInd = find(contains(trainG,'AS'));
    ASTestInd = find(contains(testG,'AS'));
    nTrainPreBal(Ni, 3) = length(ASTrainInd);
    nTestPreBal(Ni, 3) = length(ASTestInd);

    SCDTrainInd = find(contains(trainG,'SCD'));
    SCDTestInd = find(contains(testG,'SCD'));
    nTrainPreBal(Ni, 4) = length(SCDTrainInd);
    nTestPreBal(Ni, 4) = length(SCDTestInd);

    % fprintf('Training data numbers ..... \n');
    % fprintf('\t Number of AA in training data = %d \n', length(AATrainInd));
    % fprintf('\t Number of ABeta in training data = %d \n', length(ABetaTrainInd));
    % fprintf('\t Number of AS in training data = %d \n', length(ASTrainInd));
    % fprintf('\t Number of SCD in training data = %d \n', length(SCDTrainInd));
    fprintf('\t Training data numbers: AA = %d; ABeta = %d; AS = %d; SCD = %d \n', length(AATrainInd), length(ABetaTrainInd), length(ASTrainInd), length(SCDTrainInd));

    % fprintf('Testing data numbers ..... \n');
    % fprintf('\t Number of AA in testing data = %d \n', length(AATestInd));
    % fprintf('\t Number of ABeta in testing data = %d \n', length(ABetaTestInd));
    % fprintf('\t Number of AS in testing data = %d \n', length(ASTestInd));
    % fprintf('\t Number of SCD in testing data = %d \n', length(SCDTestInd));
    fprintf('\t Testing data numbers: AA = %d; ABeta = %d; AS = %d; SCD = %d \n', length(AATestInd), length(ABetaTestInd), length(ASTestInd), length(SCDTestInd));

    nTrain = length(AATrainInd)+length(ABetaTrainInd)+length(ASTrainInd)+length(SCDTrainInd);
    nTest = length(AATestInd) + length(ABetaTestInd) + length(ASTestInd) + length(SCDTestInd);
    % fprintf('Training data percentages ..... \n');
    % fprintf('\t Number of AA in training data = %.1f%% \n', length(AATrainInd)*100/nTrain);
    % fprintf('\t Number of ABeta in training data = %.1f%% \n', length(ABetaTrainInd)*100/nTrain);
    % fprintf('\t Number of AS in training data = %.1f%% \n', length(ASTrainInd)*100/nTrain);
    % fprintf('\t Number of SCD in training data = %.1f%% \n', length(SCDTrainInd)*100/nTrain);
    % 
    % fprintf('Testing data percentages ..... \n');
    % fprintf('\t Number of AA in testing data = %.1f%% \n', length(AATestInd)*100/nTest);
    % fprintf('\t Number of ABeta in testing data = %.1f%% \n', length(ABetaTestInd)*100/nTest);
    % fprintf('\t Number of AS in testing data = %.1f%% \n', length(ASTestInd)*100/nTest);
    % fprintf('\t Number of SCD in testing data = %.1f%% \n', length(SCDTestInd)*100/nTest);


    %% 6) Balance the groups by upsampling the minority data

    balancedtrainDataCN_4groups_30b19ndmp = trainDataCN_4groups_30b19ndmp;
    balancedtestDataCN_4groups_30b19ndmp = testDataCN_4groups_30b19ndmp;

    % Find which class has maximum (AS or SCD)
    if length(ASTrainInd) > length(SCDTrainInd)
        nMaxTrain = length(ASTrainInd);

    else
        nMaxTrain = length(SCDTrainInd);
    end

    if length(ASTestInd) > length(SCDTestInd)
        nMaxTest = length(ASTestInd);

    else
        nMaxTest = length(SCDTestInd);
    end

    % AA
    % Balance training set
    % Randomly select and copy a row for AA at the end of the table
    %fprintf('Upsampling AA class by randomly repeating %d rows from the training set \n', nMaxTrain - length(AATrainInd));
    for ii = 1:1:(nMaxTrain - length(AATrainInd))
        % Randomly choose an index from the AA list of indices
        randInd = randsample(AATrainInd,1);

        % New table (with one row) corresponding to the randomly sampled index
        t1 = trainDataCN_4groups_30b19ndmp(randInd,:);

        % Concatenate new table
        balancedtrainDataCN_4groups_30b19ndmp = [balancedtrainDataCN_4groups_30b19ndmp; t1];

    end

    % Balance testing set
    % Randomly select and copy a row for AA at the end of the table
    %fprintf('Upsampling AA class by randomly repeating %d rows from the testing set \n', nMaxTest - length(AATestInd));
    for ii = 1:1:(nMaxTest - length(AATestInd))
        % Randomly choose an index from the AA list of indices
        randInd = randsample(AATestInd,1);

        % New table (with one row) corresponding to the randomly sampled index
        t1 = testDataCN_4groups_30b19ndmp(randInd,:);

        % Concatenate new table
        balancedtestDataCN_4groups_30b19ndmp = [balancedtestDataCN_4groups_30b19ndmp; t1];

    end

    % ABeta
    % Balance training set
    % Randomly select and copy a row for ABeta at the end of the table
    %fprintf('Upsampling ABeta class by randomly repeating %d rows from the training set \n', nMaxTrain - length(ABetaTrainInd));
    for ii = 1:1:(nMaxTrain - length(ABetaTrainInd))
        % Randomly choose an index from the ABeta list of indices
        randInd = randsample(ABetaTrainInd,1);

        % New table (with one row) corresponding to the randomly sampled index
        t1 = trainDataCN_4groups_30b19ndmp(randInd,:);

        % Concatenate new table
        balancedtrainDataCN_4groups_30b19ndmp = [balancedtrainDataCN_4groups_30b19ndmp; t1];

    end

    % Balance testing set
    % Randomly select and copy a row for ABeta at the end of the table
    %fprintf('Upsampling ABeta class by randomly repeating %d rows from the testing set \n', nMaxTest - length(ABetaTestInd));
    for ii = 1:1:(nMaxTest - length(ABetaTestInd))
        % Randomly choose an index from the SBeta list of indices
        randInd = randsample(ABetaTestInd,1);

        % New table (with one row) corresponding to the randomly sampled index
        t1 = testDataCN_4groups_30b19ndmp(randInd,:);

        % Concatenate new table
        balancedtestDataCN_4groups_30b19ndmp = [balancedtestDataCN_4groups_30b19ndmp; t1];

    end

    if length(ASTrainInd) > length(SCDTrainInd)
        % SCD
        % Balance training set
        % Randomly select and copy a row for SCD at the end of the table
        %fprintf('Upsampling SCD class by randomly repeating %d rows from the training set \n', nMaxTrain - length(SCDTrainInd));
        for ii = 1:1:(nMaxTrain - length(SCDTrainInd))
            % Randomly choose an index from the SS list of indices
            randInd = randsample(SCDTrainInd,1);

            % New table (with one row) corresponding to the randomly sampled index
            t1 = trainDataCN_4groups_30b19ndmp(randInd,:);

            % Concatenate new table
            balancedtrainDataCN_4groups_30b19ndmp = [balancedtrainDataCN_4groups_30b19ndmp; t1];

        end
    else
        % AS
        % Balance training set
        % Randomly select and copy a row for SCD at the end of the table
        %fprintf('Upsampling AS class by randomly repeating %d rows from the training set \n', nMaxTrain - length(ASTrainInd));
        for ii = 1:1:(nMaxTrain - length(ASTrainInd))
            % Randomly choose an index from the SS list of indices
            randInd = randsample(ASTrainInd,1);

            % New table (with one row) corresponding to the randomly sampled index
            t1 = trainDataCN_4groups_30b19ndmp(randInd,:);

            % Concatenate new table
            balancedtrainDataCN_4groups_30b19ndmp = [balancedtrainDataCN_4groups_30b19ndmp; t1];

        end
    end


    if length(ASTestInd) > length(SCDTestInd)
        % Balance testing set
        % Randomly select and copy a row for SCD at the end of the table
        %fprintf('Upsampling SCD class by randomly repeating %d rows from the testing set \n', nMaxTest - length(SCDTestInd));
        for ii = 1:1:(nMaxTest - length(SCDTestInd))
            % Randomly choose an index from the SS list of indices
            randInd = randsample(SCDTestInd,1);

            % New table (with one row) corresponding to the randomly sampled index
            t1 = testDataCN_4groups_30b19ndmp(randInd,:);

            % Concatenate new table
            balancedtestDataCN_4groups_30b19ndmp = [balancedtestDataCN_4groups_30b19ndmp; t1];

        end
    else
        % Balance testing set
        % Randomly select and copy a row for AS at the end of the table
        %fprintf('Upsampling AS class by randomly repeating %d rows from the testing set \n', nMaxTest - length(ASTestInd));
        for ii = 1:1:(nMaxTest - length(ASTestInd))
            % Randomly choose an index from the SS list of indices
            randInd = randsample(ASTestInd,1);

            % New table (with one row) corresponding to the randomly sampled index
            t1 = testDataCN_4groups_30b19ndmp(randInd,:);

            % Concatenate new table
            balancedtestDataCN_4groups_30b19ndmp = [balancedtestDataCN_4groups_30b19ndmp; t1];

        end
    end



    %% 7) Check the number of SS and SBeta in the balanced set

    trainG1 = cellstr(balancedtrainDataCN_4groups_30b19ndmp.Group);
    testG1 = cellstr(balancedtestDataCN_4groups_30b19ndmp.Group);

    AATrainInd1 = find(contains(trainG1,'AA'));
    AATestInd1 = find(contains(testG1,'AA'));
    nTrainPostBal(Ni, 1) = length(AATrainInd1);
    nTestPostBal(Ni, 1) = length(AATestInd1);

    ABetaTrainInd1 = find(contains(trainG1,'ABeta'));
    ABetaTestInd1 = find(contains(testG1,'ABeta'));
    nTrainPostBal(Ni, 2) = length(ABetaTrainInd1);
    nTestPostBal(Ni, 2) = length(ABetaTestInd1);
    

    ASTrainInd1 = find(contains(trainG1,'AS'));
    ASTestInd1 = find(contains(testG1,'AS'));
    nTrainPostBal(Ni, 3) = length(ASTrainInd1);
    nTestPostBal(Ni, 3) = length(ASTestInd1);

    SCDTrainInd1 = find(contains(trainG1,'SCD'));
    SCDTestInd1 = find(contains(testG1,'SCD'));
    nTrainPostBal(Ni, 4) = length(SCDTrainInd1);
    nTestPostBal(Ni, 4) = length(SCDTestInd1);

    % fprintf('Training data numbers after balancing classes ..... \n');
    % fprintf('\t Number of AA in training data = %d \n', length(AATrainInd1));
    % fprintf('\t Number of ABeta in training data = %d \n', length(ABetaTrainInd1));
    % fprintf('\t Number of AS in training data = %d \n', length(ASTrainInd1));
    % fprintf('\t Number of SCD in training data = %d \n', length(SCDTrainInd1));
    fprintf('\t Training data numbers after balancing classes: AA = %d; ABeta = %d; AS = %d; SCD = %d \n', length(AATrainInd1), length(ABetaTrainInd1), length(ASTrainInd1), length(SCDTrainInd1));


    % fprintf('Testing data numbers after balancing classes ..... \n');
    % fprintf('\t Number of AA in testing data = %d \n', length(AATestInd1));
    % fprintf('\t Number of ABeta in testing data = %d \n', length(ABetaTestInd1));
    % fprintf('\t Number of AS in testing data = %d \n', length(ASTestInd1));
    % fprintf('\t Number of SCD in testing data = %d \n', length(SCDTestInd1));
    fprintf('\t Testing data numbers after balancing classes: AA = %d; ABeta = %d; AS = %d; SCD = %d \n', length(AATestInd1), length(ABetaTestInd1), length(ASTestInd1), length(SCDTestInd1));

    toc

    %% 8) Save the balanced classes
    %{
tic

fprintf("Saving table balancedtrainDataCN_4groups_30b19ndmp.... \n");
save('balancedtrainDataCN_4groups_30b19ndmp.mat', 'balancedtrainDataCN_4groups_30b19ndmp');

fprintf("Saving table balancedtestDataCN_4groups_30b19ndmp.... \n");
save('balancedtestDataCN_4groups_30b19ndmp.mat', 'balancedtestDataCN_4groups_30b19ndmp');

toc
    %}

    % ------------------------------------------------------------------------
    %                               TRAIN DATA
    % ------------------------------------------------------------------------

    %% Train data using Cubic SVM classifier
    tic
    fprintf('\t Training data using classifier: %s .....', classifierName);
    %[trainedClassifier, accVal(Ni)] = trainClassifier12CubicSVM(balancedtrainDataCN_4groups_30b19ndmp);
    [trainedClassifier, accVal(Ni)] = eval(strcat('trainClassifier',classifierName,'(balancedtrainDataCN_4groups_30b19ndmp);'));
    %accVal(Ni) = accVal(Ni)*100; % percentage
    toc

    % ------------------------------------------------------------------------
    %                               TEST DATA
    % ------------------------------------------------------------------------

    %% Test data using the trained classifier
    tic
    labelsGT{Ni} = balancedtestDataCN_4groups_30b19ndmp.Group;
    fprintf('\t Testing data .....');
    [labelsPred{Ni},scores{Ni}] = trainedClassifier.predictFcn(balancedtestDataCN_4groups_30b19ndmp);
    toc

    %% Calculate the testing accuracy of the classifier for this iteration
    pred = char(labelsPred{Ni});
    iscorrect = pred == (char(balancedtestDataCN_4groups_30b19ndmp.Group));
    iscorrect1 = iscorrect(:,2);
    accTest(Ni) = sum(iscorrect1)/(length(balancedtestDataCN_4groups_30b19ndmp.Group));

    fprintf('Accuracy in iteration number %d = %.2f \n', Ni,accTest(Ni));

    % Calculate and store confusion matrix in a 3D matrix cMat
    cMat(:,:,Ni) = confusionmat(labelsGT{Ni},labelsPred{Ni});

    % Calculate evaluation metrics using the confusion matrix
    % The evaluation metrics are calculated using one-vs-all
    % The following metrics are calculated from the confusion matrix for EACH
    % class
    for i = 1:1:length(group)
        diagConf = diag(cMat(:,:,Ni));  % Diagnoal elements
        sumRows = sum(cMat(:,:,Ni),2);
        sumColumns = sum(cMat(:,:,Ni),1);
        totalSumConf = sum(sumRows); % Same as sum(sumColumns)

        tp(Ni,i) = diagConf(i); % Diagonal element corresponding to the class
        fp(Ni,i) = sumColumns(i)-tp(Ni,i);
        fn(Ni,i) = sumRows(i)-tp(Ni,i);
        tn(Ni,i) = totalSumConf - tp(Ni,i) - fp(Ni,i) - fn(Ni,i);

        accOVAClass(Ni,i) = (tp(Ni,i)+tn(Ni,i))/(tp(Ni,i)+tn(Ni,i)+fp(Ni,i)+fn(Ni,i));
        sensClass(Ni,i) = tp(Ni,i)/(tp(Ni,i)+fn(Ni,i));
        specClass(Ni,i) = tn(Ni,i)/(tn(Ni,i)+fp(Ni,i));
        ppvClass(Ni,i) = tp(Ni,i)/(tp(Ni,i)+fp(Ni,i));
        npvClass(Ni,i) = tn(Ni,i)/(tn(Ni,i)+fn(Ni,i));
        f1scoreClass(Ni,i) = 2/((1/ppvClass(Ni,i))+(1/sensClass(Ni,i)));

        % Calculating the AUROC per class, but not saving the other
        % variables
        [xROC,yROC,tROC, aurocClass(Ni,i)] = perfcurve(labelsGT{Ni},scores{Ni}(:,i),group{i});


    end

    % Find macro-averaged values
    accOVAMacroAvg(Ni) = mean(accOVAClass(Ni,:));
    sensMacroAvg(Ni) = mean(sensClass(Ni,:));
    specMacroAvg(Ni) = mean(specClass(Ni,:));
    ppvMacroAvg(Ni) = mean(ppvClass(Ni,:));
    npvMacroAvg(Ni) = mean(npvClass(Ni,:));
    f1scoreMacroAvg(Ni) = mean(f1scoreClass(Ni,:));
    aurocMacroAvg(Ni) = mean(aurocClass(Ni,:));

    %{
    % For calculating ROC curve and AUC for each iteration (and each class)
    for iClass = 1:1:length(group)
        [xROCn{Ni}{iClass},yROCn{Ni}{iClass},tROCn{Ni}{iClass}, AUCn(Ni,iClass)] = perfcurve(labelsGT{Ni},scores{Ni}(:,iClass),group{iClass});
    end
    %}

    % Merge all labels, labelsPred{Ni} and scores in larger arrays/matrices 
    % concatenating the results from all the runs
    if Ni ==1
        labelsGTMerged = labelsGT{Ni};
        scoresMerged = scores{Ni};
        labelsPredMerged = labelsPred{Ni};
    else
        labelsGTMerged = [labelsGTMerged; labelsGT{Ni}];
        scoresMerged = [scoresMerged; scores{Ni}];
        labelsPredMerged = [labelsPredMerged; labelsPred{Ni}];
    end

    
end
%%

% Calculate mean and confidence intervals (assuming normal distribution)

pdAccVal = fitdist(accVal,'Normal');
accValMean = mean(accVal);
ciAccVal = paramci(pdAccVal);

fprintf('Validation accuracy after %d iterations = %.2f (%.2f, %.2f) \n', nIterations,accValMean,ciAccVal(1,1), ciAccVal(2,1));

pdAccTest = fitdist(accTest,'Normal');
accTestMean = mean(accTest);
ciAccTest = paramci(pdAccTest);

fprintf('Testing accuracy after %d iterations = %.2f (%.2f, %.2f) \n', nIterations,accTestMean,ciAccTest(1,1), ciAccTest(2,1));

pdauroc = fitdist(aurocMacroAvg,'Normal');
aurocMean = mean(aurocMacroAvg);
ciAUROC = paramci(pdauroc);

fprintf('Macro-averaged AUROC after %d iterations = %.2f (%.2f, %.2f) \n', nIterations,aurocMean,ciAUROC(1,1), ciAUROC(2,1));

%% Create a matrix and table summarizing the evaluation metrics

metricNames = ["AUROC", "F1-Score","Accuracy (OvA)", "Sensitivity", "Specificity", "Positive Predictive Value", "Negative Predictive Value"];
metricVariables = ["aurocClass","f1scoreClass", "accOVAClass", "sensClass", "specClass", "ppvClass", "npvClass"];
metricVariablesMacroAvg = ["aurocMacroAvg","f1scoreMacroAvg", "accOVAMacroAvg", "sensMacroAvg", "specMacroAvg", "ppvMacroAvg", "npvMacroAvg"];

columnNames = ["mean_AA","median_AA","ciL_AA", "ciU_AA", "mean_ABeta","median_ABeta","ciL_ABeta", "ciU_ABeta", "mean_AS","median_AS","ciL_AS", "ciU_AS", "mean_SCD","median_SCD","ciL_SCD", "ciU_SCD", "mean_macroAvg","median_macroAvg","ciL_macroAvg", "ciU_macroAvg"];

evalMetricsMatrix = zeros(length(metricNames), length(columnNames));

for rowI = 1:1:length(metricVariables)
    
    % Calculate evaluation metrics per group, saved according to
    % columnNames variable above
    for i = 1:1:length(group)
        % e.g. classVar = aurocClass(:,i);
        classVar = eval(strcat(metricVariables(rowI),'(:,i)'));
        evalMetricsMatrix(rowI,4*(i-1)+1) = mean(classVar);
        evalMetricsMatrix(rowI,4*(i-1)+2) = median(classVar);

        pd = fitdist(classVar,'Normal');
        ci = paramci(pd);
        evalMetricsMatrix(rowI,4*(i-1)+3) = ci(1,1);
        evalMetricsMatrix(rowI,4*(i-1)+4) = ci(2,1);

    end
    % Calculate macro-averaged evaluation metrics
    % e.g. macroAvgVar = aurocMacroAvg;
    macroAvgVar = eval(metricVariablesMacroAvg(rowI));

    evalMetricsMatrix(rowI,17) = mean(macroAvgVar);
    evalMetricsMatrix(rowI,18) = median(macroAvgVar);
    
    pd = fitdist(macroAvgVar,'Normal');
    ci = paramci(pd);
    evalMetricsMatrix(rowI,19) = ci(1,1);
    evalMetricsMatrix(rowI,20) = ci(2,1);

end

% Create table
evalMetricsSummary = array2table(evalMetricsMatrix);
evalMetricsSummary.Properties.VariableNames = columnNames;
evalMetricsSummary.Properties.RowNames = metricNames;
disp(evalMetricsSummary); % Print in command line (logging)

%% Create a matrix and table summarizing the overall model metrics

metricNames = ["Validation accuracy", "Testing accuracy", "AUROC (Macro-averaged)", "F1-Score (Macro-averaged)","Accuracy OvA (Macro-averaged)", "Sensitivity (Macro-averaged)", "Specificity (Macro-averaged)", "Positive Predictive Value (Macro-averaged)", "Negative Predictive Value (Macro-averaged)"];
metricVariables = ["accVal", "accTest", "aurocMacroAvg","f1scoreMacroAvg", "accOVAMacroAvg", "sensMacroAvg", "specMacroAvg", "ppvMacroAvg", "npvMacroAvg"];

columnNames = ["Mean","Median","CI lower", "CI upper"];

overallMetricsMatrix = zeros(length(metricNames), length(columnNames));

for rowI = 1:1:length(metricVariables)
    
    % Calculate evaluation metrics per group, saved according to
    % columnNames variable above
    modelVar = eval(metricVariables(rowI));
    overallMetricsMatrix(rowI,1) = mean(modelVar);
    overallMetricsMatrix(rowI,2) = median(modelVar);

    pd = fitdist(modelVar,'Normal');
    ci = paramci(pd);
    overallMetricsMatrix(rowI,3) = ci(1,1);
    overallMetricsMatrix(rowI,4) = ci(2,1);


end

% Create table
overallMetricsSummary = array2table(overallMetricsMatrix);
overallMetricsSummary.Properties.VariableNames = columnNames;
overallMetricsSummary.Properties.RowNames = metricNames;
disp(overallMetricsSummary); % Print in command line (logging)

%% Save variables

tic

% Create a folder with the name of the model
if not(isfolder(strcat("/",classifierName)))
    mkdir(classifierName);
end

fprintf("Saving the general model results as modelResults.mat .... \n");
save(strcat(classifierName,'/modelResults.mat'), "accTest", "accVal", "labelsPred", "labelsGT", "scores", "cMat", "nTrainPreBal", "nTrainPostBal", "nTestPreBal", "nTestPostBal",'-v7.3');

fprintf("Saving the evaluation metrics as evalMetrics.mat .... \n");
save(strcat(classifierName,'/evalMetrics.mat'), "tp", "fp", "tn", "tp", "accOVAClass", "sensClass", "specClass", "aurocClass", "ppvClass", "npvClass", "f1scoreClass", "accOVAMacroAvg", "sensMacroAvg", "specMacroAvg", "aurocMacroAvg", "ppvMacroAvg", "npvMacroAvg", "f1scoreMacroAvg","accTestMean", "ciAccTest", "aurocMean", "ciAUROC","accValMean", "ciAccVal", "evalMetricsMatrix", "evalMetricsSummary", "overallMetricsMatrix","overallMetricsSummary",'-v7.3');

% Save tables
fprintf("Saving tables as .csv files .... \n");
writetable(overallMetricsSummary, strcat(classifierName,'/overallMetricsSummary.csv'),'WriteRowNames',true);
writetable(evalMetricsSummary, strcat(classifierName,'/evalMetricsSummary.csv'), 'WriteRowNames',true);

toc


%% Save confusion matrix merged and ROC curves
tic

cMatMerged = confusionmat(labelsGTMerged,labelsPredMerged);

xROCmerged = cell(1,length(group)); yROCmerged = cell(1,length(group));
tROCmerged = cell(1,length(group));
aurocMerged = zeros(1,length(group));
for iClass = 1:1:length(group)
    [xROCmerged{iClass},yROCmerged{iClass},tROCmerged{iClass}, aurocMerged(iClass)] = perfcurve(labelsGTMerged,scoresMerged(:,iClass),group{iClass});
end
aurocMergedMacroAvg = mean(aurocMerged);
fprintf('Average AUC of MERGED ROC curves for all groups after %d iterations = %.2f%% \n', nIterations,aurocMergedMacroAvg);

% Save merged parameters
fprintf("Saving the merged model results as modelResultsMerged.mat .... \n");
save(strcat(classifierName,'/modelResultsMerged.mat'),  "labelsPredMerged", "labelsGTMerged", "scoresMerged", "cMatMerged", "xROCmerged", "yROCmerged", "tROCmerged", "aurocMerged", "aurocMergedMacroAvg",'-v7.3'); 
toc
%% Plot merged confusion matrix
tic

if not(isfolder("/Figures"))
    mkdir("Figures");
end

fprintf("Plotting confusion matrices .... \n");

figure;
confusionchart(labelsGTMerged,labelsPredMerged,'RowSummary','row-normalized','ColumnSummary','column-normalized');
%set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
set(findall(gcf,'-property','FontSize'),'FontSize',8)
fnamePNG = strcat('Figures/1_','confusionMatrixMerged_Coverslip_',classifierName,'.png');
fnameFig = strcat('Figures/1_','confusionMatrixMerged_Coverslip_',classifierName,'.fig');
print(gcf,fnamePNG,'-dpng','-r600');
saveas(gcf,fnameFig);
close;

figure;
confusionchart(labelsGTMerged,labelsPredMerged, 'Normalization','row-normalized');
%set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
%set(findall(gcf,'-property','FontSize'),'FontSize',8)
fnamePNG = strcat('Figures/2_','confusionMatrixMerged_rowNormalized_Coverslip_',classifierName,'.png');
fnameFig = strcat('Figures/2_','confusionMatrixMerged__rowNormalized_Coverslip_',classifierName,'.fig');
print(gcf,fnamePNG,'-dpng','-r600');
saveas(gcf,fnameFig);
close;
toc
% Can save both versions

%% Plot merged ROC curve
tic
fprintf("Plotting ROC curves .... \n");
figure;
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [231/255 41/255 138/255], [166/255 118/255 29/255],[35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2

hold on;
% Plot all ROC curves of different classes
for i = 1:1:length(group)
    plot(xROCmerged{i}, yROCmerged{i}, 'color', colors{i})
end

% Plot macro-averaged ROC curve
rocObj = rocmetrics(labelsGTMerged,scoresMerged,group);
[xROCavg,yROCavg,tROCavg,aurocAvg] = average(rocObj,"macro");
plot(xROCavg,yROCavg, 'color',[0.5 0.5 0.5], 'LineWidth',1.5);

plot([0,1],[0,1],"k--");
hold off;
xlabel("1 - Specificity");
ylabel("Sensitivity");
box on;
%axis padded;
legend(join([group(1) ," (AUROC = ", sprintf('%.2f',aurocMerged(1)), ")"]), join([group(2) ," (AUROC = ", sprintf('%.2f',aurocMerged(2)), ")"]), join([group(3) ," (AUROC = ", sprintf('%.2f',aurocMerged(3)), ")"]), join([group(4) ," (AUROC = ", sprintf('%.2f',aurocMerged(4)), ")"]), join(["Macro-averaged (AUROC = ", sprintf('%.2f',aurocMerged(1)), ")"]), 'Location', 'southeast'); % Legends with AUC curve
%set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
%set(findall(gcf,'-property','FontSize'),'FontSize',8)
fnamePNG = strcat('Figures/3_','ROCcurveMerged_Coverslip_',classifierName,'.png');
fnameFig = strcat('Figures/3_','ROCcurveMerged_Coverslip_',classifierName,'.fig');
print(gcf,fnamePNG,'-dpng','-r600');
saveas(gcf,fnameFig);
close;

%{
% Automatically curve ROC curve based on rocmetrics
rocObj = rocmetrics(labelsGTMerged,scoresMerged,group);
figure; 
plot(rocObj,AverageROCType="macro");
%}


toc




% carpenter30b19ndmpTrainTest3groupsCoverslipL_singleSplit.m
% Program to read normalized frequency distribution data (combined data for
% Nepal and Canada), select only non-dimensional morphological parameters,
% create 80:20 training:testing split, select features based on statistical
% analysis (done previously) and balance classes

% Note: The training/testing data split is done participant or donor-wise
% ps UBC 2023

% Full list of morphological parameters: 
% Area,Perimeter,Angle,Major,Minor,Height,Width,Mean,StdDev,Median,Skewness,Kurtosis,Feret,FeretAngle,MinFeret,FeretX,FeretY,Circularity,AR,Roundness,ConvexHullArea,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,ConvexPerimeter,Convexity,FibreLength,FibreWidth,Curl,CurlW,AreaNorm,PerimeterNorm,MajorNorm,MinorNorm,HeightNorm,WidthNorm,FeretNorm,MinFeretNorm

% List of 19 non-dimensional morphological parameters: 
% Circularity,AR,Roundness,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,Convexity,Curl,CurlW,AreaNorm,PerimeterNorm,MajorNorm,MinorNorm,HeightNorm,WidthNorm,FeretNorm,MinFeretNorm

%% 1) Import histData for processing of morphological parameters
% Note: Here, the histDataCS (or coverslip-wise) histcount data is selected
% (For finer analysis, image-wise data or groups of image-wise data can be
% selected, for coarser analysis, donor-wise data can be selected

% Load all histData - [Note: Need to modify path based on where the data is
% stored in local machine or cluster]
tic 
fprintf('Reading data - generalInfo.mat and histData.mat \n');

group = {'AA', 'ABeta','AS','SCD'}; %Group name for reading data
% Actual group name = {'AA-ABeta','AS','SCD'};

for i = 1:1:length(group)
    % Run similar command for all groups to read generalInfo.mat and
    % histData.mat
    % Example: AA.info = load('../../../1DCharacterization/t02hrs/cnCombined/30b40mp/AA/generalInfo.mat');
    % Example: AA.hD = load('../../../1DCharacterization/t02hrs/cnCombined/30b40mp/AA/histData.mat');
    eval(strcat(group{i},".info = load('../../../../Plots/1DCharacterization/t02hrs/cnCombined/30b40mp/",group{i},"/generalInfo.mat');"));
    eval(strcat(group{i},".hD = load('../../../../Plots/1DCharacterization/t02hrs/cnCombined/30b40mp/",group{i},"/histData.mat');"));

end

toc


%% 2) Randomly separate data into training and testing sets
% Training:testing data split approximately 80:20 
% It is important that the split is participant or donor wise, such that
% all the data (e.g. coverslips) from one participant only appears in
% either trainig or testing data

% The numbers of all groups are "Group (Total) = Training:Testing"
% The following example is for a trainRatio of 0.8 (train:test = 80:20)
% AA (30) + ABeta (23) = 42:11
% AS (45) = 36:9
% SCD (40) = 32:8

tic

% Training dataset ratio
trainRatio = 0.8;
fprintf('Randomly separating data, with training:testing = %.f:%.f \n',trainRatio*100, (1-trainRatio)*100);

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

toc


%% 3) Store training data into table
% Load histData for each morphological parameter as a feature (in the form
% of an array) and the group information related to that cover-slip or
% donor. If each morphological parameter has 'nbins' bins, then each
% parameter occupies 'nbins' rows, e.g. for Area: Area1, Area2, ....
tic

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

group = {'AA-ABeta','AS','SCD'}; %Group names used for analysis

fprintf('Storing training data into table \n');
% Fill in the first row with the information from AA (first coverslip only)
% The morphological parameters and the groupName are added as columns
trainDataCN_3groups_30b19ndmp = cell2table([trainAA{1}{1}, char(group{1})]);

% Add the remaining rows
% Loop through all donors or UIDs in the group AA
fprintf('\t Storing data from AA \n');
for uidI = 1:1:length(trainAA)
    % Loop through all coverslips of the donor
    for csI = 1:1:length(trainAA{1,uidI})
        if tableRow > 1
            trainDataCN_3groups_30b19ndmp = [trainDataCN_3groups_30b19ndmp; cell2table([trainAA{1,uidI}{1,csI}, char(group{1})])];
        end
        tableRow = tableRow + 1;
    end
end

% Loop through all donors or UIDs in the group ABeta
fprintf('\t Storing data from ABeta \n');
for uidI = 1:1:length(trainABeta)
    % Loop through all coverslips of the donor
    for csI = 1:1:length(trainABeta{1,uidI})
        trainDataCN_3groups_30b19ndmp = [trainDataCN_3groups_30b19ndmp; cell2table([trainABeta{1,uidI}{1,csI}, char(group{1})])];
        tableRow = tableRow + 1;
    end
end

% Loop through all donors or UIDs in the group AS
fprintf('\t Storing data from AS \n');
for uidI = 1:1:length(trainAS)
    % Loop through all coverslips of the donor
    for csI = 1:1:length(trainAS{1,uidI})
        trainDataCN_3groups_30b19ndmp = [trainDataCN_3groups_30b19ndmp; cell2table([trainAS{1,uidI}{1,csI}, char(AS.info.groupName)])];
        tableRow = tableRow + 1;
    end
end

% Loop through all donors or UIDs in the group SS
fprintf('\t Storing data from SCD (SBeta and SS) \n');
for uidI = 1:1:length(trainSCD)
    % Loop through all coverslips of the donor
    for csI = 1:1:length(trainSCD{1,uidI})
        trainDataCN_3groups_30b19ndmp = [trainDataCN_3groups_30b19ndmp; cell2table([trainSCD{1,uidI}{1,csI}, char(SCD.info.groupName)])];
        tableRow = tableRow + 1;
    end
end

% After all additions, make tableRow equal to number of rows
tableRow = tableRow -1;

% Replace variable names with morphParamNameC 
trainDataCN_3groups_30b19ndmp.Properties.VariableNames = columnNames;
% Make Group variable categorical
trainDataCN_3groups_30b19ndmp.Group = categorical(trainDataCN_3groups_30b19ndmp.Group);
% Add description
trainDataCN_3groups_30b19ndmp.Properties.Description = 'Table with histDataCS30b40mp for coverslips from training set and 11 non-dimensional morphological parameters';

% Delete the parameters that are dimensional or intensity based
for ii = 1:1:length(delete21)
    eval(strcat('trainDataCN_3groups_30b19ndmp.',morphParamNameC{delete21(ii)},'=[];'));
end
toc

%% 4) Store testing data into table
% Load histData for each morphological parameter as a feature (in the form
% of an array) and the group information related to that cover-slip or
% donor. If each morphological parameter has 'nbins' bins, then each
% parameter occupies 'nbins' rows, e.g. for Area: Area1, Area2, ....
tic

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

fprintf('Storing testing data into table \n');
% Fill in the first row with the information from AA (first coverslip only)
% The morphological parameters and the groupName are added as columns
testDataCN_3groups_30b19ndmp = cell2table([testAA{1}{1}, char(group{1})]);

% Add the remaining rows
% Loop through all donors or UIDs in the group AA
fprintf('\t Storing data from AA \n');
for uidI = 1:1:length(testAA)
    % Loop through all coverslips of the donor
    for csI = 1:1:length(testAA{1,uidI})
        if tableRowTest > 1
            testDataCN_3groups_30b19ndmp = [testDataCN_3groups_30b19ndmp; cell2table([testAA{1,uidI}{1,csI}, char(group{1})])];
        end
        tableRowTest = tableRowTest + 1;
    end
end

% Loop through all donors or UIDs in the group ABeta
fprintf('\t Storing data from ABeta \n');
for uidI = 1:1:length(testABeta)
    % Loop through all coverslips of the donor
    for csI = 1:1:length(testABeta{1,uidI})
        testDataCN_3groups_30b19ndmp = [testDataCN_3groups_30b19ndmp; cell2table([testABeta{1,uidI}{1,csI}, char(group{1})])];
        tableRowTest = tableRowTest + 1;
    end
end

% Loop through all donors or UIDs in the group AS
fprintf('\t Storing data from AS \n');
for uidI = 1:1:length(testAS)
    % Loop through all coverslips of the donor
    for csI = 1:1:length(testAS{1,uidI})
        testDataCN_3groups_30b19ndmp = [testDataCN_3groups_30b19ndmp; cell2table([testAS{1,uidI}{1,csI}, char(AS.info.groupName)])];
        tableRowTest = tableRowTest + 1;
    end
end

% Loop through all donors or UIDs in the group SS
fprintf('\t Storing data from SCD (SBeta and SS) \n');
for uidI = 1:1:length(testSCD)
    % Loop through all coverslips of the donor
    for csI = 1:1:length(testSCD{1,uidI})
        testDataCN_3groups_30b19ndmp = [testDataCN_3groups_30b19ndmp; cell2table([testSCD{1,uidI}{1,csI}, char(SCD.info.groupName)])];
        tableRowTest = tableRowTest + 1;
    end
end

% After all additions, make tableRow equal to number of rows
tableRowTest = tableRowTest -1;

% Replace variable names with morphParamNameC 
testDataCN_3groups_30b19ndmp.Properties.VariableNames = columnNames;
% Make Group variable categorical
testDataCN_3groups_30b19ndmp.Group = categorical(testDataCN_3groups_30b19ndmp.Group);
% Add description
testDataCN_3groups_30b19ndmp.Properties.Description = 'Table with histDataCS30b40mp for coverslips from testing set and 11 non-dimensional morphological parameters';

% Delete the parameters that are dimensional or intensity based
for ii = 1:1:length(delete21)
    eval(strcat('testDataCN_3groups_30b19ndmp.',morphParamNameC{delete21(ii)},'=[];'));
end
toc

%% 5) Save tables with training and testing data
%{
tic
roundN = 'round3';

fprintf("Saving table trainDataCN_3groups_30b19ndmp.... \n");
save(strcat('randomSplits3groups/',roundN,'/trainDataCN_3groups_30b19ndmp.mat'), 'trainDataCN_3groups_30b19ndmp');

fprintf("Saving table testDataCN_3groups_30b19ndmp.... \n");
save(strcat('randomSplits3groups/',roundN,'/testDataCN_3groups_30b19ndmp.mat'), 'testDataCN_3groups_30b19ndmp');

toc
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          BALANCE CLASSES                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 5) Find number of all classes

trainG = cellstr(trainDataCN_3groups_30b19ndmp.Group);
testG = cellstr(testDataCN_3groups_30b19ndmp.Group);

AA_ABetaTrainInd = find(contains(trainG,'AA-ABeta'));
AA_ABetaTestInd = find(contains(testG,'AA-ABeta'));

ASTrainInd = find(contains(trainG,'AS'));
ASTestInd = find(contains(testG,'AS'));

SCDTrainInd = find(contains(trainG,'SCD'));
SCDTestInd = find(contains(testG,'SCD'));

fprintf('Training data numbers ..... \n');
fprintf('\t Number of AA-ABeta in training data = %d \n', length(AA_ABetaTrainInd));
fprintf('\t Number of AS in training data = %d \n', length(ASTrainInd));
fprintf('\t Number of SCD in training data = %d \n', length(SCDTrainInd));

fprintf('Testing data numbers ..... \n');
fprintf('\t Number of AA-ABeta in testing data = %d \n', length(AA_ABetaTestInd));
fprintf('\t Number of AS in testing data = %d \n', length(ASTestInd));
fprintf('\t Number of SCD in testing data = %d \n', length(SCDTestInd));

nTrain = length(AA_ABetaTrainInd)+length(ASTrainInd)+length(SCDTrainInd);
nTest = length(AA_ABetaTestInd) + length(ASTestInd) + length(SCDTestInd);
fprintf('Training data percentages ..... \n');
fprintf('\t Number of AA-ABeta in training data = %.1f%% \n', length(AA_ABetaTrainInd)*100/nTrain);
fprintf('\t Number of AS in training data = %.1f%% \n', length(ASTrainInd)*100/nTrain);
fprintf('\t Number of SCD in training data = %.1f%% \n', length(SCDTrainInd)*100/nTrain);

fprintf('Testing data percentages ..... \n');
fprintf('\t Number of AA-ABeta in testing data = %.1f%% \n', length(AA_ABetaTestInd)*100/nTest);
fprintf('\t Number of AS in testing data = %.1f%% \n', length(ASTestInd)*100/nTest);
fprintf('\t Number of SCD in testing data = %.1f%% \n', length(SCDTestInd)*100/nTest);


%% 6) Balance the groups by upsampling the minority data

balancedtrainDataCN_3groups_30b19ndmp = trainDataCN_3groups_30b19ndmp;
balancedtestDataCN_3groups_30b19ndmp = testDataCN_3groups_30b19ndmp;

nMaxTrain = length(AA_ABetaTestInd);
nMaxTest = length(AA_ABetaTestInd);


% AS
% Balance training set
% Randomly select and copy a row for AS at the end of the table
fprintf('Upsampling AS class by randomly repeating %d rows from the training set \n', nMaxTrain - length(ASTrainInd));
for ii = 1:1:(nMaxTrain - length(ASTrainInd))
    % Randomly choose an index from the AS list of indices
    randInd = randsample(ASTrainInd,1);

    % New table (with one row) corresponding to the randomly sampled index
    t1 = trainDataCN_3groups_30b19ndmp(randInd,:);

    % Concatenate new table 
    balancedtrainDataCN_3groups_30b19ndmp = [balancedtrainDataCN_3groups_30b19ndmp; t1];

end

% Balance testing set
% Randomly select and copy a row for AS at the end of the table
fprintf('Upsampling AS class by randomly repeating %d rows from the testing set \n', nMaxTest - length(ASTestInd));
for ii = 1:1:(nMaxTest - length(ASTestInd))
    % Randomly choose an index from the AS list of indices
    randInd = randsample(ASTestInd,1);

    % New table (with one row) corresponding to the randomly sampled index
    t1 = testDataCN_3groups_30b19ndmp(randInd,:);

    % Concatenate new table 
    balancedtestDataCN_3groups_30b19ndmp = [balancedtestDataCN_3groups_30b19ndmp; t1];

end


% SCD
% Balance training set
% Randomly select and copy a row for SCD at the end of the table
fprintf('Upsampling SCD class by randomly repeating %d rows from the training set \n', nMaxTrain - length(SCDTrainInd));
for ii = 1:1:(nMaxTrain - length(SCDTrainInd))
    % Randomly choose an index from the SCD list of indices
    randInd = randsample(SCDTrainInd,1);

    % New table (with one row) corresponding to the randomly sampled index
    t1 = trainDataCN_3groups_30b19ndmp(randInd,:);

    % Concatenate new table 
    balancedtrainDataCN_3groups_30b19ndmp = [balancedtrainDataCN_3groups_30b19ndmp; t1];

end

% Balance testing set
% Randomly select and copy a row for SCD at the end of the table
fprintf('Upsampling SCD class by randomly repeating %d rows from the testing set \n', nMaxTest - length(SCDTestInd));
for ii = 1:1:(nMaxTest - length(SCDTestInd))
    % Randomly choose an index from the AS list of indices
    randInd = randsample(SCDTestInd,1);

    % New table (with one row) corresponding to the randomly sampled index
    t1 = testDataCN_3groups_30b19ndmp(randInd,:);

    % Concatenate new table 
    balancedtestDataCN_3groups_30b19ndmp = [balancedtestDataCN_3groups_30b19ndmp; t1];

end



%% 7) Check the number of SS and SBeta in the balanced set

trainG1 = cellstr(balancedtrainDataCN_3groups_30b19ndmp.Group);
testG1 = cellstr(balancedtestDataCN_3groups_30b19ndmp.Group);

AA_ABetaTrainInd1 = find(contains(trainG1,'AA-ABeta'));
AA_ABetaTestInd1 = find(contains(testG1,'AA-ABeta'));

ASTrainInd1 = find(contains(trainG1,'AS'));
ASTestInd1 = find(contains(testG1,'AS'));

SCDTrainInd1 = find(contains(trainG1,'SCD'));
SCDTestInd1 = find(contains(testG1,'SCD'));

fprintf('Training data numbers after balancing classes ..... \n');
fprintf('\t Number of AA-ABeta in training data = %d \n', length(AA_ABetaTrainInd1));
fprintf('\t Number of AS in training data = %d \n', length(ASTrainInd1));
fprintf('\t Number of SCD in training data = %d \n', length(SCDTrainInd1));

fprintf('Testing data numbers after balancing classes ..... \n');
fprintf('\t Number of AA-ABeta in testing data = %d \n', length(AA_ABetaTestInd1));
fprintf('\t Number of AS in testing data = %d \n', length(ASTestInd1));
fprintf('\t Number of SCD in testing data = %d \n', length(SCDTestInd1));

%% 8) Save the balanced classes 
tic

fprintf("Saving table balancedtrainDataCN_3groups_30b19ndmp.... \n");
save('balancedtrainDataCN_3groups_30b19ndmp.mat', 'balancedtrainDataCN_3groups_30b19ndmp');

fprintf("Saving table balancedtestDataCN_3groups_30b19ndmp.... \n");
save('balancedtestDataCN_3groups_30b19ndmp.mat', 'balancedtestDataCN_3groups_30b19ndmp');

toc

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          FEATURE SELECTION                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 9) Import indices of minimum p values from statistical analysis
tic 
fprintf('Reading pValuesEtAl.mat \n');
load('../../../../Plots/statisticalDifferences/30b40mp/coverslipLevel/pValuesEtAl.mat');
toc

%% 10) Create table with selected features

tic
% Iterate through all non-dimensional morphological parameters
for k = 1:1:length(pMinUniqueNDMP)

    % Iterate through all indices in the morphological parameter
    for kP = 1:1:length(pMinUniqueNDMP{k})
        % For the first value, create a table, for the rest concatenate
        if k ==1 && kP ==1
            featSelBaltestDataCN_3groups_30b19ndmp = eval(strcat('array2table(balancedtestDataCN_3groups_30b19ndmp.',morphParamName{keep19(k)},'(:,',num2str(pMinUniqueNDMP{k}(kP)),'))'));
            featSelBaltrainDataCN_3groups_30b19ndmp = eval(strcat('array2table(balancedtrainDataCN_3groups_30b19ndmp.',morphParamName{keep19(k)},'(:,',num2str(pMinUniqueNDMP{k}(kP)),'))'));

            varName = strcat(morphParamName{keep19(k)},num2str(pMinUniqueNDMP{k}(kP)));
            featSelBaltestDataCN_3groups_30b19ndmp = renamevars(featSelBaltestDataCN_3groups_30b19ndmp,"Var1",varName);
            featSelBaltrainDataCN_3groups_30b19ndmp = renamevars(featSelBaltrainDataCN_3groups_30b19ndmp,"Var1",varName);
        else
            featSelBaltestDataCN_3groups_30b19ndmp = [featSelBaltestDataCN_3groups_30b19ndmp,eval(strcat('array2table(balancedtestDataCN_3groups_30b19ndmp.',morphParamName{keep19(k)},'(:,',num2str(pMinUniqueNDMP{k}(kP)),'))'))];
            featSelBaltrainDataCN_3groups_30b19ndmp = [featSelBaltrainDataCN_3groups_30b19ndmp,eval(strcat('array2table(balancedtrainDataCN_3groups_30b19ndmp.',morphParamName{keep19(k)},'(:,',num2str(pMinUniqueNDMP{k}(kP)),'))'))];

            varName = strcat(morphParamName{keep19(k)},num2str(pMinUniqueNDMP{k}(kP)));
            featSelBaltestDataCN_3groups_30b19ndmp = renamevars(featSelBaltestDataCN_3groups_30b19ndmp,"Var1",varName);
            featSelBaltrainDataCN_3groups_30b19ndmp = renamevars(featSelBaltrainDataCN_3groups_30b19ndmp,"Var1",varName);
        end

    end
end
% Concatenate the group information
% Make Group variable categorical
testDataCN_3groups_30b19ndmp.Group = categorical(testDataCN_3groups_30b19ndmp.Group);

featSelBaltestDataCN_3groups_30b19ndmp.Group = categorical(balancedtestDataCN_3groups_30b19ndmp.Group);
featSelBaltrainDataCN_3groups_30b19ndmp.Group = categorical(balancedtrainDataCN_3groups_30b19ndmp.Group);

% Can generate column names based on morphological parameter
toc

%% 11) Save the balanced classes 
tic

fprintf("Saving table featSelBaltestDataCN_3groups_30b19ndmp.... \n");
save('featSelBaltestDataCN_3groups_30b19ndmp.mat', 'featSelBaltestDataCN_3groups_30b19ndmp');

fprintf("Saving table balancedtestDataCN_3groups_30b19ndmp.... \n");
save('featSelBaltrainDataCN_3groups_30b19ndmp.mat', 'featSelBaltrainDataCN_3groups_30b19ndmp');

toc
%}
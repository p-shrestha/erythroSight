% morphParamCombine.m
% Program to combine morphological parameters from Nepal and Canada data
% into 4 groups - AA, ABeta, AS, SCD
% ps UBC 2023

%% 0) Different groups 

% I) AA (Nepal) + N (Canada) = AA
% ABeta (Nepal) = ABeta (no combining required)
% II) AS (Nepal) + SCT (Canada) = AS
% III) SBeta (Nepal) + SS (Nepal) + SCD (Canada) = SCD

% -------------------------------------------------------------------------
%             I of III)  AA (Nepal) + N (Canada) = AA
% -------------------------------------------------------------------------
%% 1) Load data - morphParamData and generalInfo

% morphParam is a 2D cell with morphological parameters for all
% the data. The two levels are 1) donor or UID, and 2)
% morphological parameter or k

tic
fprintf('Reading data - morphParamData.mat and generalInfo.mat for AA and N \n');

groupNameR1 = 'AA'; % Group name to read 
groupNameR2 = 'N'; % Group name to read 

R1.info = load(strcat("../cnSeparated/",groupNameR1,"/generalInfo.mat"));
R1.mP = load(strcat("../cnSeparated/",groupNameR1,"/morphParamData.mat"));
R2.info = load(strcat("../cnSeparated/",groupNameR2,"/generalInfo.mat"));
R2.mP = load(strcat("../cnSeparated/",groupNameR2,"/morphParamData.mat"));

toc


%% 2) Combine morphParam and general information data

tic
fprintf('Combining morphological parameter data and general information for Canada and Nepal data \n');
groupNameW = 'AA'; % Group name to write

morphParam = cat(2, R1.mP.morphParam, R2.mP.morphParam); % Concatenate data from Nepal and Canada

% Print size before and after concatenation
fprintf('Length of group 1 = %d before concatenation.\n',length(R1.mP.morphParam));
fprintf('Length of group 2 = %d before concatenation.\n',length(R2.mP.morphParam));
fprintf('Length of combined group = %d AFTER concatenation.\n',length(morphParam));
fprintf('Difference between lengths before and after concatenation = %d.\n',(length(R1.mP.morphParam)+length(R2.mP.morphParam)-length(morphParam)));

% Combine variables for general info
groupName = groupNameW;
tName = R1.info.tName;
morphParamName = R1.info.morphParamName;
morphParamAxis = R1.info.morphParamAxis;
nImagesGroup = R1.info.nImagesGroup + R2.info.nImagesGroup;
nCellsGroup = R1.info.nCellsGroup + R2.info.nCellsGroup;

nCells = cat(2,R1.info.nCells, R2.info.nCells); 
nImages = cat(2,R1.info.nImages, R2.info.nImages);

fprintf('Number of donors in nCells = %d.\n', length(nCells));
fprintf('Difference between number of donors in nCells before and after concatenation = %d.\n \n',(length(R1.info.nCells)+length(R2.info.nCells)-length(nCells)));

fprintf('Number of donors in nImages = %d.\n', length(nImages));
fprintf('Difference between number of donors in nImages before and after concatenation = %d.\n',(length(R1.info.nImages)+length(R2.info.nImages)-length(nImages)));

toc


%% 3) Save variables for morphParam and general information
% Create a folder with the groupName if one does not exist
if not(isfolder(groupName))
    mkdir(groupName)
end
 
tic 
fprintf("Saving general information variables.... \n");
save(strcat(groupName,"/generalInfo.mat"), "groupName", "tName", "morphParamName", "morphParamAxis", "nImages","nCells", "nImagesGroup", "nCellsGroup", '-v7.3');
 
fprintf("Saving variables for morphParam.... \n");
save(strcat(groupName, "/morphParamData.mat"),"morphParam", '-v7.3');

toc


% -------------------------------------------------------------------------
%             II of III)  AS (Nepal) + SCT (Canada) = AS
% -------------------------------------------------------------------------
%% 4) Load data - morphParamData and generalInfo

% morphParam is a 2D cell with morphological parameters for all
% the data. The two levels are 1) donor or UID, and 2)
% morphological parameter or k

tic
% Clear all variables for next iteration
fprintf('Clearing all variables for AA and N \n');
clear;

fprintf('Reading data - morphParamData.mat and generalInfo.mat for AS and SCT \n');

groupNameR1 = 'AS'; % Group name to read 
groupNameR2 = 'SCT'; % Group name to read 

R1.info = load(strcat("../cnSeparated/",groupNameR1,"/generalInfo.mat"));
R1.mP = load(strcat("../cnSeparated/",groupNameR1,"/morphParamData.mat"));
R2.info = load(strcat("../cnSeparated/",groupNameR2,"/generalInfo.mat"));
R2.mP = load(strcat("../cnSeparated/",groupNameR2,"/morphParamData.mat"));

toc


%% 5) Combine morphParam and general information data

tic
fprintf('Combining morphological parameter data and general information for Canada and Nepal data \n');
groupNameW = 'AS'; % Group name to write

morphParam = cat(2, R1.mP.morphParam, R2.mP.morphParam); % Concatenate data from Nepal and Canada

% Print size before and after concatenation
fprintf('Length of group 1 = %d before concatenation.\n',length(R1.mP.morphParam));
fprintf('Length of group 2 = %d before concatenation.\n',length(R2.mP.morphParam));
fprintf('Length of combined group = %d AFTER concatenation.\n',length(morphParam));
fprintf('Difference between lengths before and after concatenation = %d.\n',(length(R1.mP.morphParam)+length(R2.mP.morphParam)-length(morphParam)));

% Combine variables for general info
groupName = groupNameW;
tName = R1.info.tName;
morphParamName = R1.info.morphParamName;
morphParamAxis = R1.info.morphParamAxis;
nImagesGroup = R1.info.nImagesGroup + R2.info.nImagesGroup;
nCellsGroup = R1.info.nCellsGroup + R2.info.nCellsGroup;

nCells = cat(2,R1.info.nCells, R2.info.nCells); 
nImages = cat(2,R1.info.nImages, R2.info.nImages);

fprintf('Number of donors in nCells = %d.\n', length(nCells));
fprintf('Difference between number of donors in nCells before and after concatenation = %d.\n \n',(length(R1.info.nCells)+length(R2.info.nCells)-length(nCells)));

fprintf('Number of donors in nImages = %d.\n', length(nImages));
fprintf('Difference between number of donors in nImages before and after concatenation = %d.\n',(length(R1.info.nImages)+length(R2.info.nImages)-length(nImages)));

toc


%% 6) Save variables for morphParam and general information
% Create a folder with the groupName if one does not exist
if not(isfolder(groupName))
    mkdir(groupName)
end
 
tic 
fprintf("Saving general information variables.... \n");
save(strcat(groupName,"/generalInfo.mat"), "groupName", "tName", "morphParamName", "morphParamAxis", "nImages","nCells", "nImagesGroup", "nCellsGroup", '-v7.3');
 
fprintf("Saving variables for morphParam.... \n");
save(strcat(groupName, "/morphParamData.mat"),"morphParam", '-v7.3');

toc

% -------------------------------------------------------------------------
%     III of III)  SBeta (Nepal) + SS (Nepal) + SCD (Canada) = SCD
% -------------------------------------------------------------------------
%% 7) Load data - morphParamData and generalInfo

% morphParam is a 2D cell with morphological parameters for all
% the data. The two levels are 1) donor or UID, and 2)
% morphological parameter or k

tic
% Clear all variables for next iteration
fprintf('Clearing all variables for AS and SCT \n');
clear;

fprintf('Reading data - morphParamData.mat and generalInfo.mat for SBeta, SS, and SCD \n');

groupNameR1 = 'SBeta'; % Group name to read 
groupNameR2 = 'SS'; % Group name to read 
groupNameR3 = 'SCD'; % Group name to read 

R1.info = load(strcat("../cnSeparated/",groupNameR1,"/generalInfo.mat"));
R1.mP = load(strcat("../cnSeparated/",groupNameR1,"/morphParamData.mat"));
R2.info = load(strcat("../cnSeparated/",groupNameR2,"/generalInfo.mat"));
R2.mP = load(strcat("../cnSeparated/",groupNameR2,"/morphParamData.mat"));
R3.info = load(strcat("../cnSeparated/",groupNameR3,"/generalInfo.mat"));
R3.mP = load(strcat("../cnSeparated/",groupNameR3,"/morphParamData.mat"));

toc


%% 8) Combine morphParam and general information data

tic
fprintf('Combining morphological parameter data and general information for Canada and Nepal data \n');
groupNameW = 'SCD'; % Group name to write

morphParam = cat(2, R1.mP.morphParam, R2.mP.morphParam, R3.mP.morphParam); % Concatenate data from Nepal and Canada

% Print size before and after concatenation
fprintf('Length of group 1 = %d before concatenation.\n',length(R1.mP.morphParam));
fprintf('Length of group 2 = %d before concatenation.\n',length(R2.mP.morphParam));
fprintf('Length of group 3 = %d before concatenation.\n',length(R3.mP.morphParam));
fprintf('Length of combined group = %d AFTER concatenation.\n',length(morphParam));
fprintf('Difference between lengths before and after concatenation = %d.\n',(length(R1.mP.morphParam)+length(R2.mP.morphParam)+length(R3.mP.morphParam)-length(morphParam)));

% Combine variables for general info
groupName = groupNameW;
tName = R1.info.tName;
morphParamName = R1.info.morphParamName;
morphParamAxis = R1.info.morphParamAxis;
nImagesGroup = R1.info.nImagesGroup + R2.info.nImagesGroup + R3.info.nImagesGroup;
nCellsGroup = R1.info.nCellsGroup + R2.info.nCellsGroup + R3.info.nCellsGroup;

nCells = cat(2,R1.info.nCells, R2.info.nCells, R3.info.nCells); 
nImages = cat(2,R1.info.nImages, R2.info.nImages, R3.info.nImages);

fprintf('Number of donors in nCells = %d.\n', length(nCells));
fprintf('Difference between number of donors in nCells before and after concatenation = %d.\n \n',(length(R1.info.nCells)+length(R2.info.nCells)+length(R3.info.nCells)-length(nCells)));

fprintf('Number of donors in nImages = %d.\n', length(nImages));
fprintf('Difference between number of donors in nImages before and after concatenation = %d.\n',(length(R1.info.nImages)+length(R2.info.nImages)+length(R3.info.nImages)-length(nImages)));

toc


%% 9) Save variables for morphParam and general information
% Create a folder with the groupName if one does not exist
if not(isfolder(groupName))
    mkdir(groupName)
end
 
tic 
fprintf("Saving general information variables.... \n");
save(strcat(groupName,"/generalInfo.mat"), "groupName", "tName", "morphParamName", "morphParamAxis", "nImages","nCells", "nImagesGroup", "nCellsGroup", '-v7.3');
 
fprintf("Saving variables for morphParam.... \n");
save(strcat(groupName, "/morphParamData.mat"),"morphParam", '-v7.3');

toc

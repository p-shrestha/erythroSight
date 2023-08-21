% rmSBeta003Bi.m
% Program to remove SBeta_003_Bi because the coverslip was wrongly labelled
% and was likley a coverslip for a normal participant

groupName = 'SCD'; % Group name to read; The first of SCD is SBeta


load(strcat("../",groupName,"/generalInfo.mat"));
load(strcat("../",groupName,"/morphParamData.mat"));


%% Fix morphParam


tic

% Iterate through all files within groupName
fprintf('Fixing SBeta-003 for %s.\n',groupName);
uidI = 3; 

% Iterate through all coverslips in the UID or donor
for csI = 1:1:length(nCells{1,uidI})
    % cellIndexI and cellIndexF are the initial and final indices for
    % the cells for the entire coverslip
    if csI ==1
        cellIndexI = 1;
        cellIndexF = sum(nCells{1,uidI}{1,csI});
    else
        cellIndexI = cellIndexI + sum(nCells{1,uidI}{1,(csI-1)});
        cellIndexF = cellIndexI + sum(nCells{1,uidI}{1,csI}) - 1;
    end

    %Iterate through all morphological parameters (columns) within each file
    for k = 1:1:length(morphParamName)
        if csI == 3 && length(nCells{1,uidI}) == 6
            fprintf('For UID  = %d , deleting this coverslip, with coverslip ID = %d \n', uidI, csI);
            morphParam{1,uidI}{1,k}(cellIndexI:cellIndexF) = [];
        else
            fprintf('For UID  = %d , NOT deleting this coverslip, with coverslip ID = %d \n', uidI, csI);
        end
    end
    

end
    



%% Fix generalInfo

nCellsGroup = nCellsGroup - sum(nCells{1,3}{3});
nImagesGroup = nImagesGroup - nImages{1,3}(3);
nCells{1,3}(3) = [];
nImages{1,3}(3) = [];


%% 9) Save variables for morphParam and general information
% Create a folder with the groupName if one does not exist

tic 
fprintf("Saving general information variables.... \n");
save(strcat("../",groupName,"/generalInfo.mat"), "groupName", "tName", "morphParamName", "morphParamAxis", "nImages","nCells", "nImagesGroup", "nCellsGroup", '-v7.3');
 
fprintf("Saving variables for morphParam.... \n");
save(strcat("../",groupName, "/morphParamData.mat"),"morphParam", '-v7.3');

toc
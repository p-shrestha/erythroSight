% normalizeBasicParam.m
% Program to normalize some basic morphological parameters by their average
% values per image. Selected basic morphological parameters: 1. Area, 2.
% Perimeter, 4. Major, 5. Minor, 6. Height, 7. Width, 13. Feret, 15.
% MinFeret (8 parameters selected) 
% ps UBC 2023

%% 1) Load data - morphParamData and generalInfo
% groupName = 'AA', 'ABeta', 'AS', 'SCD'; 
% The groupName variable is inputted in the shell script while running the
% MATLAB program in the cluster

% morphParam is a 2D cell with morphological parameters for all
% the data. The two levels are 1) donor or UID, and 2)
% morphological parameter or k

tic
fprintf('Reading data - morphParamData.mat and generalInfo.mat for groupname %s \n', groupName);

info = load(strcat("../",groupName,"/generalInfo.mat"));
mP = load(strcat("../",groupName,"/morphParamData.mat"));

% Load variable names as local variables (for saving later)
% All values stay the same except "morphParamName" and "morphParamAxis",
% which increase from having 32 to 40 names
%groupName = info.groupName;
tName = info.tName;
%morphParamName = info.morphParamName;
%morphParamAxis = info.morphParamAxis;
nImagesGroup = info.nImagesGroup;
nCellsGroup = info.nCellsGroup;
nCells = info.nCells; 
nImages = info.nImages;

morphParamName = ["Area","Perimeter","Angle","Major","Minor","Height","Width","Mean","StdDev","Median","Skewness","Kurtosis","Feret","FeretAngle","MinFeret","FeretX","FeretY","Circularity","AR","Roundness","ConvexHullArea","Solidity","Eccentricity","ESF","ElongationFW","ElongationFmF","ConvexPerimeter","Convexity","FibreLength","FibreWidth","Curl","CurlW", "AreaNorm","PerimeterNorm","MajorNorm","MinorNorm","HeightNorm","WidthNorm", "FeretNorm","MinFeretNorm"];
morphParamAxis = ["Area (pixel^2)","Perimeter (pixel)","Angle (degree)","Major axis (pixel)","Minor axis (pixel)","Height (pixel)","Width (pixel)","Mean intensity","Standard deviation intensity","Median intensity","Skewness","Kurtosis","Feret (pixel)","Feret Angle (degree)","Minimum Feret (pixel)","FeretX (pixel)","FeretY (pixel)","Circularity","Aspect Ratio","Roundness","Convex hull area (pixel^2)","Solidity","Eccentricity","Elliptical shape factor","Elongation [Width/Feret]","Elongation [MinFeret/Feret]","Convex perimeter (pixel)","Convexity","Fibre length (pixel)","Fibre width (pixel)","Curl","CurlW", "Normalized Area","Normalized Perimeter","Normalized Major Axis","Normalized Minor Axis","Normalized Height","Normalized Width", "Normalized Feret","Normalized Minimum Feret"];

morphParam = mP.morphParam;

toc

%% 2) Normalize selected basic parameters
% 1. Area, 2. Perimeter, 4. Major, 5. Minor, 6. Height, 7. Width, 13. Feret, 15. MinFeret

tic
fprintf('\n Normalizing basic parameters for groupname %s \n', groupName);
fprintf('\t Number of cells in morphological parameter %s for UID 2 (arbitrary) = %d \n', morphParamName{23}, length(morphParam{1,2}{1,23}));
kList = [1,2,4,5,6,7,13,15];
kU = length(info.morphParamName); % Original number of morphological parameters

for kInd = 1:1:length(kList) % Iterate through all the parameters to normalize
    for uidI = 1:1:length(morphParam) % Iterate through all unique IDs or donors
        %imInd = 0; % Index of image (considering all images in a donor)
        % Iterate through all coverslips in the UID or donor
        for csI = 1:1:length(nCells{1,uidI})

            % Iterate through all images in the coverslip csI of uid
            for iI = 1:1:length(nCells{1,uidI}{1,csI})
                %imInd = imInd + 1;
                nCellsInImage = nCells{1,uidI}{1,csI}(iI);

                % cellIndI and cellIndF are the initial and final indices for
                % the cells in morphParam for the entire image
                if csI ==1 && iI == 1
                    cellIndI = 1;
                    cellIndF = nCellsInImage;
                    % Add new morphological parameter in kU + kInd place
                    morphParam{1,uidI}{1,(kU+kInd)} = morphParam{1,uidI}{1,kList(kInd)}(cellIndI:cellIndF)./mean(morphParam{1,uidI}{1,kList(kInd)}(cellIndI:cellIndF));
                else
                    cellIndI = cellIndI + nCellsInPrevImage;
                    cellIndF = cellIndI + nCellsInImage - 1;
                    % Concatenate to the newly created morphological
                    % parameter
                    morphParam{1,uidI}{1,(kU+kInd)} = cat(1,morphParam{1,uidI}{1,(kU+kInd)}, morphParam{1,uidI}{1,kList(kInd)}(cellIndI:cellIndF)./mean(morphParam{1,uidI}{1,kList(kInd)}(cellIndI:cellIndF)));
                end

                nCellsInPrevImage = nCellsInImage;
                
            end
        end
    end
    fprintf('\t Number of cells in morphological parameter %s for UID 2 (arbitrary) = %d \n', morphParamName{(kU+kInd)}, length(morphParam{1,2}{1,(kU+kInd)}));
end
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

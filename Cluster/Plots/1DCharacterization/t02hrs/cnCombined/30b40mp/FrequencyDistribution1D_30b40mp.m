% FrequencyDistribution1D_30b40mp.m
% Program to calculate frequency distribution of different morphological
% identifiers - 30 bins and 40 morphological parameters
% 2023 ps UBC
 
% Full list of morphological parameters: 
% Area,Perimeter,Angle,Major,Minor,Height,Width,Mean,StdDev,Median,Skewness,Kurtosis,Feret,FeretAngle,MinFeret,FeretX,FeretY,Circularity,AR,Roundness,ConvexHullArea,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,ConvexPerimeter,Convexity,FibreLength,FibreWidth,Curl,CurlW,AreaNorm,PerimeterNorm,MajorNorm,MinorNorm,HeightNorm,WidthNorm,FeretNorm,MinFeretNorm

% Can be run in linux or the cluster with the following command (example): 
% matlab -nodisplay -r "groupName = \"AA\";FrequencyDistribution1D_30b40mp"


% NOTE: The first 3 sections runs normalizeBasicParam in this program

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
 
tic 
% This implementation replaces the old variables with updated ones
fprintf("Saving general information variables.... \n");
save(strcat("../",groupName,"/generalInfo.mat"), "groupName", "tName", "morphParamName", "morphParamAxis", "nImages","nCells", "nImagesGroup", "nCellsGroup", '-v7.3');
 
fprintf("Saving variables for morphParam.... \n");
save(strcat("../",groupName, "/morphParamData.mat"),"morphParam", '-v7.3');

toc


% ------------------------------------------------------------------------
% HISTDATA CALCULATION
% ------------------------------------------------------------------------

%% 1) Load data - morphParamData, generalIno, histData
% groupName = 'AA'; % Example of groupName, which is inputted with the
% program in the shell script (in the cluster)

% NOTE: Uncomment if normalizeBasicParam is not run in this program
%{
tic
load(strcat("../",groupName,"/generalInfo.mat"));
load(strcat("../",groupName,"/morphParamData.mat"));
%load(strcat("../",groupName,"/histData.mat"));
toc
%}


%% 2) Find minimum and maximum values of all morphological parameters
tic
% Min and Max store minimum and maximum values for all morphological
% parameters
%{
Min = zeros(1,length(morphParamName));
Max = zeros(1,length(morphParamName));

for k = 1:1:length(morphParamName)
   
    clear("MinT"); clear("MaxT");
    %MinT and MaxT store minimum and maximum values of each file for
    %the specific morphological parameter
    MinT = zeros(1,length(morphParam));
    MaxT = zeros(1,length(morphParam));
    
    for j = 1:1:length(morphParam)
        MinT(j) = min(morphParam{1,j}{1,k});
        MaxT(j) = max(morphParam{1,j}{1,k});
    end

    Min(k) = min(MinT);
    Max(k) = max(MaxT);

end
%}

% Note: Min and Max are chosen based on values of all morphological
% parameters from Canada + Nepal data
% Area,Perimeter,Angle,Major,Minor,Height,Width,Mean,StdDev,Median,Skewness,Kurtosis,Feret,FeretAngle,MinFeret,FeretX,FeretY,Circularity,AR,Roundness,ConvexHullArea,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,ConvexPerimeter,Convexity,FibreLength,FibreWidth,Curl,CurlW,AreaNorm,PerimeterNorm,MajorNorm,MinorNorm,HeightNorm,WidthNorm,FeretNorm,MinFeretNorm
%Min = [0 0 0 0 0 0 0 0 0 0 -6 -5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%Max = [2000 250 180 120 45 115 115 250 105 255 6 30 115 180 55 3000 3000 1 30 1 2500 1 1 1 1 1 250 1 40 40 1.6 1.6];

Min = [0 40 0 10 5 5 5 100 5 90 -2.5 -2.5 10 0 5 0 0 0.2 0 0.2 0 0.5 0 0.2 0.2 0.2 40 0.7 5 10 1 0.5 0 0 0 0 0 0 0 0];
Max = [1200 160 180 60 35 65 55 170 65 170 2.5 5 60 180 45 3000 3000 1 5 1 1200 1.5 1 1 1.5 1 160 1.1 30 60 3 1.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5];

% Make a variable edges with equally spaced edges for histocount for all
% morphological parameters using max and min values above
nbins = 30;

edges = cell(1, length(morphParamName));
midBin = cell(1, length(morphParamName));
for k = 1:1:length(morphParamName)
    % edges includes edges for bins
    edges{k} = linspace(Min(k),Max(k),(nbins+1)); 
    
    % midBin is the middle value of each bin (used for plotting)
    midBin{k} = edges{k}(2:end) - (edges{k}(2)-edges{k}(1))/2;
    
end
toc

%% 3a) Calculate frequency distribution
% Save histocount data for each donor and each morphological identifier
% Donor-wise: histcount is calculated for ALL the data from the donor
% (including all coverslips)
% histDataD is a 2D cell, 1) donor, 2) morphological parameter
tic

histDataD30b40mp = cell(1, length(morphParam));
% Iterate through all files within groupName
fprintf('Calculating histcountsD for all 40 morophological parameters for group %s.\n',groupName);
for uidI = 1:1:length(morphParam)
    
    %Iterate through all morphological parameters (columns) within each file
    for k = 1:1:length(morphParamName)
        % histDataD is a 2D cell with histogram counts for all
        % the data. The two levels are 1) donor, and 2)
        % morphological parameter
        [histDataD30b40mp{1,uidI}{1,k}, edgesTemp] = histcounts(morphParam{1,uidI}{1,k}, edges{k}, 'Normalization','probability');
        
    end
    
end

toc


%% 3b) Calculate frequency distribution
% Save histocount data for each coverslip and each morphological identifier
% Coverslip-wise: histcount is calculated for each coverslip from the donor
% histDataCS is a 3D cell, 1) donor, 3) coverslip, 2) morphological 
% parameter
tic

histDataCS30b40mp = cell(1, length(morphParam));  % Length equal to number of donors or UIDs in this group
% Iterate through all files within groupName
fprintf('Calculating histcountsCS for all 40 morophological parameters for group %s.\n',groupName);
for uidI = 1:1:length(morphParam)

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
            % histDataD is a 2D cell with histogram counts for all
            % the data. The two levels are 1) donor, and 2)
            % morphological parameter
            [histDataCS30b40mp{1,uidI}{1,csI}{1,k}, edgesTemp] = histcounts(morphParam{1,uidI}{1,k}(cellIndexI:cellIndexF), edges{k}, 'Normalization','probability');

        end

    end
    
end


%% 3c) Calculate frequency distribution
% Save histocount data for each morphological identifier (for ALL data from
% all donors of a group)
% Group-wise: histcount is calculated for 
% histDataG is an array with histcount data for the group with groupName
tic

% Iterate through all files within groupName
fprintf('Calculating histcountsG for all 40 morophological parameters for group %s.\n',groupName);
histDataG30b40mp = cell(1, length(morphParamName));
morphParamG30b40mp = cell(1, length(morphParamName));
% Concatenate all the data into one cell 
for k = 1:1:length(morphParamName)
    morphParamG30b40mp{k} = morphParam{1,1}{1,k};
    for uidI = 2:1:length(morphParam)
        morphParamG30b40mp{k} = cat(1, morphParamG30b40mp{k}, morphParam{1,uidI}{1,k});
    end
    histDataG30b40mp{k} = histcounts(morphParamG30b40mp{1,k}, edges{k}, 'Normalization','probability');
end

toc


%% 4) Calculate mean and percentiles or std. dev. for all frequency
% distributions. Note, that they are averaged over different donors, and
% percentiles are for different donors
tic
meanHistD30b40mp = cell(1, length(morphParamName));
medianHistD30b40mp = cell(1, length(morphParamName));
perc25D30b40mp = cell(1, length(morphParamName));
perc75D30b40mp = cell(1, length(morphParamName));
perc5D30b40mp = cell(1, length(morphParamName));
perc95D30b40mp = cell(1, length(morphParamName));
stdHistD30b40mp = cell(1, length(morphParamName));
for k = 1:1:length(morphParamName)
    fprintf('%d out of %d) Calculating statistics for morphological parameter %s.\n',k,length(morphParamName),morphParamName(k));


    for p = 1:1:length(histDataD30b40mp{1,1}{1,k})
        clear('jTemp');
        jTemp = zeros(1,length(morphParam));
        
        for j = 1:1:length(morphParam)
            jTemp(j) = histDataD30b40mp{1,j}{1,k}(p);
            
        end
        meanHistD30b40mp{k}(p) = mean(jTemp);
        medianHistD30b40mp{k}(p) = median(jTemp);
        perc25D30b40mp{k}(p) = prctile(jTemp,25);
        perc75D30b40mp{k}(p) = prctile(jTemp,75);
        perc5D30b40mp{k}(p) = prctile(jTemp,5);
        perc95D30b40mp{k}(p) = prctile(jTemp,95);
        stdHistD30b40mp{k}(p) = std(jTemp);
        
    end

end
toc


%% 5) Save variables for histData
tic

% Create a folder with the groupName if one does not exist
if not(isfolder(groupName))
    mkdir(groupName)
end

fprintf("Saving histData.... \n");
save(strcat(groupName,"/histData.mat"), "edges", "midBin", "histDataD30b40mp", "histDataCS30b40mp", "histDataG30b40mp", "meanHistD30b40mp", "medianHistD30b40mp", "perc25D30b40mp", "perc75D30b40mp", "perc5D30b40mp", "perc95D30b40mp", "stdHistD30b40mp", '-v7.3');

fprintf("Saving general information variables.... \n");
save(strcat(groupName,"/generalInfo.mat"), "groupName", "tName", "morphParamName", "morphParamAxis", "nImages","nCells", "nImagesGroup", "nCellsGroup", '-v7.3');

toc

%% 6a) Calculate frequency distribution
% Save histocount data for each image and each morphological identifier
% Image-wise: histcount is calculated for each image from the donor
% histDataI is a 3D cell, 1) donor, 3) image, 2) morphological 
% parameter
tic

histDataI30b40mp = cell(1, length(morphParam));  % Length equal to number of donors or UIDs in this group
% Iterate through all files within groupName
fprintf('Calculating histcountsI for all 40 morophological parameters for group %s.\n',groupName);
for uidI = 1:1:length(morphParam)
    imInd = 0; % Index of image (considering all images in a donor)
    % Iterate through all coverslips in the UID or donor
    for csI = 1:1:length(nCells{1,uidI})

        % Iterate through all images in the coverslip csI of uid 
        for iI = 1:1:length(nCells{1,uidI}{1,csI})
            imInd = imInd + 1;
            nCellsInImage = nCells{1,uidI}{1,csI}(iI);

            % cellIndI and cellIndF are the initial and final indices for
            % the cells in morphParam for the entire image
            if csI ==1 && iI == 1
                cellIndI = 1;
                cellIndF = nCellsInImage;
            else
                cellIndI = cellIndI + nCellsInPrevImage;
                cellIndF = cellIndI + nCellsInImage - 1;
            end

            nCellsInPrevImage = nCellsInImage;

            %Iterate through all morphological parameters (columns) within each file
            for k = 1:1:length(morphParamName)
                % histDataD is a 2D cell with histogram counts for all
                % the data. The three levels are 1) donor, 3) image and 2)
                % morphological parameter
                [histDataI30b40mp{1,uidI}{1,imInd}{1,k}, edgesTemp] = histcounts(morphParam{1,uidI}{1,k}(cellIndI:cellIndF), edges{k}, 'Normalization','probability');
            end
        end

    end

end

%% 6b) Save variables for histDataI
tic

fprintf("Saving histDataI.... \n");
save(strcat(groupName,"/histDataI.mat"), "histDataI30b40mp", '-v7.3');

toc
 
%% 7_0) Figures
 
if not(isfolder("/Figures"))
    mkdir("Figures");
end

if not(isfolder(strcat("Figures/",groupName)))
    mkdir(strcat("Figures/",groupName))
end
 
%% 7a) Plot frequency distribution for all donors and group data 
tic 
fprintf('Plotting frequency distribution donor-wise %s.\n',groupName);
%kValue = [18, 23]; %18: Circularity; 23: Eccentricity
% color scheme: 
% [27/255 158/255 119/255]: AA 
% [217/255 95/255 2/255]: ABeta
% [117/255 112/255 179/255]: AS
% [231/255 41/255 138/255]: SCD
pGroups = ["AA", "ABeta", "AS", "SCD"];
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [231/255 41/255 138/255], [166/255 118/255 29/255],[35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2

% Find the index for the pGroups, which defines the index for the color
colorI = find(contains(pGroups, groupName));
 
%for kk = 1:1:length(kValue)
for k = 1:1:length(morphParamName)
    %k = kValue(kk);
    figure; hold on
    for j = 1:1:length(histDataD30b40mp)
        plot(midBin{k}, histDataD30b40mp{1,j}{1,k}, 'color', [.5 .5 .5 0.3]);
        xlabel(morphParamAxis{k}); ylabel("Normalized frequency");
        
        title(groupName);
        
    end
    plot(midBin{k}, histDataG30b40mp{k}, 'color', colors{colorI}, 'LineWidth', 1.5);
    
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat('Figures/',groupName,'/2_Donorwise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat('Figures/',groupName,'/2_Donorwise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc
%% 7b) Plot frequency distribution for all coverslip and group data
tic 
fprintf('Plotting frequency distribution coverslip-wise %s.\n',groupName);
%kValue = [18, 23]; %18: Circularity; 23: Eccentricity
 
% color scheme: 
% [27/255 158/255 119/255]: AA 
% [217/255 95/255 2/255]: ABeta
% [117/255 112/255 179/255]: AS
% [231/255 41/255 138/255]: SCD
pGroups = ["AA", "ABeta", "AS", "SCD"];
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [231/255 41/255 138/255], [166/255 118/255 29/255],[35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2

% Find the index for the pGroups, which defines the index for the color
colorI = find(contains(pGroups, groupName));
 
%for kk = 1:1:length(kValue)
for k = 1:1:length(morphParamName)
    %k = kValue(kk);
    figure; hold on
    for uidI = 1:1:length(histDataCS30b40mp)
        for csI = 1:1:length(histDataCS30b40mp{1,uidI})
            plot(midBin{k}, histDataCS30b40mp{1,uidI}{1,csI}{1,k}, 'color', [.5 .5 .5 0.3]);
            xlabel(morphParamAxis{k}); ylabel("Normalized frequency");
 
            title(groupName);
        end
        
    end
    plot(midBin{k}, histDataG30b40mp{k}, 'color', colors{colorI}, 'LineWidth', 1.5);
    
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat('Figures/',groupName,'/3_Coverslip-wise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat('Figures/',groupName,'/3_Coverslip-wise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc

 
%% 7c) Compare frequency distribution for mean and overall histcount for a few morphological parameters
tic 
fprintf('Plotting frequency distribution (only mean and overall or group-wise) for %s.\n',groupName);
kValue = [1,18,23]; %18: Circularity; 23: Eccentricity
 
pGroups = ["AA", "ABeta", "AS", "SBeta", "SS", "N", "SCT", "SCD", "Test2"];
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [166/255 118/255 29/255], [231/255 41/255 138/255], [35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
 
% Find the index for the pGroups, which defines the index for the color
colorI = find(contains(pGroups, groupName));
 
for kk = 1:1:length(kValue)
    k = kValue(kk);
    figure; hold on
    plot(midBin{k}, histDataG30b40mp{k}, 'color', colors{colorI});
    plot(midBin{k}, meanHistD30b40mp{k}, '--','color', colors{colorI});
    plot(midBin{k}, medianHistD30b40mp{k}, '-.','color', colors{colorI});
    title(groupName);
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat('Figures/',groupName,'/0_Mean-vs-overallHist-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat('Figures/',groupName,'/0_Mean-vs-overallHist-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc

%% 7d) Plot frequency distribution for all image and group data
tic 
fprintf('Plotting frequency distribution image-wise %s.\n',groupName);
%kValue = [18, 23]; %18: Circularity; 23: Eccentricity
 
% color scheme: 
% [27/255 158/255 119/255]: AA 
% [217/255 95/255 2/255]: ABeta
% [117/255 112/255 179/255]: AS
% [231/255 41/255 138/255]: SCD
pGroups = ["AA", "ABeta", "AS", "SCD"];
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [231/255 41/255 138/255], [166/255 118/255 29/255],[35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2

% Find the index for the pGroups, which defines the index for the color
colorI = find(contains(pGroups, groupName));
 
%for kk = 1:1:length(kValue)
for k = 1:1:length(morphParamName)
    %k = kValue(kk);
    figure; hold on
    for uidI = 1:1:length(histDataI30b40mp)
        for iI = 1:1:length(histDataI30b40mp{1,uidI})
            plot(midBin{k}, histDataI30b40mp{1,uidI}{1,iI}{1,k}, 'color', [.5 .5 .5 0.3]);
            xlabel(morphParamAxis{k}); ylabel("Normalized frequency");
 
            title(groupName);
        end
        
    end
    plot(midBin{k}, histDataG30b40mp{k}, 'color', colors{colorI}, 'LineWidth', 1.5);
    
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat('Figures/',groupName,'/4_Image-wise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat('Figures/',groupName,'/4_Image-wise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc

% FrequencyDistribution1Dp.m
% Program to calculate frequency distribution of different morphological
% identifiers
% 2023 ps UBC
 
% Full list of morphological parameters: 
% Area,Perimeter,Angle,Major,Minor,Height,Width,Mean,StdDev,Median,Skewness,Kurtosis,Feret,FeretAngle,MinFeret,FeretX,FeretY,Circularity,AR,Roundness,ConvexHullArea,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,ConvexPerimeter,Convexity,FibreLength,FibreWidth,Curl,CurlW
 
 
%% 1) Read .csv files with measured parameters from ImageJ or FIJI
tic
 
%dataPath = "../../../MorphologicalParameters/Nepal/";
 
% Different classes/categories/groups
%groupName = "Test2";
%tName = "t02hrs";
 
% List of relevant morphological parameters or identifiers from ImageJ or
% FIJI 
morphParamName = ["Area","Perimeter","Angle","Major","Minor","Height","Width","Mean","StdDev","Median","Skewness","Kurtosis","Feret","FeretAngle","MinFeret","FeretX","FeretY","Circularity","AR","Roundness","ConvexHullArea","Solidity","Eccentricity","ESF","ElongationFW","ElongationFmF","ConvexPerimeter","Convexity","FibreLength","FibreWidth","Curl","CurlW"];
morphParamAxis = ["Area (pixel^2)","Perimeter (pixel)","Angle (degree)","Major (pixel)","Minor (pixel)","Height (pixel)","Width (pixel)","Mean intensity","Standard deviation intensity","Median intensity","Skewness","Kurtosis","Feret (pixel)","FeretAngle (degree)","MinFeret (pixel)","FeretX (pixel)","FeretY (pixel)","Circularity","Aspect Ratio","Roundness","Convex hull area (pixel^2)","Solidity","Eccentricity","Elliptical shape factor","Elongation [Width/Feret]","Elongation [MinFeret/Feret]","Convex perimeter (pixel)","Convexity","Fibre length (pixel)","Fibre width (pixel)","Curl","CurlW"];
 
 
fprintf('Reading data from group %s and tname %s \n',groupName, tName);
 
%clear('uidName');
% uidName is a list of all UIDs or unique IDs of group(i)
% e.g. uidName could include AA-001, AA-002, ..., if group(i) = "AA"
uidName = dir(strcat(dataPath,groupName,'*')); 
 
% Loop through all the UIDs or unique IDs for the group groupName
parfor uidI=1:1:length(uidName)
    
    [morphParamT, nImagesT, nCellsT, strM] = readCSV(tName, uidName, uidI, dataPath, morphParamName);
    
    % morphParam is a 2D cell with morphological parameters for all
    % the data. The two levels are 1) donor or UID, and 2)
    % morphological parameter or k
    
    morphParam{1,uidI} = morphParamT;
    nImages{1,uidI} = nImagesT;
    nCells{1,uidI} = nCellsT;
    
    fprintf(strM);
    
end
 
toc
 
%% 2) Calculate and print total number of cells and images for the entire group
nImagesGroup = 0;
nCellsGroup = 0;
 
for jj = 1:1:length(nImages)
    nImagesGroup = nImagesGroup + sum(nImages{1,jj});
    
    for kk = 1:1:length(nCells{1,jj})
        nCellsGroup =  nCellsGroup + sum(nCells{1,jj}{1,kk});
        
    end
end
fprintf('\nTotal number of images for group %s = %d \n', groupName, nImagesGroup);
fprintf('Total number of cells for group %s = %d \n', groupName, nCellsGroup);
 
 
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
 
%% 4) Find minimum and maximum values of all morphological parameters
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
% parameters from Canada data
% Area,Perimeter,Angle,Major,Minor,Height,Width,Mean,StdDev,Median,Skewness,Kurtosis,Feret,FeretAngle,MinFeret,FeretX,FeretY,Circularity,AR,Roundness,ConvexHullArea,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,ConvexPerimeter,Convexity,FibreLength,FibreWidth,Curl,CurlW
Min = [0 0 0 0 0 0 0 0 0 0 -6 -5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
Max = [2000 250 180 120 45 115 115 250 105 255 6 30 115 180 55 3000 3000 1 30 1 2500 1 1 1 1 1 250 1 40 40 1.6 1.6];
 
% Make a variable edges with equally spaced edges for histocount for all
% morphological parameters using max and min values above
nbins = 100;
 
edges = cell(1, length(morphParamName));
midBin = cell(1, length(morphParamName));
for k = 1:1:length(morphParamName)
    % edges includes edges for bins
    edges{k} = linspace(Min(k),Max(k),(nbins+1)); 
    
    % midBin is the middle value of each bin (used for plotting)
    midBin{k} = edges{k}(2:end) - (edges{k}(2)-edges{k}(1))/2;
    
end
toc
 
%% 5) Calculate frequency distribution
% Save histocount data for each donor and each morphological identifier
% Donor-wise: histcount is calculated for ALL the data from the donor
% (including all coverslips)
% histDataD is a 2D cell, 1) donor, 2) morphological parameter
tic
 
histDataD = cell(1, length(uidName));
% Iterate through all files within groupName
fprintf('Calculating histcountsD for group %s.\n',groupName);
for uidI = 1:1:length(morphParam)
    
    %Iterate through all morphological parameters (columns) within each file
    for k = 1:1:length(morphParamName)
        % histDataD is a 2D cell with histogram counts for all
        % the data. The two levels are 1) donor, and 2)
        % morphological parameter
        [histDataD{1,uidI}{1,k}, edgesTemp] = histcounts(morphParam{1,uidI}{1,k}, edges{k}, 'Normalization','probability');
        
    end
    
end
 
toc
 
%% 5b) Calculate frequency distribution
% Save histocount data for each coverslip and each morphological identifier
% Coverslip-wise: histcount is calculated for each coverslip from the donor
% histDataCS is a 3D cell, 1) donor, 3) coverslip, 2) morphological 
% parameter
tic
 
histDataCS = cell(1, length(uidName));  % Length equal to number of donors or UIDs in this group
% Iterate through all files within groupName
fprintf('Calculating histcountsCS for group %s.\n',groupName);
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
            [histDataCS{1,uidI}{1,csI}{1,k}, edgesTemp] = histcounts(morphParam{1,uidI}{1,k}(cellIndexI:cellIndexF), edges{k}, 'Normalization','probability');
 
        end
 
    end
    
end
 
toc
 
%% 5c) Calculate frequency distribution
% Save histocount data for each morphological identifier (for ALL data from
% all donors of a group)
% Group-wise: histcount is calculated for 
% histDataG is an array with histcount data for the group with groupName
tic
 
% Iterate through all files within groupName
fprintf('Calculating histcountsG for group %s.\n',groupName);
histDataG = cell(1, length(morphParamName));
morphParamG = cell(1, length(morphParamName));
% Concatenate all the data into one cell 
for k = 1:1:length(morphParamName)
    morphParamG{k} = morphParam{1,1}{1,k};
    for uidI = 2:1:length(morphParam)
        morphParamG{k} = cat(1, morphParamG{k}, morphParam{1,uidI}{1,k});
    end
    histDataG{k} = histcounts(morphParamG{k}, edges{k}, 'Normalization','probability');
end
 
toc
 
 
%% 6) Calculate mean and percentiles or std. dev. for all frequency
% distributions. Note, that they are averaged over different donors, and
% percentiles are for different donors
tic
meanHistD = cell(1, length(morphParamName));
perc25D = cell(1, length(morphParamName));
perc75D = cell(1, length(morphParamName));
perc5D = cell(1, length(morphParamName));
perc95D = cell(1, length(morphParamName));
stdHistD = cell(1, length(morphParamName));
for k = 1:1:length(morphParamName)
    fprintf('%d out of %d) Calculating statistics for morphological parameter %s.\n',k,length(morphParamName),morphParamName(k));
 
 
    for p = 1:1:length(histDataD{1,1}{1,k})
        clear('jTemp');
        jTemp = zeros(1,length(morphParam));
        
        for j = 1:1:length(morphParam)
            jTemp(j) = histDataD{1,j}{1,k}(p);
            
        end
        meanHistD{k}(p) = mean(jTemp);
        perc25D{k}(p) = prctile(jTemp,25);
        perc75D{k}(p) = prctile(jTemp,75);
        perc5D{k}(p) = prctile(jTemp,5);
        perc95D{k}(p) = prctile(jTemp,95);
        stdHistD{k}(p) = std(jTemp);
        
    end
 
end
toc 
 
%% 7) Save variables for histData
tic
 
fprintf("Saving histData.... \n");
save(strcat(groupName,"/histData.mat"), "edges", "midBin", "histDataD", "histDataCS", "histDataG", "meanHistD", "perc25D", "perc75D", "perc5D", "perc95D", "stdHistD", '-v7.3');
toc
 
%% 8) Figures
 
if not(isfolder(strcat(groupName,"/Figures")))
    mkdir(strcat(groupName,"/Figures"))
end
 
%% 8a) Plot frequency distribution for all donors and group data 
 
%kValue = [18, 23]; %18: Circularity; 23: Eccentricity
% color scheme: 
% [27/255 158/255 119/255]: AA 
% [217/255 95/255 2/255]: ABeta
% [117/255 112/255 179/255]: AS
% [166/255 118/255 29/255]: SBeta
% [231/255 41/255 138/255]: SS
% [35/255 110/255 90/255]: N
% [115/255 85/255 140/255]: SCT
% [180/255 40/255 90/255]: SCD
pGroups = ["AA", "ABeta", "AS", "SBeta", "SS", "N", "SCT", "SCD", "Test2"];
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [166/255 118/255 29/255], [231/255 41/255 138/255], [35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
 
% Find the index for the pGroups, which defines the index for the color
colorI = find(contains(pGroups, groupName));
 
%for kk = 1:1:length(kValue)
for k = 1:1:length(morphParamName)
    %k = kValue(kk);
    figure; hold on
    for j = 1:1:length(histDataD)
        plot(midBin{k}, histDataD{1,j}{1,k}, 'color', [.5 .5 .5]);
        xlabel(morphParamAxis{k}); ylabel("Normalized frequency");
        
        title(groupName);
        
    end
    plot(midBin{k}, histDataG{k}, 'color', colors{colorI}, 'LineWidth', 1.5);
    
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat(groupName,'/Figures/2_Donorwise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat(groupName,'/Figures/2_Donorwise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
 
%% 8b) Plot frequency distribution for all coverslip and group data
 
%kValue = [18, 23]; %18: Circularity; 23: Eccentricity
 
pGroups = ["AA", "ABeta", "AS", "SBeta", "SS", "N", "SCT", "SCD", "Test2"];
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [166/255 118/255 29/255], [231/255 41/255 138/255], [35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
 
% Find the index for the pGroups, which defines the index for the color
colorI = find(contains(pGroups, groupName));
 
%for kk = 1:1:length(kValue)
for k = 1:1:length(morphParamName)
    %k = kValue(kk);
    figure; hold on
    for uidI = 1:1:length(histDataCS)
        for csI = 1:1:length(histDataCS{1,uidI})
            plot(midBin{k}, histDataCS{1,uidI}{1,csI}{1,k}, 'color', [.5 .5 .5]);
            xlabel(morphParamAxis{k}); ylabel("Normalized frequency");
 
            title(groupName);
        end
        
    end
    plot(midBin{k}, histDataG{k}, 'color', colors{colorI}, 'LineWidth', 1.5);
    
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat(groupName,'/Figures/1_Coverslip-wise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat(groupName,'/Figures/1_Coverslip-wise-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
 
%% 8c) Compare frequency distribution for mean and overall histcount for a few morphological parameters
 
kValue = [8, 11, 12, 13, 18, 19, 20, 22, 23, 24, 28]; %18: Circularity; 23: Eccentricity
 
pGroups = ["AA", "ABeta", "AS", "SBeta", "SS", "N", "SCT", "SCD", "Test2"];
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [166/255 118/255 29/255], [231/255 41/255 138/255], [35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
 
% Find the index for the pGroups, which defines the index for the color
colorI = find(contains(pGroups, groupName));
 
for kk = 1:1:length(kValue)
    k = kValue(kk);
    figure; hold on
    plot(midBin{k}, histDataG{k}, 'color', colors{colorI});
    plot(midBin{k}, meanHistD{k}, '--','color', colors{colorI});
    title(groupName);
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat(groupName,'/Figures/0_Mean-vs-overallHist-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat(groupName,'/Figures/0_Mean-vs-overallHist-',morphParamName{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
 
% FUNCTIONS
 
 
function [morphParamT, nImagesT, nCellsT, strM] = readCSV(tName, uidName, uidI, dataPath, morphParamName)
 
    timeStart = tic;
    %fprintf('\t %d / %d) Reading data from UID %s\n',uidI,length(uidName),uidName(uidI).name);
    strM = strcat('\t ', num2str(uidI), '/', num2str(length(uidName)), ') Completed reading data from UID ', uidName(uidI).name, '\n');
 
    %clear('csidFullL');
    % csidFullL is a list of all the coverslip IDs or CSIDs for the
    % particular UID
    csidFullL = dir(strcat(dataPath, uidName(uidI).name, '/',tName,'/*_*'));
 
    % Remove CSID starting from TS - for temperature study
    TSi = 0; % Number of files with TS name
    TSList = zeros(1,length(csidFullL));
    for ii = 1:1:length(csidFullL)
        if(contains(csidFullL(ii).name,"TS"))
            fprintf(strcat("\t \t *** Not including CSID: ", csidFullL(ii).name, '\n'));
            TSi = TSi + 1;
            TSList(TSi) = ii;
 
        end
    end
    csidFullL(nonzeros(TSList)) = []; %Remove element with TS
 
    % nImages stores the number of images in each coverslips
    morphParamT = cell(1, length(csidFullL));
    nImagesT = zeros(1,length(csidFullL));
    %nImagesT = cell(1, length(uidName));
    nCellsT = cell(1, length(csidFullL));
 
    % Read files (morphological parameters per image) for each CSID
    % Loop through all CSIDs or coverslip IDs in the UID
    for csidI = 1:1:length(csidFullL)
        %clear('mpList');
        % mpList is a list of all the files containing morphological
        % parameter data (.csv files) for this particular CSID
        mpList =  dir(strcat(dataPath, uidName(uidI).name, '/',tName,'/', csidFullL(csidI).name,'/*.csv'));
        strM = strcat(strM, '\t \t Number of files in ', csidFullL(csidI).name, ' = ', num2str(length(mpList)), '\n');
 
        nImagesT(csidI) = length(mpList);
 
        nCellsT{1,csidI} = zeros(1,length(mpList));
 
        % Loop through all the files in the CSID
        for mpI = 1:1:length(mpList)
            %clear('data');
            warning('off','MATLAB:table:ModifiedAndSavedVarnames');
            fileName = strcat(dataPath, uidName(uidI).name, '/t02hrs/', csidFullL(csidI).name,'/', mpList(mpI).name);
            data = readtable(fileName);
 
            nCellsT{1,csidI}(mpI) = length(data.Var1);
 
            %Iterate through all morphological parameters (columns) within each file
            for k = 1:1:length(morphParamName)
                % morphParam is a 2D cell with morphological parameters for all
                % the data. The two levels are 1) donor or UID, and 2)
                % morphological parameter or k
                
                % morphParamT is a temporary cell used only in this
                % function to save the values for all the morphological
                % parameters in this specific UID or donor
 
                % For the first file of the first coverslip of each
                % donor, add the morphological parameter data, while
                % for all the rest, concatenate the data
 
                if csidI ==1 && mpI == 1
                    %morphParam{1,i}{1,uidI}{1,k}= eval(strcat('data.',morphParamName(k))); %Save morphological parameter value (array) into a 3D cell
                    morphParamT{1,k}= eval(strcat('data.',morphParamName(k))); %Save morphological parameter value (array) into a 3D cell
                else
                    %morphParam{1,i}{1,uidI}{1,k} = cat(1,morphParam{1,i}{1,uidI}{1,k}, eval(strcat('data.',morphParamName(k))));  % Concatenate the arrays along direction 1 (keep adding rows - in column 1)
                    morphParamT{1,k} = cat(1,morphParamT{1,k}, eval(strcat('data.',morphParamName(k))));  % Concatenate the arrays along direction 1 (keep adding rows - in column 1)
                end
 
            end
 
        end
 
    end
    timeUID = toc(timeStart);
    %fprintf('\t Time for reading all data in UID %s = %d seconds \n', uidName(uidI).name, timeUID);
    strM = strcat(strM,  '\t Time for reading all data in UID', uidName(uidI).name, ' = ', num2str(timeUID), '\n');
 
end


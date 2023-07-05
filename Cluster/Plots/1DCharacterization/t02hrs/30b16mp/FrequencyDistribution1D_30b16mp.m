% FrequencyDistribution1D_30b16mp.m
% Program to calculate frequency distribution of different morphological
% identifiers - 30 bins and 16 morphological parameters
% 2023 ps UBC
 
% Full list of morphological parameters: 
% Area,Perimeter,Angle,Major,Minor,Height,Width,Mean,StdDev,Median,Skewness,Kurtosis,Feret,FeretAngle,MinFeret,FeretX,FeretY,Circularity,AR,Roundness,ConvexHullArea,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,ConvexPerimeter,Convexity,FibreLength,FibreWidth,Curl,CurlW
 
% Selected list of 16 morphological parameters: 
% Area,Perimeter,Angle,Major,Minor,Height,Width,Feret,MinFeret,Circularity,Roundness,Solidity,Eccentricity,ESF,ElongationFmF,Convexity

% Can be run in linux or the cluster with the following command (example): 
% matlab -nodisplay -r "groupName = \"AA\";tName = \"t02hrs\";FrequencyDistribution1D_30b16mp"


%% 0) List of 16 morphological parameters
% List of 16 morphological parameters or identifiers from ImageJ or FIJI
morphParamName16 = ["Area","Perimeter","Angle","Major","Minor","Height","Width","Feret","MinFeret","Circularity","Roundness","Solidity","Eccentricity","ESF","ElongationFmF","Convexity"];
morphParamAxis16 = ["Area (pixel^2)","Perimeter (pixel)","Angle (degree)","Major (pixel)","Minor (pixel)","Height (pixel)","Width (pixel)","Feret (pixel)","MinFeret (pixel)","Circularity","Roundness","Solidity","Eccentricity","Elliptical shape factor","Elongation [MinFeret/Feret]","Convexity"];

% Testing parameters
%groupName = 'Test2';
%tName = 't02hrs';

%% 1) Load data - morphParamData, generalIno, histData
% groupName = 'AA'; % Example of groupName, which is inputted with the
% program in the shell script (in the cluster)
tic
load(strcat("../",groupName,"/generalInfo.mat"));
load(strcat("../",groupName,"/morphParamData.mat"));
load(strcat("../",groupName,"/histData.mat"));
toc

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
% parameters from Canada data
% Area,Perimeter,Angle,Major,Minor,Height,Width,Mean,StdDev,Median,Skewness,Kurtosis,Feret,FeretAngle,MinFeret,FeretX,FeretY,Circularity,AR,Roundness,ConvexHullArea,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,ConvexPerimeter,Convexity,FibreLength,FibreWidth,Curl,CurlW
Min = [0 0 0 0 0 0 0 0 0 0 -6 -5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
Max = [2000 250 180 120 45 115 115 250 105 255 6 30 115 180 55 3000 3000 1 30 1 2500 1 1 1 1 1 250 1 40 40 1.6 1.6];

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

histDataD30b32mp = cell(1, length(morphParam));
% Iterate through all files within groupName
fprintf('Calculating histcountsD for all 32 morophological parameters for group %s.\n',groupName);
for uidI = 1:1:length(morphParam)
    
    %Iterate through all morphological parameters (columns) within each file
    for k = 1:1:length(morphParamName)
        % histDataD is a 2D cell with histogram counts for all
        % the data. The two levels are 1) donor, and 2)
        % morphological parameter
        [histDataD30b32mp{1,uidI}{1,k}, edgesTemp] = histcounts(morphParam{1,uidI}{1,k}, edges{k}, 'Normalization','probability');
        
    end
    
end

toc


%% 3b) Calculate frequency distribution
% Save histocount data for each coverslip and each morphological identifier
% Coverslip-wise: histcount is calculated for each coverslip from the donor
% histDataCS is a 3D cell, 1) donor, 3) coverslip, 2) morphological 
% parameter
tic

histDataCS30b32mp = cell(1, length(morphParam));  % Length equal to number of donors or UIDs in this group
% Iterate through all files within groupName
fprintf('Calculating histcountsCS for all 32 morophological parameters for group %s.\n',groupName);
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
            [histDataCS30b32mp{1,uidI}{1,csI}{1,k}, edgesTemp] = histcounts(morphParam{1,uidI}{1,k}(cellIndexI:cellIndexF), edges{k}, 'Normalization','probability');

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
fprintf('Calculating histcountsG for all 32 morophological parameters for group %s.\n',groupName);
histDataG30b32mp = cell(1, length(morphParamName));
morphParamG30b32mp = cell(1, length(morphParamName));
% Concatenate all the data into one cell 
for k = 1:1:length(morphParamName)
    morphParamG30b32mp{k} = morphParam{1,1}{1,k};
    for uidI = 2:1:length(morphParam)
        morphParamG30b32mp{k} = cat(1, morphParamG30b32mp{k}, morphParam{1,uidI}{1,k});
    end
    histDataG30b32mp{k} = histcounts(morphParamG30b32mp{1,k}, edges{k}, 'Normalization','probability');
end

toc

%% 3d) Retain frequency distribution for 16 morphological parameters
% Save histcount data for 16 morphological parameters at all levels (cover
% slip, donor and group)
tic
fprintf('Calculating histcounts for 16 morophological parameters for group %s.\n',groupName);
k16 = [1,2,3,4,5,6,7,13,15,18,20,22,23,24,26,28]; % Area,Perimeter,Angle,Major,Minor,Height,Width,Feret,MinFeret,Circularity,Roundness,Solidity,Eccentricity,ESF,ElongationFmF,Convexity

% Initialize cells for histData at different levels
histDataCS30b16mp = cell(1, length(morphParam));;
histDataD30b16mp = cell(1, length(morphParam));;
histDataG30b16mp = cell(1, length(morphParam));;
edges16 = cell(1, length(k16));
midBin16 = cell(1, length(k16));

for ki = 1:1:length(k16)
    for uidI = 1:1:length(morphParam)
        
        histDataD30b16mp{1,uidI}{1,ki} = histDataD30b32mp{1,uidI}{1, k16(ki)};
        
        for csI = 1:1:length(nCells{1,uidI})
            histDataCS30b16mp{1,uidI}{1,csI}{1,ki}= histDataCS30b32mp{1,uidI}{1,csI}{1, k16(ki)};
        end

    end
    histDataG30b16mp{1,ki} = histDataG30b32mp{1, k16(ki)};
    edges16{ki} = edges{k16(ki)};
    midBin16{ki} = midBin{k16(ki)};
end
toc
%% 4) Calculate mean and percentiles or std. dev. for all frequency
% distributions. Note, that they are averaged over different donors, and
% percentiles are for different donors
tic
meanHistD30b32mp = cell(1, length(morphParamName));
perc25D30b32mp = cell(1, length(morphParamName));
perc75D30b32mp = cell(1, length(morphParamName));
perc5D30b32mp = cell(1, length(morphParamName));
perc95D30b32mp = cell(1, length(morphParamName));
stdHistD30b32mp = cell(1, length(morphParamName));
for k = 1:1:length(morphParamName)
    fprintf('%d out of %d) Calculating statistics for morphological parameter %s.\n',k,length(morphParamName),morphParamName(k));


    for p = 1:1:length(histDataD30b32mp{1,1}{1,k})
        clear('jTemp');
        jTemp = zeros(1,length(morphParam));
        
        for j = 1:1:length(morphParam)
            jTemp(j) = histDataD30b32mp{1,j}{1,k}(p);
            
        end
        meanHistD30b32mp{k}(p) = mean(jTemp);
        perc25D30b32mp{k}(p) = prctile(jTemp,25);
        perc75D30b32mp{k}(p) = prctile(jTemp,75);
        perc5D30b32mp{k}(p) = prctile(jTemp,5);
        perc95D30b32mp{k}(p) = prctile(jTemp,95);
        stdHistD30b32mp{k}(p) = std(jTemp);
        
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
save(strcat(groupName,"/histData.mat"), "edges", "midBin", "histDataD30b32mp", "histDataCS30b32mp", "histDataG30b32mp", "meanHistD30b32mp", "perc25D30b32mp", "perc75D30b32mp", "perc5D30b32mp", "perc95D30b32mp", "stdHistD30b32mp", "edges16", "midBin16", "histDataD30b16mp", "histDataCS30b16mp", "histDataG30b16mp", '-v7.3');

fprintf("Saving general information variables.... \n");
save(strcat(groupName,"/generalInfo.mat"), "groupName", "tName", "morphParamName", "morphParamAxis", "nImages","nCells", "nImagesGroup", "nCellsGroup", "morphParamName16", "morphParamAxis16", '-v7.3');

toc
 
%% 6) Figures
 
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
for k = 1:1:length(morphParamName16)
    %k = kValue(kk);
    figure; hold on
    for j = 1:1:length(histDataD30b16mp)
        plot(midBin16{k}, histDataD30b16mp{1,j}{1,k}, 'color', [.5 .5 .5]);
        xlabel(morphParamAxis16{k}); ylabel("Normalized frequency");
        
        title(groupName);
        
    end
    plot(midBin16{k}, histDataG30b16mp{k}, 'color', colors{colorI}, 'LineWidth', 1.5);
    
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat('Figures/',groupName,'/2_Donorwise-',morphParamName16{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat('Figures/',groupName,'/2_Donorwise-',morphParamName16{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc
%% 7b) Plot frequency distribution for all coverslip and group data
tic 
fprintf('Plotting frequency distribution coverslip-wise %s.\n',groupName);
%kValue = [18, 23]; %18: Circularity; 23: Eccentricity
 
pGroups = ["AA", "ABeta", "AS", "SBeta", "SS", "N", "SCT", "SCD", "Test2"];
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [166/255 118/255 29/255], [231/255 41/255 138/255], [35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
 
% Find the index for the pGroups, which defines the index for the color
colorI = find(contains(pGroups, groupName));
 
%for kk = 1:1:length(kValue)
for k = 1:1:length(morphParamName16)
    %k = kValue(kk);
    figure; hold on
    for uidI = 1:1:length(histDataCS30b16mp)
        for csI = 1:1:length(histDataCS30b16mp{1,uidI})
            plot(midBin16{k}, histDataCS30b16mp{1,uidI}{1,csI}{1,k}, 'color', [.5 .5 .5]);
            xlabel(morphParamAxis16{k}); ylabel("Normalized frequency");
 
            title(groupName);
        end
        
    end
    plot(midBin16{k}, histDataG30b16mp{k}, 'color', colors{colorI}, 'LineWidth', 1.5);
    
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat('Figures/',groupName,'/1_Coverslip-wise-',morphParamName16{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat('Figures/',groupName,'/1_Coverslip-wise-',morphParamName16{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc
 
%% 7c) Compare frequency distribution for mean and overall histcount for a few morphological parameters
tic 
fprintf('Plotting frequency distribution (only mean and overall or group-wise) for %s.\n',groupName);
kValue = [10,11,12,13,14]; %18: Circularity; 23: Eccentricity
 
pGroups = ["AA", "ABeta", "AS", "SBeta", "SS", "N", "SCT", "SCD", "Test2"];
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [166/255 118/255 29/255], [231/255 41/255 138/255], [35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
 
% Find the index for the pGroups, which defines the index for the color
colorI = find(contains(pGroups, groupName));
 
for kk = 1:1:length(kValue)
    k = kValue(kk);
    figure; hold on
    plot(midBin16{k}, histDataG30b16mp{k}, 'color', colors{colorI});
    plot(midBin16{k}, meanHistD30b32mp{k16(k)}, '--','color', colors{colorI});
    title(groupName);
    hold off
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    fnamePNG = strcat('Figures/',groupName,'/0_Mean-vs-overallHist-',morphParamName16{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.png');
    fnameFig = strcat('Figures/',groupName,'/0_Mean-vs-overallHist-',morphParamName16{k}, '-', groupName,'-Ncells-', num2str(nCellsGroup),'-Nimages-', num2str(nImagesGroup),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc


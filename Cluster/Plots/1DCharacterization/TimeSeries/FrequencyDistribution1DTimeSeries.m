% FrequencyDistribution1DTimeSeries.m
% Program to calculate frequency distribution of different morphological
% identifiers for the Time Series data (Canada)
% 2023 ps UBC
 
% Full list of morphological parameters: 
% Area,Perimeter,Angle,Major,Minor,Height,Width,Mean,StdDev,Median,Skewness,Kurtosis,Feret,FeretAngle,MinFeret,FeretX,FeretY,Circularity,AR,Roundness,ConvexHullArea,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,ConvexPerimeter,Convexity,FibreLength,FibreWidth,Curl,CurlW
 


%% 1) Read .csv files with measured parameters from ImageJ or FIJI
tic

% File structure for data: 
% [main]/MorphologicalParameters/Canada/TimeSeries/[UID]/[.csv files]
% File structure for code: 
% [main]/Plots/1DCharacterization/TimeSeries/FrequencyDistribution1DTimeSeries.m
%dataPath = "../../../MorphologicalParameters/Canada/TimeSeries/";
 
% Different unique IDs or UIDs with time series data (e.g. SCD-001,
% SCT-003, etc.)
%uidName = "SCD-001"; %uidName is the unique ID for the data being analyzed
 
% List of relevant morphological parameters or identifiers from ImageJ or
% FIJI 
morphParamName = ["Area","Perimeter","Angle","Major","Minor","Height","Width","Mean","StdDev","Median","Skewness","Kurtosis","Feret","FeretAngle","MinFeret","FeretX","FeretY","Circularity","AR","Roundness","ConvexHullArea","Solidity","Eccentricity","ESF","ElongationFW","ElongationFmF","ConvexPerimeter","Convexity","FibreLength","FibreWidth","Curl","CurlW"];
morphParamAxis = ["Area (pixel^2)","Perimeter (pixel)","Angle (degree)","Major (pixel)","Minor (pixel)","Height (pixel)","Width (pixel)","Mean intensity","Standard deviation intensity","Median intensity","Skewness","Kurtosis","Feret (pixel)","FeretAngle (degree)","MinFeret (pixel)","FeretX (pixel)","FeretY (pixel)","Circularity","Aspect Ratio","Roundness","Convex hull area (pixel^2)","Solidity","Eccentricity","Elliptical shape factor","Elongation [Width/Feret]","Elongation [MinFeret/Feret]","Convex perimeter (pixel)","Convexity","Fibre length (pixel)","Fibre width (pixel)","Curl","CurlW"];

fprintf('Reading data from UID %s \n',uidName);
 
% csidFull is the ID of the coverslip
% List all the directories in UID, which include csidFull
csidFullL = dir(strcat(dataPath,uidName,'/*_*'));

% morphParam contains the information about the morphological parameters
% for each cell of all the images in the time series
% morphParam is a 3D cell of levels 1) CSID, 2) Image (as a function of
% time) and 3), morphological parameter
% morphParam{1,csidI}{1,imI}{1,mpI}

morphParam = cell(1, length(csidFullL));

% nImages is an array containing the number of images in each CSID
nImages = zeros(1,length(csidFullL));
% nCells is a cell with arrays containing the number of cells in each image
% in each CSID 
nCells = cell(1,length(csidFullL));

% dt is an array containing the time interval between images (dt) for each
% CSID (in seconds)
dt = zeros(1,length(csidFullL));

% t is a cell with arrays containing the time (based on dt) for each CSID
t = cell(1,length(csidFullL));

% Loop through all CSIDs for the UID
for csidI = 1:1:length(csidFullL)
    clear("mpList");
    mpList =  dir(strcat(dataPath, uidName, '/', csidFullL(csidI).name,'/*.csv')); % List of all .csv files with morphological parameters
    nImages(csidI) = length(mpList);
    fprintf('\t \t Number of files in %s = %d \n',csidFullL(csidI).name, length(mpList)); 
    
    % Note: mpList is not organized in ascending order, but rather and
    % order similar to: 0, 100, 101 ....109, 10, 11, 12 .....
    % The following is used to sort it in ascending order
    
    % Extract the number between "i_" and "_dt" and store it in a separate field
    pattern = 'i_(\d+)_dt';
    for i = 1:1:length(mpList)
        match = regexp(mpList(i).name, pattern, 'tokens', 'once');
        if ~isempty(match)
            mpList(i).number = str2double(match{1});
        else
            mpList(i).number = NaN;  % or handle the case when no match is found
            fprintf("Number not found");
        end
    end

    % Sort the struct based on the extracted numbers
    [~, sortedIndices] = sort([mpList.number]);

    % Reorder the struct elements based on the sorted indices
    mpList = mpList(sortedIndices);
    
    % The following code is used to extract the dt or time interval
    % Extract the number (e.g. "5.0") between "_dt_" and "seconds"
    pattern2 = '_dt_(\d+\.\d+)seconds';
    match2 = regexp(mpList(1).name, pattern2, 'tokens', 'once');
    
    % Check if a match was found and extract the number
    if ~isempty(match2)
        dt(csidI) = str2double(match2{1});
        fprintf('\t \t Time interval between images (dt) = %d seconds \n', dt(csidI));
    else
        disp('No match found.');
    end
    
    % Read the csv files containing the morphological parameters
    nCells{1,csidI} = zeros(1,length(mpList));  % Initialize array with size of the total number of files
    t{1,csidI} = zeros(1,length(mpList));       % Initialize array with size of the total number of files
    
    % Loop through all the files containing morphological parameter
    % information for each image
    for imI = 1:1:length(mpList)
        t{1,csidI}(imI) = dt(csidI)*mpList(imI).number;
        
        warning('off','MATLAB:table:ModifiedAndSavedVarnames');
        fileName = strcat(dataPath, uidName, '/', csidFullL(csidI).name,'/', mpList(imI).name);
        data = readtable(fileName);
        
        nCells{1,csidI}(imI) = length(data.Var1);
        
        %Iterate through all morphological parameters (columns) within each file
        for mpI = 1:1:length(morphParamName)
            % morphParam contains the information about the morphological parameters
            % for each cell of all the images in the time series
            % morphParam is a 3D cell of levels 1) CSID, 2) Image (as a function of
            % time) and 3), morphological parameter
            % morphParam{1,csidI}{1,imN}{1,mpI}
            
            morphParam{1,csidI}{1,imI}{1,mpI} = eval(strcat('data.',morphParamName(mpI))); %Save morphological parameter value (array) into a 3D cell
            
        end
    end
end

toc


%% 2) Find minimum and maximum values of all morphological parameters
tic
% Min and Max store minimum and maximum values for all morphological
% parameters
Min = zeros(1,length(morphParamName));
Max = zeros(1,length(morphParamName));

for k = 1:1:length(morphParamName)
    clear("MinG"); clear("MaxG");
    MinG = zeros(1,length(csidFullL));
    MaxG = zeros(1,length(csidFullL));
    
    for i = 1:1:length(csidFullL)
        
        clear("MinT"); clear("MaxT");
        %MinT and MaxT store minimum and maximum values of each file for
        %the specific morphological parameter
        MinT = zeros(1,length(morphParam{1,i}));
        MaxT = zeros(1,length(morphParam{1,i}));
        
        for j = 1:1:length(morphParam{1,i})
            MinT(j) = min(morphParam{1,i}{1,j}{1,k});
            MaxT(j) = max(morphParam{1,i}{1,j}{1,k});       
        end
        
        MinG(i) = min(MinT);
        MaxG(i) = max(MaxT);
        
    end
    
    Min(k) = min(MinG);
    Max(k) = max(MaxG);
    
    
    
end

%{
% Note: Min and Max are chosen based on values of all morphological
% parameters from Canada data
% Area,Perimeter,Angle,Major,Minor,Height,Width,Mean,StdDev,Median,Skewness,Kurtosis,Feret,FeretAngle,MinFeret,FeretX,FeretY,Circularity,AR,Roundness,ConvexHullArea,Solidity,Eccentricity,ESF,ElongationFW,ElongationFmF,ConvexPerimeter,Convexity,FibreLength,FibreWidth,Curl,CurlW
Min = [0 0 0 0 0 0 0 0 0 0 -6 -5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
Max = [2000 250 180 120 45 115 115 250 105 255 6 30 115 180 55 3000 3000 1 30 1 2500 1 1 1 1 1 250 1 40 40 1.6 1.6];
%}

% Make a variable edges with equally spaced edges for histocount for all
% morphological parameters using max and min values above
nbins = 30; % 30 is chosen because histData is per image (~3000-7000 cells/image)
 
edges = cell(1, length(morphParamName));
midBin = cell(1, length(morphParamName));
for k = 1:1:length(morphParamName)
    % edges includes edges for bins
    edges{k} = linspace(Min(k),Max(k),(nbins+1)); 
    
    % midBin is the middle value of each bin (used for plotting)
    midBin{k} = edges{k}(2:end) - (edges{k}(2)-edges{k}(1))/2;
    
end
toc

%% 3) Calculate frequency distribution and emperical cumulative distribution
% Save histocount data for each image and each morphological identifier
% histData is a 3D cell of levels 1) CSID, 2) Image (as a function of
% time) and 3), morphological parameter
% histData{1,csidI}{1,imI}{1,mpI}

% ecdfY contains emperical cumulative distribution for all data
% ecdfX contains X axis for emperical cumulative distribution for all data
% ecdfY{1,csidI}{1,imI}{1,mpI}
% ecdfX{1,csidI}{1,imI}{1,mpI}
tic
 
histData = cell(1, length(csidFullL));
ecdfY = cell(1, length(csidFullL));
ecdfX = cell(1, length(csidFullL));
% Iterate through all files within groupName
fprintf('Calculating histcounts for UID %s.\n',uidName);

for csidI = 1:1:length(csidFullL)    
    
    for imI = 1:1:length(morphParam{1,csidI})
       
        for mpI = 1:1:length(morphParamName)
            [histData{1,csidI}{1,imI}{1,mpI}, edgesTemp] = histcounts(morphParam{1,csidI}{1,imI}{1,mpI}, edges{mpI}, 'Normalization','probability');
            [ecdfY{1,csidI}{1,imI}{1,mpI}, ecdfX{1,csidI}{1,imI}{1,mpI}] = ecdf(morphParam{1,csidI}{1,imI}{1,mpI});
        end

    end

end
toc
%% 3) Save variables for morphParam, general information, and histData
% Note: These variables are per UID, and store information for all CSID

% Create a folder with the uidpName if one does not exist
if not(isfolder(uidName))
    mkdir(uidName)
end
 
tic 
fprintf("Saving general information variables.... \n");
save(strcat(uidName,"/generalInfo.mat"), "uidName", "morphParamName", "morphParamAxis", "nImages","nCells", "dt", "t", "csidFullL", '-v7.3');
 
fprintf("Saving variables for morphParam.... \n");
save(strcat(uidName, "/morphParamData.mat"),"morphParam", '-v7.3');

fprintf("Saving histData.... \n");
save(strcat(uidName,"/histData.mat"), "edges", "midBin", "histData", "ecdfY", "ecdfX", '-v7.3');

toc

%% 4) Figures

for csidI = 1:1:length(histData)
    if not(isfolder(strcat(uidName,"/",csidFullL(csidI).name,"/Figures")))
        mkdir(strcat(uidName,"/",csidFullL(csidI).name,"/Figures"))
    end
end

%% 5a) Plot 1D frequency distribution - individually for morphological parameters

tic 
% Iterate through all coverslip IDs
for csidI = 1:1:length(histData)
    fprintf('\t Creating and saving 1D frequency distribution plots for %s \n',csidFullL(csidI).name);
    % Itearate through all morphological parameters
    for mpI = 1:1:length(morphParamName)
        figure;
        hold on;
        % Create colorbar
        c = colorbar;
        c.Label.String = 'Time (s)';
        cmap2 = colormap(flipud(parula));
        cmap = colormap(cmap2(round(length(cmap2)/6):length(cmap2),:)); % Omit part of the colormap that is lighter

        mint = t{1,csidI}(1); maxt = t{1,csidI}(end); % Minimum and maximum time
        caxis([mint maxt]);
        
        % Iterate through all images or files and plot the temporal
        % variation of the morphological parameter
        for imI = 1:1:length(histData{1,csidI})
            iCmap = round(interp1([mint maxt],[1 length(cmap)],t{1,csidI}(imI)));

            %set the line color using the colormap
            color = cmap(iCmap,:);
            plot(midBin{mpI}, histData{1,csidI}{1,imI}{1,mpI}, 'color',color, 'LineWidth', 0.15);

        end
        hold off
        box on
        xlabel(morphParamAxis{mpI}); ylabel("Normalized frequency");
        set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
        set(findall(gcf,'-property','FontSize'),'FontSize',8);
        fnamePNG = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/1_TimeSeries_', uidName, '-', morphParamName{mpI}, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)), '-Nbins-', num2str(nbins),'.png');
        fnameFig = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/1_TimeSeries_', uidName, '-', morphParamName{mpI}, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)), '-Nbins-', num2str(nbins),'.fig');
        print(gcf,fnamePNG,'-dpng','-r600');
        saveas(gcf,fnameFig);
        close;
    end
    
end
toc

%% 5b) Plot and save the number of cells per image (over time)
tic
for csidI = 1:1:length(csidFullL)
    fprintf('\t Creating and saving # of images plot for %s \n',csidFullL(csidI).name);
    figure; 
    plot(t{1,csidI},nCells{1, csidI});
    xlabel('Time (s)'); ylabel('Number of cells');
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8);
    fnamePNG = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/0_TimeSeries_NumberofCells-', uidName, '-', '.png');
    fnameFig = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/0_TimeSeries_NumberofCells-', uidName, '-', '.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc

%% 5c) Plot empirical cumulative distribution 

tic
% Iterate through all coverslip IDs
for csidI = 1:1:length(histData)
    fprintf('\t Creating and saving empirical cumulative distribution plots for %s \n',csidFullL(csidI).name);
    % Itearate through all morphological parameters
    for mpI = 1:1:length(morphParamName)
        figure;
        hold on;
        %create colorbar
        c = colorbar;
        c.Label.String = 'Time (s)';
        cmap2 = colormap(flipud(parula));
        cmap = colormap(cmap2(round(length(cmap2)/6):length(cmap2),:)); % Omit part of the colormap that is lighter

        mint = t{1,csidI}(1); maxt = t{1,csidI}(end);   % Minimum and maximum time (for plotting)
        caxis([mint maxt]);
        % Iterate through all images or files and plot the temporal
        % variation of the morphological parameter
        for imI = 1:1:length(histData{1,csidI})
            iCmap = round(interp1([mint maxt],[1 length(cmap)],t{1,csidI}(imI)));

            %set the line color using the colormap
            color = cmap(iCmap,:);

            %[f, x] = ecdf(morphParam{1,csidI}{1,imI}{1,mpI});
            %plot(x, f, 'color',color, 'LineWidth', 0.15);
            plot(ecdfX{1,csidI}{1,imI}{1,mpI}, ecdfY{1,csidI}{1,imI}{1,mpI}, 'color',color, 'LineWidth', 0.15);

        end
        hold off; box on;
        xlabel(morphParamAxis{mpI}); ylabel("Empirical cumulative distribution");
        set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 2.163], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 2.163]);
        set(findall(gcf,'-property','FontSize'),'FontSize',8);

        fnamePNG = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/2_TimeSeries_ECDF_', uidName, '-', morphParamName{mpI}, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)),'.png');
        fnameFig = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/2_TimeSeries__ECDF_', uidName, '-', morphParamName{mpI}, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)),'.fig');
        print(gcf,fnamePNG,'-dpng','-r600');
        saveas(gcf,fnameFig);
        close;
    end
    
end
toc

%% 5d) Plot frequency distribution - ALL at once
% tiledlayout version (Available from 2019 matlab version onwards)

tic
% Iterate through all coverslip IDs
for csidI = 1:1:length(histData)
    fprintf('\t Creating and saving combined plot for 1D frequency distribution for %s \n',csidFullL(csidI).name);
    figure;
    tL = tiledlayout('flow');
    %create colorbar
    c = colorbar;
    c.Label.String = 'Time (s)';
    cmap2 = colormap(flipud(parula));
    cmap = colormap(cmap2(round(length(cmap2)/6):length(cmap2),:)); % Omit part of the colormap that is lighter
    
    mint = t{1,csidI}(1); maxt = t{1,csidI}(end);   % Minimum and maximum time (for plotting)
    caxis([mint maxt]);
    % Itearate through all morphological parameters
    for mpI = 1:1:length(morphParamName)
        nexttile
        hold on;
        
        for imI = 1:1:length(histData{1,csidI})
            iCmap = round(interp1([mint maxt],[1 length(cmap)],t{1,csidI}(imI)));

            %set the line color using the colormap
            color = cmap(iCmap,:);

            plot(midBin{mpI}, histData{1,csidI}{1,imI}{1,mpI}, 'color',color, 'LineWidth', 0.15);

        end
        hold off; box on;
        xlabel(morphParamAxis{mpI}); 
        
    end
    tL.TileSpacing = 'tight';
    tL.Padding = 'compact';
    ylabel(tL, 'Frequency distribution');
    c.Layout.Tile = 'east'; % Place color bar to the east of the tiled layout
    
    fnamePNG = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/3_TimeSeries_FreqDist_Combined', uidName, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)),'.png');
    fnameFig = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/3_TimeSeries_FreqDist_Combined', uidName, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc

%% 5e) Plot empirical cumulative distribution - ALL at once
% tiledlayout version (Available from 2019 matlab version onwards)

tic
% Iterate through all coverslip IDs
for csidI = 1:1:length(histData)
    fprintf('\t Creating and saving combined plot for empirical cumulative distribution for %s \n',csidFullL(csidI).name);
    figure;
    tL = tiledlayout('flow');
    %create colorbar
    c = colorbar;
    c.Label.String = 'Time (s)';
    cmap2 = colormap(flipud(parula));
    cmap = colormap(cmap2(round(length(cmap2)/6):length(cmap2),:)); % Omit part of the colormap that is lighter
    
    mint = t{1,csidI}(1); maxt = t{1,csidI}(end);   % Minimum and maximum time (for plotting)
    caxis([mint maxt]);
    % Itearate through all morphological parameters
    for mpI = 1:1:length(morphParamName)
        nexttile
        hold on;
        
        for imI = 1:1:length(histData{1,csidI})
            iCmap = round(interp1([mint maxt],[1 length(cmap)],t{1,csidI}(imI)));

            %set the line color using the colormap
            color = cmap(iCmap,:);

            %[f, x] = ecdf(morphParam{1,csidI}{1,imI}{1,mpI});
            %plot(x, f, 'color',color, 'LineWidth', 0.15);
            plot(ecdfX{1,csidI}{1,imI}{1,mpI}, ecdfY{1,csidI}{1,imI}{1,mpI}, 'color',color, 'LineWidth', 0.15);

        end
        hold off; box on;
        xlabel(morphParamAxis{mpI}); 
        
    end
    tL.TileSpacing = 'tight';
    tL.Padding = 'compact';
    ylabel(tL, 'Empirical cumulative distribution');
    c.Layout.Tile = 'east'; % Place color bar to the east of the tiled layout
    
    fnamePNG = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/4_TimeSeries_ECDF_Combined', uidName, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)),'.png');
    fnameFig = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/4_TimeSeries_ECDF_Combined', uidName, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
toc

%{
%% Plot emperical cumulative distribution - ALL at once
% subplot version

% Iterate through all coverslip IDs
for csidI = 1:1:length(histData)
    figure;
    %create colorbar
    c = colorbar;
    c.Label.String = 'Time (s)';
    cmap2 = colormap(flipud(parula));
    cmap = colormap(cmap2(round(length(cmap2)/6):length(cmap2),:)); % Omit part of the colormap that is lighter
    
    mint = t{1,csidI}(1); maxt = t{1,csidI}(end);   % Minimum and maximum time (for plotting)
    caxis([mint maxt]);
    % Itearate through all morphological parameters
    for mpI = 1:1:length(morphParamName)
        subplot(8,4,mpI)
        hold on;
        
        for imI = 1:1:length(histData{1,csidI})
            iCmap = round(interp1([mint maxt],[1 length(cmap)],t{1,csidI}(imI)));

            %set the line color using the colormap
            color = cmap(iCmap,:);

            %[f, x] = ecdf(morphParam{1,csidI}{1,imI}{1,mpI});
            %plot(x, f, 'color',color, 'LineWidth', 0.15);
            plot(ecdfX{1,csidI}{1,imI}{1,mpI}, ecdfY{1,csidI}{1,imI}{1,mpI}, 'color',color, 'LineWidth', 0.15);

        end
        hold off; box on;
        xlabel(morphParamAxis{mpI}); 
        % Include a y label only on the left most column
        if mod(mpI,4) == 1
            ylabel("ECDF");
        end
        
    end
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5*5, 2.163*5], 'PaperUnits', 'Inches', 'PaperSize', [3.5*5, 2.163*5]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8);
     
    fnamePNG = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/3_TimeSeries_ECDF_Combined', uidName, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)),'.png');
    fnameFig = strcat(uidName,'/',csidFullL(csidI).name,'/Figures/3_TimeSeries__ECDF_Combined', uidName, '-dt-', num2str(dt(csidI)),'-seconds-Nt-', num2str(nImages(csidI)),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFig);
    close;
end
%}
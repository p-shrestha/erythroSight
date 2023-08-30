% statDiff_cnComb_4groups_30b40mp_coverslipL.m
% Program to measure and plot statistical diffenreces between frequency
% distributions of morphological parameters (30 bin 40 morph Param). 
% Coverslip level 
% ps UBC 2023

%% 1) Read the frequency distribution data (30 bin) or histData

% Note: Here, the histDataCS (or coverslip-wise) histcount data is selected
% (For finer analysis, image-wise data or groups of image-wise data can be
% selected, for coarser analysis, donor-wise data can be selected
tic 
fprintf('Reading data - generalInfo.mat and histData.mat \n');

group = {'AA', 'ABeta','AS','SCD'};

for i = 1:1:length(group)
    % Run similar command for all groups to read generalInfo.mat and
    % histData.mat
    % Example: AA.info = load('../../../1DCharacterization/t02hrs/cnCombined/30b40mp/AA/generalInfo.mat');
    % Example: AA.hD = load('../../../1DCharacterization/t02hrs/cnCombined/30b40mp/AA/histData.mat');
    eval(strcat(group{i},".info = load('../../../1DCharacterization/t02hrs/cnCombined/30b40mp/",group{i},"/generalInfo.mat');"));
    eval(strcat(group{i},".hD = load('../../../1DCharacterization/t02hrs/cnCombined/30b40mp/",group{i},"/histData.mat');"));

end

toc


%% 2) Concatenate data for Canada and Nepal data so that all coverslips for each group are together

tic

%group = ["AA", "ABeta", "AS", "SCD"]; % Combined data from Nepal and Canada; SCD = SBeta + SS
fprintf('Concatenating data for Canada and Nepal data so that all coverslips for each group are together \n');

% histDataCS30b40mp{i}{j}{k} includes data for i = donor, j = coverslip, and
% k = morphological parameter. This needs to be adjusted so that all the
% coverslips for a particular group are in the same group in histData
% histData{i}{k} includes coverslip level data of i = groups, k =
% morphological parameters

% Concatenate data such that all coverslips for one group is stored in the
% i level of histData{i}{k}
for i = 1:1:length(group)
    % Store the first set of all coverslips from the first donor
    % Similar to command: histData{1} = AA.hD.histDataCS30b40mp;
    eval(strcat('histData{',num2str(i),'}=',group{i},'.hD.histDataCS30b40mp{1};'));

    dNumber = eval(strcat('length(',group{i},'.hD.histDataCS30b40mp', ')'));
    for iD = 2:1:dNumber
        % Concatenate the coverslips for the rest of the donors
        % Similar to histData{1} = cat(2, histData{1},  AA.hD.histDataCS30b40mp{iD});

        eval(strcat('histData{',num2str(i),'}=cat(2,histData{',num2str(i),'},',group{i},'.hD.histDataCS30b40mp{',num2str(iD),'});'))

    end
end

toc

%% 3) Calculate the mean and percentiles for all the data
% These are averaged coverslip-wise 

tic
morphParamName = AA.info.morphParamName;

% meanHist{i}{k} has mean frequency distribution for group i and
% morphologial parameter k
meanHist = cell(1, length(group));
perc25 = cell(1, length(group));
perc75 = cell(1, length(group));
perc5 = cell(1, length(group));
perc95 = cell(1, length(group));
stdHist = cell(1, length(group));

for i = 1:1:length(group)
    fprintf('%d out of %d) Calculating mean, median, percentiles and std. dev. for group %s.\n',i,length(group),group{i});
    for k = 1:1:length(morphParamName)
        
        for p = 1:1:length(histData{1}{1}{1}) % Iterate through all bins
            clear('jTemp');
            jTemp = zeros(1,length(histData{i})); 
            
            for j = 1:1:length(histData{i}) % Iterate through all files or coverslip frequencies
                jTemp(j) = histData{i}{j}{k}(p);
                
            end
            meanHist{i}{k}(p) = mean(jTemp);
            medianHist{i}{k}(p) = median(jTemp);
            perc25{i}{k}(p) = prctile(jTemp,25);
            perc75{i}{k}(p) = prctile(jTemp,75);
            perc5{i}{k}(p) = prctile(jTemp,5);
            perc95{i}{k}(p) = prctile(jTemp,95);
            stdHist{i}{k}(p) = std(jTemp);
            
        end
    
    end
end

toc


% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%    One-way ANOVA (comparing all 4 groups together)
%    and Multiple Comparison test (pairwise comparison of all 4 groups)
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %



% ----------------------------------------------------------------------- %
%                   All morphological parameters
% ----------------------------------------------------------------------- %


%% 4) One-way ANOVA and multiple comparison test for 4 groups (multiple implementation - optimization)
% Note: Multiple values (morphological parameter) and
% all threshold or MPt values to find combination with minimum p values
% (iterate through all midBin{k} values
% Code for single implementation also provided later

tic

if not(isfolder("/Figures"))
    mkdir("Figures");
end

midBin = AA.hD.midBin;
morphParamLegends = AA.info.morphParamName;
morphParamName = AA.info.morphParamName;
morphParamAxis = AA.info.morphParamAxis;

pCombM_Allk = cell(1,length(morphParamLegends));
pAvgMin_Allk = zeros(1,length(morphParamLegends));
pAvgMinInd_Allk = zeros(1,length(morphParamLegends));
pIndMinInd_Allk = zeros(1,length(morphParamLegends));

% The following include a matrix for all morphological parameters (columns)
% and all combinations or pairs (rows)
% pCombM_oneMatrix contains the p values for the combined minimum
% pIndM_oneMatrix contains the p values for the individual minimum
pCombM_oneMatrix = zeros(6,length(morphParamLegends));
pIndM_oneMatrix = zeros(6,length(morphParamLegends));

%k = 23; % Eccentricity
%MPt = 0.5; % Threshold for morphological parameter

for k = 1:1:length(morphParamLegends)
    fprintf('Calculating statistical differences for morphological parameter %s \n', morphParamName{k})
    pAvgL = zeros(1,length(midBin{k})); % Array storing average p values
    pCombM = zeros(6,length(midBin{k}));    % Matrix storing all p values
    for MPtInd = 1:1:length(midBin{k})

        MPt = midBin{k}(MPtInd);    % Assign value for MPt

        % Find value of the morphological parameter at the threshold MPt for all
        % files or coverslips
        fMPt = cell(1,length(group));
        for i = 1:1:length(group)
            fMPt{i} = zeros(1,length(histData{1,i}));

            for j = 1:1:length(histData{1,i})
                [ d, ix ] = min( abs( midBin{k}- MPt ) );
                fMPt{i}(j) = histData{1,i}{1,j}{1,k}(ix);
            end
        end

        % Combine all 4 groups
        Xbox = [fMPt{1} fMPt{2} fMPt{3} fMPt{4}];
        grp = [ones(size(fMPt{1})) 2.*ones(size(fMPt{2})) 3.*ones(size(fMPt{3})) 4.*ones(size(fMPt{4}))];

        % One-way ANOVA
        % https://www.mathworks.com/help/stats/one-way-anova.html
        [p,tbl,stats] = anova1(Xbox, grp);

        % Kruskal-Wallis test - non-parametric version
        % https://www.mathworks.com/help/stats/kruskalwallis.html
        %[p,tbl,stats] = kruskalwallis(Xbox, grp);

        close; close; % Close the two plots from the One-way ANOVA

        % Multiple compare
        figure;
        results = multcompare(stats, "CriticalValueType","scheffe");
        close; % Close plot for multiple comparison test

        pCombM(:,MPtInd) = results(:,6); % Stores the p values for all the combinations
        %pAvgL(MPtInd) = mean(results(:,6)); % Stores the average of p values for all combinations
        pAvgL(MPtInd) = geomean(results(:,6)); % Stores the geometric mean of p values for all combinations

    end

    pCombM(isnan(pCombM))=1; % Replace all NaN with 1 (NaN occurs when comparisons cannot be made; 1 indicates no difference)



    [pAvgMin, pAvgMinInd] = min(pAvgL); % Find the minimum averaged p value and index
    MPt = midBin{k}(pAvgMinInd);    % Find the corresponding MPt with the minimum average p value

    pAvgMin_Allk(k) = pAvgMin;
    pAvgMinInd_Allk(k) = pAvgMinInd;

    fMPt = cell(1,length(group));
    for i = 1:1:length(group)
        fMPt{i} = zeros(1,length(histData{1,i}));

        for j = 1:1:length(histData{1,i})
            [ d, ix ] = min( abs( midBin{k}- MPt ) );
            fMPt{i}(j) = histData{1,i}{1,j}{1,k}(ix);
        end
    end


    % Plotting

    % Combine all 4 groups
    Xbox = [fMPt{1} fMPt{2} fMPt{3} fMPt{4}];
    grp = [ones(size(fMPt{1})) 2.*ones(size(fMPt{2})) 3.*ones(size(fMPt{3})) 4.*ones(size(fMPt{4}))];
    results = multcompare(stats, "CriticalValueType","scheffe");    % Run multiple comparison test
    close;  % Close the default plot 
    pComb = results(:,6); % Stores the p values for all the combinations
    pMinAvg = mean(results(:,6)); % Stores the average of p values for all combinations


    % --- Figure 1 --------------------

    figure
    tL = tiledlayout('flow', 'Padding','tight', 'TileSpacing','compact');
    nexttile([1 2]);
    %nexttile;
    imagesc(pCombM);
    %heatmap(pCombM);

    set(gca,'ColorScale','log');
    yticks([1,2,3,4,5,6]); yticklabels({'AA-ABeta', 'AA-AS','AA-SCD', 'ABeta-AS', 'ABeta-SCD', 'AS-SCD'});
    cmap1 = colormap(sky);
    axis image; % so that the boxes in heatmaps are squares


    pCombM(pCombM == 0) = 1e-300;   % For any p value = 0, make it equate 1*10^-350, which is less than lowest recorded p (~10^-300)
    nM =floor( log10(min(min(pCombM))));
    nMe = 2*floor(nM/2) + 10;
    if nMe < -300
        nMe = -300;  % If limit is less than -300, the labelling of limits in the colorbar is in weird notation, so instead show 10^-300 location instead of max
    end
    nMid = (-2 + nMe)/2;
    colormap(flipud(cmap1));
    clim([min(min(pCombM)) 0.05]);
    %clim([min(min(pCombM)) max(max(pCombM))]);

    cbar = colorbar;
    cbar.Label.String = "p-value";

    % Ticks for colorbar - select 3 ticks: max at 10^-2; min at 10^-nMe, where
    % nMe is 10 more than (so that label is not hidden) than exponent nM
    % rounded to nearest even number; and one between the two

    

    % Set three ticks in the colorbar, similar to the command:
    % cbar.Ticks = [1e-180, 1e-91, 1e-2];
    % cbar.Ticks = [nMe, nMid, 1e-2];
    eval(strcat('cbar.Ticks = [1e',num2str(nMe), ',1e',num2str(nMid), ',1e-2];' ));

    % Set a threshold for p value and color all above this threshold light
    % orange (so all combinations that are not significantly different are
    % colored orange).

    pTh = 0.05; % 5% threshold
    [xM, yM] = size(pCombM);
    for ii = 1:1:xM
        for jj = 1:1:yM
            if pCombM(ii,jj) >= 0.05
                rectangle('Position',[jj-0.5 ii-0.5 1 1],'FaceColor',[253/255, 208/255, 162/255],'LineStyle','none'); % Create a rectangle
            end
        end
    end

    % Highlight overall minimum at index pAvgMinInd 116/255, 196/255, 118/255
    %rectangle('Position',[pAvgMinInd-0.5 1-0.5 1 6],'EdgeColor',[49/255, 163/255, 84/255]);%, 'Curvature',0.5); % Create a rectangle
    rectangle('Position',[pAvgMinInd-0.5 1-0.5 1 6],'EdgeColor',[180/255, 180/255, 180/255]);



    pCombDiffBin = zeros(6,1); % Stores minimum of each row or combination (can be different bin)
    pCombDiff_binIndex = zeros(6,1);
    % Highlight the minimum for each row (between each combination)
    for ii = 1:1:xM
        [xMMin, xMi] = min(pCombM(ii,:));
        rectangle('Position',[xMi-0.5 ii-0.5 1 1],'EdgeColor',[0/255, 90/255, 50/255], 'Curvature',0.5); % Create a rectangle
        pCombDiffBin(ii) = xMMin;
        pCombDiff_binIndex(ii) = xMi;

        pIndMinInd_Allk(ii,k) = xMi;

        pIndM_oneMatrix(ii,k) = xMMin; % Store value of the individual minimum in the matrix storing all values from all morphological parameters
    end
    dRound = 2;
    xticklabels(round([midBin{k}(5), midBin{k}(10), midBin{k}(15), midBin{k}(20), midBin{k}(25), midBin{k}(30)],dRound));
    xlabel(morphParamAxis{k});

    % --- Figure 2 --------------------

    colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [231/255 41/255 138/255], [166/255 118/255 29/255],[35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
    %figure;
    nexttile
    hold on;
    for i = 1:1:length(group)
        plot(midBin{k},meanHist{i}{k}, 'Color',colors{i});
        xlabel(morphParamAxis{k});
        ylabel('Normalized frequency');
    end

    % Percentiles plotted in a different loop to keep ordering of legends
    for i = 1:1:length(group)
        x2 = [midBin{k}, fliplr(midBin{k})];
        inBetween = [perc25{i}{k}, fliplr(perc75{i}{k})];
        fill(x2, inBetween, colors{i},'FaceAlpha',0.1, 'LineStyle','none');
    end
    %legend('AA', 'ABeta','AS','SCD', 'Location', 'northoutside','Orientation','horizontal', 'box', 'off');
    hold off
    % Highlight overall minimum at index pAvgMinInd
    %dx = (midBin{k}(ix) - midBin{k}(ix-1))/2;
    dx = (midBin{k}(2) - midBin{k}(1))/2;
    yL = ylim;
    rectangle('Position',[midBin{k}(ix-1)+dx yL(1) dx*2 yL(2)-yL(1)],'EdgeColor',[180/255, 180/255, 180/255]);%, 'Curvature',0.5); % Create a rectangle
    box on;
    ylim([yL(1), yL(2)]);

    % --- Figure 3 --------------------


    nexttile;
    % For plotting box plot
    %figure;
    %boxplot(Xbox,grp,'Notch', 'on');
    boxplot(Xbox,grp,'Notch', 'on', 'Colors',[27/255 158/255 119/255; 217/255 95/255 2/255; 117/255 112/255 179/255; 231/255 41/255 138/255],'Symbol','.' ); %'PlotStyle','compact'
    xticklabels({'AA', 'ABeta', 'AS', 'SCD'});
    %ylabel(strcat('Normalized frequency at MP_t =',{' '},num2str(MPt)));
    ylabel('Normalized frequency');
    %title(morphParamLegends{k});

    % One-way ANOVA
    % https://www.mathworks.com/help/stats/one-way-anova.html
    [p,tbl,stats] = anova1(Xbox, grp);

    % Kruskal-Wallis test - non-parametric version
    % https://www.mathworks.com/help/stats/kruskalwallis.html
    %[p,tbl,stats] = kruskalwallis(Xbox, grp);

    close; close; % Close the two plots from the One-way ANOVA

    % --- Figure 4 --------------------

    % Plot mean and standard deviation
    %figure;
    nexttile;
    hold on;
    eb1 = errorbar(1,mean(fMPt{1}),std(fMPt{1}), 'o', 'MarkerSize', 3, 'MarkerEdgeColor',[19/255 129/255 96/255], 'MarkerFaceColor',[27/255 158/255 119/255]);
    eb2 = errorbar(2,mean(fMPt{2}),std(fMPt{2}), 'o', 'MarkerSize', 3, 'MarkerEdgeColor',[179/255 78/255 1/255],'MarkerFaceColor',[217/255 95/255 2/255]);
    eb3 = errorbar(3,mean(fMPt{3}),std(fMPt{3}), 'o', 'MarkerSize', 3, 'MarkerEdgeColor',[94/255 90/255 142/255],'MarkerFaceColor',[117/255 112/255 179/255]);
    eb4 = errorbar(4,mean(fMPt{4}),std(fMPt{4}), 'o', 'MarkerSize', 3, 'MarkerEdgeColor',[186/255 33/255 110/255],'MarkerFaceColor',[231/255 41/255 138/255]);
    hold off;
    eb1.Color = [19/255 129/255 96/255];
    eb2.Color = [179/255 78/255 1/255];
    eb3.Color = [94/255 90/255 142/255];
    eb4.Color = [186/255 33/255 110/255];
    xticks([1 2 3 4]);
    xticklabels({'AA', 'ABeta', 'AS','SCD'});
    xlim([0.5 4.5]);
    ylabel(strcat('Normalized frequency'));


    %{
    % Multiple compare
    results = multcompare(stats, "CriticalValueType","scheffe");
    % Edit/label plots
    %xlabel(strcat('Mean of frequency at MP_t =',{' '},num2str(MPt)));
    xlabel(strcat('Mean normalized frequency'));
    title('');
    %ylabel('Label Y');
    yticklabels({'SCD', 'AS', 'ABeta','AA'});
    set ( gca, 'ydir', 'reverse' , 'xdir', 'reverse');
    ax2 = gca; set(ax2, 'YAxisLocation', 'right');
    camroll(270);
    pComb = results(:,6); % Stores the p values for all the combinations
    pMinAvg = mean(results(:,6)); % Stores the average of p values for all combinations
    %}

    % --- Figure 5 --------------------

    % Visualize p value
    %figure;
    nexttile;
    semilogy(pCombM(:,pAvgMinInd),'.', 'MarkerSize', 12);
    hold on; plot([0 7], [0.05 0.05], '--');
    semilogy([1 2 3 4 5 6], pCombDiffBin,'o', 'Color', [0/255, 90/255, 50/255],'MarkerSize', 5);
    hold off;
    xticks([1 2 3 4 5 6]); xlim([0 7])
    xticklabels({'AA-ABeta', 'AA-AS','AA-SCD', 'ABeta-AS', 'ABeta-SCD', 'AS-SCD'});
    %ylabel(strcat('p-value at MP_t = ',{' '},num2str(MPt)));
    ylabel('p-value');
    %title(strcat(morphParamLegends{k}, {' '},'at MP_t =',num2str(MPt)));

    % Save
    %tL.Units = 'centimeters';
    %t.OuterPosition = [0.635 0.635 17 34];
    %exportgraphics(tL,'test.png','Resolution',300)

    %set(gcf, 'Units', 'centimeters', 'Position', [4, 4, 17, 10.51], 'PaperUnits', 'Inches', 'PaperSize', [17, 17]);
    set(findall(gcf,'-property','FontSize'),'FontSize',7)
    fnamePNG = strcat('Figures/1_',num2str(k),'_statDiff_',morphParamLegends{k}, '_MPt_',num2str(MPt),'.png');
    fnameFIG = strcat('Figures/1_',num2str(k),'_statDiff_',morphParamLegends{k},'_MPt_',num2str(MPt),'.fig');
    print(gcf,fnamePNG,'-dpng','-r600');
    saveas(gcf,fnameFIG);
    close;

    % Store variables in total variable - morphological parameter wise
    pCombM_Allk{k} = pCombM;

    pCombM_oneMatrix(:,k) = pCombM(:,pAvgMinInd);


end
toc

%% 5) Plot combined matrices
% Plot the combined matrix with overall minimum p-value (with minimum
% geometric mean)
tic
fprintf('Plotting combined heatmaps \n');

figure
tL = tiledlayout('flow', 'Padding','tight', 'TileSpacing','compact');
nexttile([1 2]);
%nexttile;
imagesc(pCombM_oneMatrix);
%heatmap(pCombM_oneMatrix);

set(gca,'ColorScale','log');
yticks([1,2,3,4,5,6]); yticklabels({'AA-ABeta', 'AA-AS','AA-SCD', 'ABeta-AS', 'ABeta-SCD', 'AS-SCD'});
cmap1 = colormap(sky);
axis image; % so that the boxes in heatmaps are squares

colormap(flipud(cmap1));
clim([min(min(pCombM_oneMatrix)) 0.05]);
%clim([min(min(pCombM_oneMatrix)) max(max(pCombM_oneMatrix))]);

cbar = colorbar;
cbar.Label.String = "p-value";

% Ticks for colorbar - select 3 ticks: max at 10^-2; min at 10^-nMe, where
% nMe is 10 more than (so that label is not hidden) than exponent nM
% rounded to nearest even number; and one between the two

nM =floor( log10(min(min(pCombM_oneMatrix))));
nMe = 2*floor(nM/2) + 10;
nMid = (-2 + nMe)/2;
if nMe < -300
    nMe = -300;  % If limit is less than -300, the labelling of limits in the colorbar is in weird notation, so instead show 10^-300 location instead of max
end

% Set three ticks in the colorbar, similar to the command:
% cbar.Ticks = [1e-180, 1e-91, 1e-2];
eval(strcat('cbar.Ticks = [1e',num2str(nMe), ',1e',num2str(nMid), ',1e-2];' ));

% Set a threshold for p value and color all above this threshold light
% orange (so all combinations that are not significantly different are
% colored orange).

pTh = 0.05; % 5% threshold
[xM, yM] = size(pCombM_oneMatrix);
for ii = 1:1:xM
    for jj = 1:1:yM
        if pCombM_oneMatrix(ii,jj) >= 0.05
            rectangle('Position',[jj-0.5 ii-0.5 1 1],'FaceColor',[253/255, 208/255, 162/255],'LineStyle','none'); % Create a rectangle
        end
    end
end

% Highlight overall minimum at index pAvgMinInd 116/255, 196/255, 118/255
%rectangle('Position',[pAvgMinInd-0.5 1-0.5 1 6],'EdgeColor',[49/255, 163/255, 84/255]);%, 'Curvature',0.5); % Create a rectangle
rectangle('Position',[pAvgMinInd-0.5 1-0.5 1 6],'EdgeColor',[180/255, 180/255, 180/255]);

% Highlight the minimum for each row (between each combination)
for ii = 1:1:xM
    [xMMin, xMi] = min(pCombM_oneMatrix(ii,:));
    rectangle('Position',[xMi-0.5 ii-0.5 1 1],'EdgeColor',[0/255, 90/255, 50/255], 'Curvature',0.5); % Create a rectangle

    pIndM_oneMatrix(ii,k) = xMMin; % Store value of the individual minimum in the matrix storing all values from all morphological parameters
end
dRound = 2;
%xticklabels(round([midBin{k}(5), midBin{k}(10), midBin{k}(15), midBin{k}(20), midBin{k}(25), midBin{k}(30)],dRound));
xlabel('Morphological parameters');
%title('Combined minimum for all combinations per parameter');


% Plot the combined matrix with minimum p-value from individual
% combinations

nexttile([1 2]);
%nexttile;
imagesc(pIndM_oneMatrix);
%heatmap(pIndM_oneMatrix);

set(gca,'ColorScale','log');
yticks([1,2,3,4,5,6]); yticklabels({'AA-ABeta', 'AA-AS','AA-SCD', 'ABeta-AS', 'ABeta-SCD', 'AS-SCD'});
cmap1 = colormap(sky);
axis image; % so that the boxes in heatmaps are squares

colormap(flipud(cmap1));
clim([min(min(pIndM_oneMatrix)) 0.05]);
%clim([min(min(pIndM_oneMatrix)) max(max(pIndM_oneMatrix))]);

cbar = colorbar;
cbar.Label.String = "p-value";

% Ticks for colorbar - select 3 ticks: max at 10^-2; min at 10^-nMe, where
% nMe is 10 more than (so that label is not hidden) than exponent nM
% rounded to nearest even number; and one between the two

nM =floor( log10(min(min(pIndM_oneMatrix))));
nMe = 2*floor(nM/2) + 10;
nMid = (-2 + nMe)/2;
if nMe < -300
    nMe = -300;  % If limit is less than -300, the labelling of limits in the colorbar is in weird notation, so instead show 10^-300 location instead of max
end

% Set three ticks in the colorbar, similar to the command:
% cbar.Ticks = [1e-180, 1e-91, 1e-2];
eval(strcat('cbar.Ticks = [1e',num2str(nMe), ',1e',num2str(nMid), ',1e-2];' ));

% Set a threshold for p value and color all above this threshold light
% orange (so all combinations that are not significantly different are
% colored orange).

pTh = 0.05; % 5% threshold
[xM, yM] = size(pIndM_oneMatrix);
for ii = 1:1:xM
    for jj = 1:1:yM
        if pIndM_oneMatrix(ii,jj) >= 0.05
            rectangle('Position',[jj-0.5 ii-0.5 1 1],'FaceColor',[253/255, 208/255, 162/255],'LineStyle','none'); % Create a rectangle
        end
    end
end

% Highlight overall minimum at index pAvgMinInd 116/255, 196/255, 118/255
%rectangle('Position',[pAvgMinInd-0.5 1-0.5 1 6],'EdgeColor',[49/255, 163/255, 84/255]);%, 'Curvature',0.5); % Create a rectangle
rectangle('Position',[pAvgMinInd-0.5 1-0.5 1 6],'EdgeColor',[180/255, 180/255, 180/255]);

% Highlight the minimum for each row (between each combination)
for ii = 1:1:xM
    [xMMin, xMi] = min(pIndM_oneMatrix(ii,:));
    rectangle('Position',[xMi-0.5 ii-0.5 1 1],'EdgeColor',[0/255, 90/255, 50/255], 'Curvature',0.5); % Create a rectangle

    pIndM_oneMatrix(ii,k) = xMMin; % Store value of the individual minimum in the matrix storing all values from all morphological parameters
end
dRound = 2;
%xticklabels(round([midBin{k}(5), midBin{k}(10), midBin{k}(15), midBin{k}(20), midBin{k}(25), midBin{k}(30)],dRound));
xlabel('Morphological parameters');
%title('Individual minimum per combination per parameter');

set(findall(gcf,'-property','FontSize'),'FontSize',7);
fnamePNG = strcat('Figures/0_statDiff_allMorphParam.png');
fnameFIG = strcat('Figures/0_statDiff_allMorphParam.fig');
print(gcf,fnamePNG,'-dpng','-r600');
saveas(gcf,fnameFIG);
%close;

toc

%% 6) Find indices with minimum p-values (also for feature selection)
% Find indices with minimum overall and pairwise p-values for the complete
% list of morphological parameters and for all non-dimensional
% morphological parameters
tic
fprintf('Finding indices with minimum p-values \n');

% List of normalized parameters
kNorm = [18,19,20,22,23,24,25,26,28,31,32,33,34,35,36,37,38,39,40];

% pAvgMinInd_Allk(k) stores the index for minimum average p-value
% pIndMinInd_Allk(k,com) stores the index for the minimum p-value for each
% combination (com) and each morphological parameter (k)

% p7rowMinInd_Allk contains indices of minimum of geometric averaged
% p-values and individual minimum p-values (of each combination)
p7rowMinInd_Allk = [pAvgMinInd_Allk; pIndMinInd_Allk];
pMinUniqueAllk = cell(1,length(morphParamName));
pMinUniqueNDMP = cell(1,length(kNorm));

for k = 1:1:length(morphParamName)
    pMinUniqueAllk{k} = unique(p7rowMinInd_Allk(:,k));
end

nFeatures = 0; % Count total number of unique features using non-dimensional paramters
for ki = 1:1:length(kNorm)
    k = kNorm(ki);
    pMinUniqueNDMP{ki} = unique(p7rowMinInd_Allk(:,k));
    nFeatures = nFeatures+ length(pMinUniqueNDMP{ki});
end
fprintf('Total number of unique features in %d non-dimensional parameters = %d. \n', length(kNorm), nFeatures);

toc


%% 7) Save variables
% Save matrices of minimum p-values, and related variables 
% Save list of top features (with minimum p-values)
tic

fprintf("Saving variables storing minimum p-value and related indices.... \n");
save("pValuesEtAl.mat", "pCombM_Allk", "pMinUniqueNDMP", "pMinUniqueAllk", "p7rowMinInd_Allk", "pAvgMin_Allk", "pAvgMinInd_Allk", "pIndM_oneMatrix", "pIndMinInd_Allk", "pCombM_oneMatrix", "morphParamName", "morphParamAxis", "kNorm",'-v7.3');

toc




% ----------------------------------------------------------------------- %
%     Normalized frequency of mean values for 4 groups (selected
%     parameters)
% ----------------------------------------------------------------------- %

%% Plot normalized frequency of mean (tiled for select few)

k1 = [1, 4,5,18, 20,23];
figure;
tiledlayout('flow', 'Padding','tight', 'TileSpacing','compact');
for kk = 1:1:length(k1)
    nexttile;
    k = k1(kk);
    colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [231/255 41/255 138/255], [166/255 118/255 29/255],[35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
    %figure; 
    hold on;
    for i = 1:1:length(group)
        plot(midBin{k},meanHist{i}{k}, 'Color',colors{i});
        xlabel(morphParamAxis{k});
        ylabel('Normalized frequency');
        box on;
    end
    
    % Percentiles plotted in a different loop to keep ordering of legends
    for i = 1:1:length(group)
        x2 = [midBin{k}, fliplr(midBin{k})];
        inBetween = [perc25{i}{k}, fliplr(perc75{i}{k})];
        fill(x2, inBetween, colors{i},'FaceAlpha',0.1, 'LineStyle','none');
    end
    hold off;

end
lgT = legend('AA', 'ABeta','AS','SCD', 'Location', 'northoutside','Orientation','horizontal', 'box', 'off');
lgT.Layout.Tile = 'north';

set(findall(gcf,'-property','FontSize'),'FontSize',7)
fnamePNG = strcat('Figures/0_freqDist_Tiled.png');
fnameFIG = strcat('Figures/0_freqDist_Tiled.fig');
print(gcf,fnamePNG,'-dpng','-r600');
saveas(gcf,fnameFIG);
close;

%% PLOT 2 - Non-dimensional parameters

k1 = [18,19,20,22,23,24,25,26,28,31,32,33,34,35,36,37,38,39,40];
figure;
tiledlayout('flow', 'Padding','tight', 'TileSpacing','tight');
for kk = 1:1:length(k1)
    nexttile;
    k = k1(kk);
    colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [231/255 41/255 138/255], [166/255 118/255 29/255],[35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
    %figure; 
    hold on;
    for i = 1:1:length(group)
        plot(midBin{k},meanHist{i}{k}, 'Color',colors{i});
        xlabel(morphParamAxis{k});
        ylabel('Frequency');
        box on;
    end
    
    % Percentiles plotted in a different loop to keep ordering of legends
    for i = 1:1:length(group)
        x2 = [midBin{k}, fliplr(midBin{k})];
        inBetween = [perc25{i}{k}, fliplr(perc75{i}{k})];
        fill(x2, inBetween, colors{i},'FaceAlpha',0.1, 'LineStyle','none');
    end
    hold off;

end
lgT = legend('AA', 'ABeta','AS','SCD', 'Location', 'northoutside','Orientation','horizontal', 'box', 'off');
lgT.Layout.Tile = 'north';

set(findall(gcf,'-property','FontSize'),'FontSize',7)
fnamePNG = strcat('Figures/0_freqDist_NDMP_Tiled.png');
fnameFIG = strcat('Figures/0_freqDist_NDMP_Tiled.fig');
print(gcf,fnamePNG,'-dpng','-r600');
saveas(gcf,fnameFIG);
close;


% ----------------------------------------------------------------------- %
%                           EXTRA PLOTS
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%            SINGLE IMPLEMENTATION OF MULTIPLE COMPARISON
% ----------------------------------------------------------------------- %

%{
%% One-way ANOVA and multiple comparison test for 4 groups (single implementation)
% Note: Single implementation for one k value (morphological parameter) and
% one threshold or MPt value

tic

midBin = AA.hD.midBin;
morphParamLegends = AA.info.morphParamName;
morphParamAxis = AA.info.morphParamAxis;

k = 23; % Eccentricity
MPt = 0.5; % Threshold for morphological parameter

% Find value of the morphological parameter at the threshold MPt for all
% files or coverslips
fMPt = cell(1,length(group));
for i = 1:1:length(group)
    fMPt{i} = zeros(1,length(histData{1,i}));

    for j = 1:1:length(histData{1,i})
        [ d, ix ] = min( abs( midBin{k}- MPt ) );
        fMPt{i}(j) = histData{1,i}{1,j}{1,k}(ix);
    end
end

% Combine all 4 groups
Xbox = [fMPt{1} fMPt{2} fMPt{3} fMPt{4}];
grp = [ones(size(fMPt{1})) 2.*ones(size(fMPt{2})) 3.*ones(size(fMPt{3})) 4.*ones(size(fMPt{4}))];

% One-way ANOVA 
% https://www.mathworks.com/help/stats/one-way-anova.html
[p,tbl,stats] = anova1(Xbox, grp);

% Kruskal-Wallis test - non-parametric version
% https://www.mathworks.com/help/stats/kruskalwallis.html
%[p,tbl,stats] = kruskalwallis(Xbox, grp);

close; close; % Close the two plots from the One-way ANOVA

% For plotting box plot
figure; boxplot(Xbox,grp,'Notch', 'on');
xticklabels({'AA', 'ABeta', 'AS', 'SCD'});
ylabel('Normalized frequency at MP_t');
title(morphParamLegends{k});

% Multiple compare
figure;
results = multcompare(stats, "CriticalValueType","scheffe");

pComb = results(:,6); % Stores the p values for all the combinations
pMinAvg = mean(results(:,6)); % Stores the average of p values for all combinations


% Visualize p value
figure; semilogy(pComb,'o'); 
hold on; plot([1 6], [0.05 0.05], '--'); hold off;
xticks([1 2 3 4 5 6]);
xticklabels({'AA-ABeta', 'AA-AS','AA-SCD', 'ABeta-AS', 'ABeta-SCD', 'AS-SCD'});
ylabel('p value');
title(strcat(morphParamLegends{k}, {' '},'at MP_t =',num2str(MPt)));

toc


%% One-way ANOVA and multiple comparison test for 4 groups (multiple implementation - optimization)
% Note: Single implementation for one k value (morphological parameter) and
% all threshold or MPt values to find combination with minimum p values
% (iterate through all midBin{k} values

tic

midBin = AA.hD.midBin;
morphParamLegends = AA.info.morphParamName;
morphParamName = AA.info.morphParamName;
morphParamAxis = AA.info.morphParamAxis;

k = 39; %23; % Eccentricity
%MPt = 0.5; % Threshold for morphological parameter

pAvgL = zeros(1,length(midBin{k})); % Array storing average p values
pCombM = zeros(6,length(midBin{k}));    % Matrix storing all p values
for MPtInd = 1:1:length(midBin{k})

    MPt = midBin{k}(MPtInd);    % Assign value for MPt

    % Find value of the morphological parameter at the threshold MPt for all
    % files or coverslips
    fMPt = cell(1,length(group));
    for i = 1:1:length(group)
        fMPt{i} = zeros(1,length(histData{1,i}));

        for j = 1:1:length(histData{1,i})
            [ d, ix ] = min( abs( midBin{k}- MPt ) );
            fMPt{i}(j) = histData{1,i}{1,j}{1,k}(ix);
        end
    end

    % Combine all 4 groups
    Xbox = [fMPt{1} fMPt{2} fMPt{3} fMPt{4}];
    grp = [ones(size(fMPt{1})) 2.*ones(size(fMPt{2})) 3.*ones(size(fMPt{3})) 4.*ones(size(fMPt{4}))];

    % One-way ANOVA
    % https://www.mathworks.com/help/stats/one-way-anova.html
    [p,tbl,stats] = anova1(Xbox, grp);

    % Kruskal-Wallis test - non-parametric version
    % https://www.mathworks.com/help/stats/kruskalwallis.html
    %[p,tbl,stats] = kruskalwallis(Xbox, grp);

    close; close; % Close the two plots from the One-way ANOVA

    % Multiple compare
    figure;
    results = multcompare(stats, "CriticalValueType","scheffe");
    close; % Close plot for multiple comparison test

    pCombM(:,MPtInd) = results(:,6); % Stores the p values for all the combinations
    %pAvgL(MPtInd) = mean(results(:,6)); % Stores the average of p values for all combinations
    pAvgL(MPtInd) = geomean(results(:,6)); % Stores the geometric mean of p values for all combinations

end

pCombM(isnan(pCombM))=1; % Replace all NaN with 1 (NaN occurs when comparisons cannot be made; 1 indicates no difference)

[pAvgMin, pAvgMinInd] = min(pAvgL); % Find the minimum averaged p value and index
MPt = midBin{k}(pAvgMinInd);    % Find the corresponding MPt with the minimum average p value

fMPt = cell(1,length(group));
for i = 1:1:length(group)
    fMPt{i} = zeros(1,length(histData{1,i}));

    for j = 1:1:length(histData{1,i})
        [ d, ix ] = min( abs( midBin{k}- MPt ) );
        fMPt{i}(j) = histData{1,i}{1,j}{1,k}(ix);
    end
end
toc

%% Plotting

tic
% Combine all 4 groups
Xbox = [fMPt{1} fMPt{2} fMPt{3} fMPt{4}];
grp = [ones(size(fMPt{1})) 2.*ones(size(fMPt{2})) 3.*ones(size(fMPt{3})) 4.*ones(size(fMPt{4}))];

results = multcompare(stats, "CriticalValueType","scheffe");
close; 

pComb = results(:,6); % Stores the p values for all the combinations
pMinAvg = mean(results(:,6)); % Stores the average of p values for all combinations

% --- Figure 1 --------------------

figure
tL = tiledlayout('flow', 'Padding','tight', 'TileSpacing','compact');
nexttile([1 2]);
%nexttile;
imagesc(pCombM); 
%heatmap(pCombM);

set(gca,'ColorScale','log');
yticks([1,2,3,4,5,6]); yticklabels({'AA-ABeta', 'AA-AS','AA-SCD', 'ABeta-AS', 'ABeta-SCD', 'AS-SCD'});
cmap1 = colormap(sky);
axis image; % so that the boxes in heatmaps are squares

pCombM(pCombM == 0) = 1e-300;   % For any p value = 0, make it equate 1*10^-350, which is less than lowest recorded p (~10^-300)

colormap(flipud(cmap1));
clim([min(min(pCombM)) 0.05]);
%clim([min(min(pCombM)) max(max(pCombM))]);

cbar = colorbar; 
cbar.Label.String = "p-value";

% Ticks for colorbar - select 3 ticks: max at 10^-2; min at 10^-nMe, where
% nMe is 10 more than (so that label is not hidden) than exponent nM
% rounded to nearest even number; and one between the two

nM =floor( log10(min(min(pCombM))));
nMe = 2*floor(nM/2) + 10;
if nMe < -300
    nMe = -300;  % If limit is less than -300, the labelling of limits in the colorbar is in weird notation, so instead show 10^-300 location instead of max
end
nMid = (-2 + nMe)/2;

% Set three ticks in the colorbar, similar to the command:
% cbar.Ticks = [1e-180, 1e-91, 1e-2];
eval(strcat('cbar.Ticks = [1e',num2str(nMe), ',1e',num2str(nMid), ',1e-2];' ));

% Set a threshold for p value and color all above this threshold light
% orange (so all combinations that are not significantly different are
% colored orange).

pTh = 0.05; % 5% threshold
[xM, yM] = size(pCombM);
for ii = 1:1:xM
    for jj = 1:1:yM
        if pCombM(ii,jj) >= 0.05
            rectangle('Position',[jj-0.5 ii-0.5 1 1],'FaceColor',[253/255, 208/255, 162/255],'LineStyle','none'); % Create a rectangle 
        end
    end
end

% Highlight overall minimum at index pAvgMinInd 116/255, 196/255, 118/255
%rectangle('Position',[pAvgMinInd-0.5 1-0.5 1 6],'EdgeColor',[49/255, 163/255, 84/255]);%, 'Curvature',0.5); % Create a rectangle
rectangle('Position',[pAvgMinInd-0.5 1-0.5 1 6],'EdgeColor',[180/255, 180/255, 180/255]);

pCombDiffBin = zeros(6,1); % Stores minimum of each row or combination (can be different bin)
pCombDiff_binIndex = zeros(6,1);
% Highlight the minimum for each row (between each combination)
for ii = 1:1:xM
    [xMMin, xMi] = min(pCombM(ii,:));
    rectangle('Position',[xMi-0.5 ii-0.5 1 1],'EdgeColor',[0/255, 90/255, 50/255], 'Curvature',0.5); % Create a rectangle
    pCombDiffBin(ii) = xMMin;
    pCombDiff_binIndex(ii) = xMi;
end
dRound = 2;
xticklabels(round([midBin{k}(5), midBin{k}(10), midBin{k}(15), midBin{k}(20), midBin{k}(25), midBin{k}(30)],dRound));
xlabel(morphParamAxis{k});

% --- Figure 2 --------------------

colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [231/255 41/255 138/255], [166/255 118/255 29/255],[35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
%figure; 
nexttile
hold on;
for i = 1:1:length(group)
    plot(midBin{k},meanHist{i}{k}, 'Color',colors{i});
    xlabel(morphParamAxis{k});
    ylabel('Normalized frequency');
end

% Percentiles plotted in a different loop to keep ordering of legends
for i = 1:1:length(group)
    x2 = [midBin{k}, fliplr(midBin{k})];
    %curve1 = meanHist{i}{k} + stdHist{i}{k}; curve2 = meanHist{i}{k} - stdHist{i}{k}; inBetween = [curve1, fliplr(curve2)];
    inBetween = [perc25{i}{k}, fliplr(perc75{i}{k})];
    fill(x2, inBetween, colors{i},'FaceAlpha',0.1, 'LineStyle','none');
end
%legend('AA', 'ABeta','AS','SCD', 'Location', 'northoutside','Orientation','horizontal', 'box', 'off');
hold off
% Highlight overall minimum at index pAvgMinInd 
dx = (midBin{k}(ix+1) - midBin{k}(ix))/2;
yL = ylim;
rectangle('Position',[midBin{k}(ix-1)+dx yL(1) dx*2 yL(2)-yL(1)],'EdgeColor',[180/255, 180/255, 180/255]);%, 'Curvature',0.5); % Create a rectangle
ylim([yL(1), yL(2)]);
box on;

% --- Figure 3 --------------------


nexttile;
% For plotting box plot
%figure; 
boxplot(Xbox,grp,'Notch', 'on', 'Colors',[27/255 158/255 119/255; 217/255 95/255 2/255; 117/255 112/255 179/255; 231/255 41/255 138/255],'Symbol','.' ); %'PlotStyle','compact'
xticklabels({'AA', 'ABeta', 'AS', 'SCD'});
%ylabel(strcat('Normalized frequency at MP_t =',{' '},num2str(MPt)));
ylabel('Normalized frequency');
%title(morphParamLegends{k});

% One-way ANOVA
% https://www.mathworks.com/help/stats/one-way-anova.html
[p,tbl,stats] = anova1(Xbox, grp);

% Kruskal-Wallis test - non-parametric version
% https://www.mathworks.com/help/stats/kruskalwallis.html
%[p,tbl,stats] = kruskalwallis(Xbox, grp);

close; close; % Close the two plots from the One-way ANOVA

% --- Figure 4 --------------------

% Plot mean and standard deviation
%figure;
nexttile;
hold on;
eb1 = errorbar(1,mean(fMPt{1}),std(fMPt{1}), 'o', 'MarkerSize', 3, 'MarkerEdgeColor',[19/255 129/255 96/255], 'MarkerFaceColor',[27/255 158/255 119/255]);
eb2 = errorbar(2,mean(fMPt{2}),std(fMPt{2}), 'o', 'MarkerSize', 3, 'MarkerEdgeColor',[179/255 78/255 1/255],'MarkerFaceColor',[217/255 95/255 2/255]);
eb3 = errorbar(3,mean(fMPt{3}),std(fMPt{3}), 'o', 'MarkerSize', 3, 'MarkerEdgeColor',[94/255 90/255 142/255],'MarkerFaceColor',[117/255 112/255 179/255]);
eb4 = errorbar(4,mean(fMPt{4}),std(fMPt{4}), 'o', 'MarkerSize', 3, 'MarkerEdgeColor',[186/255 33/255 110/255],'MarkerFaceColor',[231/255 41/255 138/255]);
hold off;
eb1.Color = [19/255 129/255 96/255];
eb2.Color = [179/255 78/255 1/255];
eb3.Color = [94/255 90/255 142/255];
eb4.Color = [186/255 33/255 110/255];
xticks([1 2 3 4]);
xticklabels({'AA', 'ABeta', 'AS','SCD'});
xlim([0.5 4.5]);
ylabel(strcat('Normalized frequency'));


%{
% Multiple compare
results = multcompare(stats, "CriticalValueType","scheffe");
% Edit/label plots
%xlabel(strcat('Mean of frequency at MP_t =',{' '},num2str(MPt)));
xlabel(strcat('Mean normalized frequency'));
title('');
%ylabel('Label Y');
yticklabels({'SCD', 'AS', 'ABeta','AA'});
set ( gca, 'ydir', 'reverse' , 'xdir', 'reverse');
ax2 = gca; set(ax2, 'YAxisLocation', 'right');
camroll(270);
pComb = results(:,6); % Stores the p values for all the combinations
pMinAvg = mean(results(:,6)); % Stores the average of p values for all combinations
%}


% --- Figure 5 --------------------

% Visualize p value
%figure; 
nexttile;
semilogy(pComb,'.', 'MarkerSize', 12); 
hold on; plot([0 7], [0.05 0.05], '--'); 
semilogy([1 2 3 4 5 6], pCombDiffBin,'o', 'Color', [0/255, 90/255, 50/255],'MarkerSize', 5); 
hold off;
xticks([1 2 3 4 5 6]); xlim([0 7])
xticklabels({'AA-ABeta', 'AA-AS','AA-SCD', 'ABeta-AS', 'ABeta-SCD', 'AS-SCD'});
%ylabel(strcat('p-value at MP_t = ',{' '},num2str(MPt)));
ylabel('p-value');
%title(strcat(morphParamLegends{k}, {' '},'at MP_t =',num2str(MPt)));

% Save
%tL.Units = 'centimeters';
%t.OuterPosition = [0.635 0.635 17 34];
%exportgraphics(tL,'test.png','Resolution',300)

if not(isfolder("/Figures"))
    mkdir("Figures");
end

%set(gcf, 'Units', 'centimeters', 'Position', [4, 4, 17, 10.51], 'PaperUnits', 'Inches', 'PaperSize', [17, 17]);
set(findall(gcf,'-property','FontSize'),'FontSize',7)
fnamePNG = strcat('Figures/1_',num2str(k),'_statDiff_',morphParamLegends{k}, '_MPt_',num2str(MPt),'.png');
fnameFIG = strcat('Figures/1_',num2str(k),'_statDiff_',morphParamLegends{k},'_MPt_',num2str(MPt),'.fig');
print(gcf,fnamePNG,'-dpng','-r600');
saveas(gcf,fnameFIG);
% close;
toc


% ----------------------------------------------------------------------- %
%                 INDIVIDUAL HEATMAP
% ----------------------------------------------------------------------- %

%% Plot heatmap of p-values for a single morphological parameters (all bins)
% Plot individually, after analysis (One-way ANOVA and multiple comparison
% test)

figure; 
imagesc(pCombM); 
%heatmap(pCombM);

set(gca,'ColorScale','log');
yticks([1,2,3,4,5,6]); yticklabels({'AA-ABeta', 'AA-AS','AA-SCD', 'ABeta-AS', 'ABeta-SCD', 'AS-SCD'});
cmap1 = colormap(sky);
axis image; % so that the boxes in heatmaps are squares

colormap(flipud(cmap1));
clim([min(min(pCombM)) 0.05]);
%clim([min(min(pCombM)) max(max(pCombM))]);

cbar = colorbar; 
cbar.Label.String = "p-value";

% Ticks for colorbar - select 3 ticks: max at 10^-2; min at 10^-nMe, where
% nMe is 10 more than (so that label is not hidden) than exponent nM
% rounded to nearest even number; and one between the two

nM =floor( log10(min(min(pCombM))));
nMe = 2*floor(nM/2) + 10;
nMid = (-2 + nMe)/2;

% Set three ticks in the colorbar, similar to the command:
% cbar.Ticks = [1e-180, 1e-91, 1e-2];
eval(strcat('cbar.Ticks = [1e',num2str(nMe), ',1e',num2str(nMid), ',1e-2];' ));

% Set a threshold for p value and color all above this threshold light
% orange (so all combinations that are not significantly different are
% colored orange).

pTh = 0.05; % 5% threshold
[xM, yM] = size(pCombM);
for ii = 1:1:xM
    for jj = 1:1:yM
        if pCombM(ii,jj) >= 0.05
            rectangle('Position',[jj-0.5 ii-0.5 1 1],'FaceColor',[253/255, 208/255, 162/255],'LineStyle','none'); % Create a rectangle 
        end
    end
end

% Highlight overall minimum at index pAvgMinInd 116/255, 196/255, 118/255
rectangle('Position',[pAvgMinInd-0.5 1-0.5 1 6],'EdgeColor',[49/255, 163/255, 84/255]);%, 'Curvature',0.5); % Create a rectangle

% Highlight the minimum for each row (between each combination)
for ii = 1:1:xM
    [xMMin, xMi] = min(pCombM(ii,:));
    rectangle('Position',[xMi-0.5 ii-0.5 1 1],'EdgeColor',[0/255, 90/255, 50/255], 'Curvature',0.5); % Create a rectangle
end
dRound = 2;
xticklabels(round([midBin{k}(5), midBin{k}(10), midBin{k}(15), midBin{k}(20), midBin{k}(25), midBin{k}(30)],dRound));
xlabel(morphParamAxis{k});

%}

% ----------------------------------------------------------------------- %
%        Normalized frequency - mean values for different groups
% ----------------------------------------------------------------------- %


%{
%% Plot normalized frequency of mean
k=23;
colors = {[27/255 158/255 119/255], [217/255 95/255 2/255], [117/255 112/255 179/255], [231/255 41/255 138/255], [166/255 118/255 29/255],[35/255 110/255 90/255],[115/255 85/255 140/255], [180/255 40/255 90/255], [180/255 40/255 90/255] }; % First five colors taken from ColorBrewer2 Dark2
figure; 
hold on;
for i = 1:1:length(group)
    plot(midBin{k},meanHist{i}{k}, 'Color',colors{i});
    xlabel(morphParamAxis{k});
    ylabel('Normalized frequency');
end

% Percentiles plotted in a different loop to keep ordering of legends
for i = 1:1:length(group)
    x2 = [midBin{k}, fliplr(midBin{k})];
    inBetween = [perc25{i}{k}, fliplr(perc75{i}{k})];
    fill(x2, inBetween, colors{i},'FaceAlpha',0.1, 'LineStyle','none');
end
legend('AA', 'ABeta','AS','SCD', 'Location', 'northoutside','Orientation','horizontal', 'box', 'off');
% Highlight overall minimum at index pAvgMinInd 
dx = (midBin{k}(ix+1) - midBin{k}(ix))/2;
yL = ylim;
rectangle('Position',[midBin{k}(ix-1)+dx yL(1) dx*2 yL(2)-yL(1)],'EdgeColor',[180/255, 180/255, 180/255]);%, 'Curvature',0.5); % Create a rectangle
hold off;
box on;
%}


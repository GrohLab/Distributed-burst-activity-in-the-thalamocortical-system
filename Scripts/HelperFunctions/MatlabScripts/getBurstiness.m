%% Gather data
close all; clearvars; clc

% Access correct individual
scriptFullPath = matlab.desktop.editor.getActiveFilename();
load(regexprep(scriptFullPath, 'Scripts.*', 'Scripts\\userDataPath.mat'), 'cohortPath');
try
    load(fullfile(cohortPath, 'animalData.mat'))
    load(fullfile(cohortPath, 'allFiles.mat'),'FileInfo')
    updatedFolders = fullfile(cohortPath,{FileInfo.folder});
    [FileInfo.folder] = deal(updatedFolders{:});
catch
end

fileSelection = cellfun(@(x) fullfile(fileparts(x),'intan-signals\automatedCuration'), {FileInfo.folder}, 'UniformOutput', false)';
answer = listdlg('ListString',fileSelection,...
    'PromptString','Choose sessions to plot.',...
    'ListSize',[600 350]);
fileSelection = fileSelection(answer,:);

%% Plot data

area_names = {'BC','VPM','POm','ZIv'};
totalBurstRatio = table('Size',[3,4],'VariableTypes',repmat({'cell'},1,4),'VariableNames',area_names,'RowNames',{'mean','std','no_units'});

for ses = 1:height(fileSelection)
    dataInfo = dir(fullfile(fileSelection{ses},'*all_channels.mat'));
    try
        load(fullfile(dataInfo(1).folder,dataInfo(1).name),'sortedData','fs')
        sortedData = sortedData(all(~cellfun(@isempty, sortedData),2),:);
        clInfo = getClusterInfo(fullfile(fileSelection{ses},'cluster_info.tsv'));
    catch
        try
            [sortedData, fs] = importPhyFiles(fileSelection{ses});
            sortedData = sortedData(all(~cellfun(@isempty, sortedData),2),:);
            clInfo = getClusterInfo(fullfile(fileSelection{ses},'cluster_info.tsv'));
        catch
            fprintf('\nCould not gather sortedData and clInfo.\n')
            return
        end
    end

    str_idx = regexp(fileSelection{ses},'#\d*','end');
    animalPath = fileSelection{ses}(1:str_idx);
    targetEstimate = true;
    % Include only those units, which tetrodes were within the recording area
    try
        load(fullfile(animalPath,'targetHit.mat'),'targetHit')
    catch
        fprintf(2,'\nNo targetHit.mat file. Progress with unfiltered analysis...\n')
        targetEstimate = false;
    end

    clInfo.group(cell2mat(sortedData(:,3))==1) = deal({'good'});
    clInfo.group(cell2mat(sortedData(:,3))==2) = deal({'mua'});
    clInfo.group(cell2mat(sortedData(:,3))==3) = deal({'noise'});

    include = ismember(clInfo.group,{'good','mua'}) & clInfo.isolation_distance > 15 & clInfo.isi_viol < 3 & ismember(clInfo.ch,targetHit);
    sortedData = sortedData(include,:);

    % Label units with their respective area as subscripts
    sortedData(:,4) = deal({'NaN'});
    for i = 1:size(sortedData,1)
        depth = clInfo.depth(strcmp(clInfo.id,char(sortedData{i,1})));
        switch depth
            case 1
                sortedData{i,4} = 'BC';
            case 1400
                sortedData{i,4} = 'POm';
            case 1700
                sortedData{i,4} = 'VPM';
            case 2400
                sortedData{i,4} = 'ZIv';
        end
    end
    
    % Calculate burstiness of single cells and plot CDF of all clusters
    TypeOfSpike = cell(1,height(sortedData));
    BurstStarts = cell(1,height(sortedData));
    BurstEnds = cell(1,height(sortedData));
    BurstDurations = cell(1,height(sortedData));
    NumSpikesInBursts = cell(1,height(sortedData));
    BurstFrequencies = cell(1,height(sortedData));
    numberBursts = nan(1,height(sortedData));
    burstRatio = nan(1,height(sortedData));
    
    CDFvals = 0:0.01:1;
    ttlCDFs = cell(1,length(CDFvals));
    
    figure
    hold on
    for i = 1:height(sortedData)
        spikeTimes_msec = sortedData{i,2}*1000;
        ISIs = diff(spikeTimes_msec);
        if ~isempty(ISIs)
            [f,x] = ecdf(ISIs);
            plot(x,f,'Color','#bfbfbf')
            
            for iCDF = 1:numel(CDFvals)
                pos = find(f<=CDFvals(iCDF),1,'last');
                ttlCDFs{iCDF} = [ttlCDFs{iCDF}, x(pos)];
            end
            
            [TypeOfSpike_indi,BurstStarts_indi,BurstEnds_indi,BurstDurations_indi,NumSpikesInBursts_indi,numberBursts_indi] = BurstDetect(spikeTimes_msec, 10, 15, 2);
            TypeOfSpike{i} = TypeOfSpike_indi;
            BurstStarts{i} = BurstStarts_indi;
            BurstEnds{i} = BurstEnds_indi;
            BurstDurations{i} = BurstDurations_indi;
            % FOR THE BURST FREQUENCY YOU NEED TO CALCULATE THE MEAN ISI
            % AND THEN DIVIDE OVER A SECOND
            if numberBursts_indi==0
                BurstFrequencies{i} = NaN;
            else
                BurstFrequencies{i} = nan(numberBursts_indi,1);
            end
            start_TOTAL = nan(1,numberBursts_indi);
            for burst = 1:numberBursts_indi
                start_idx = find(BurstStarts_indi(burst)==spikeTimes_msec);
                start_TOTAL(burst) = start_idx;
                mean_isi = mean(ISIs(start_idx:start_idx+NumSpikesInBursts_indi(burst)-2));
                BurstFrequencies{i}(burst) = 1000*(1/mean_isi);
            end
            NumSpikesInBursts{i} = NumSpikesInBursts_indi;
            numberBursts(i) = numberBursts_indi;
            burstRatio(i) = numberBursts_indi/(sum(~logical(TypeOfSpike_indi))+numberBursts_indi);
        end
    end
    xlim([0 5000])
    xticks([1,2,5,10,20,50,100,200,500,1000,5000])
    set(gca,'xscale','log')
    title('Cumlative distribution function of ISIs_{All Clusters}')
    ylabel('F(x)')
    xlabel('ISI_{msec}')
    
    avgCDFs = cellfun(@median, ttlCDFs);
    plot(avgCDFs,CDFvals,'k','LineWidth',2)
    hold off
    
    unitIDs = sortedData(:,1);
    % Save data in file folder
    save(fullfile(fileSelection{ses},'BurstinessData.mat'),...
        'TypeOfSpike','BurstStarts','BurstEnds','BurstDurations',...
        'NumSpikesInBursts','BurstFrequencies','numberBursts','burstRatio','unitIDs')
    
    % Plot burst data
    figure
    hold on
    for ar = 1:numel(area_names)
        plot(burstRatio(cellfun(@(x) isequal(x,area_names{ar}), sortedData(:,4))),'.--','MarkerSize',20)
        totalBurstRatio{'mean',ar}{:} = [totalBurstRatio{'mean',ar}{:},mean(burstRatio(cellfun(@(x) isequal(x,area_names{ar}), sortedData(:,4))),'omitnan')];
        totalBurstRatio{'std',ar}{:} = [totalBurstRatio{'std',ar}{:},std(burstRatio(cellfun(@(x) isequal(x,area_names{ar}), sortedData(:,4))),'omitnan')];
        totalBurstRatio{'no_units',ar}{:} = [totalBurstRatio{'no_units',ar}{:},numel(burstRatio(cellfun(@(x) isequal(x,area_names{ar}), sortedData(:,4))))];
    end
    title('Burtiness of clusters')
    xlabel('No. Clusters')
    ylabel('Burst ratio')
    legend(area_names)
    
    % Plot Cumulative distribution function (CDF) of individual areas
    avgCDFs = table('Size',[1,4],'VariableTypes',repmat({'cell'},1,4),'VariableNames',area_names,'RowNames',{'avgCDF'});
    for ar = 1:numel(area_names)
        idx = find(cellfun(@(x) isequal(x,area_names{ar}), sortedData(:,4)));
        
        CDFvals = 0:0.01:1;
        ttlCDFs = cell(1,length(CDFvals));
        
        figure
        hold on
        for i = 1:numel(idx)
            ISIs = diff(sortedData{idx(i),2}*1000);
            if ~isempty(ISIs)
                [f,x] = ecdf(ISIs);
                plot(x,f,'Color','#bfbfbf')
                
                for iCDF = 1:numel(CDFvals)  % Probably can be vectorized for speed, but this shows the idea
                    pos = find(f<=CDFvals(iCDF),1,'last');
                    ttlCDFs{iCDF} = [ttlCDFs{iCDF}, x(pos)];
                end
            end
        end
        xlim([0 5000])
        xticks([1,2,5,10,20,50,100,200,500,1000,5000])
        set(gca,'xscale','log')
        title(sprintf('Cumlative distribution function of ISIs_{%s}',area_names{ar}))
        ylabel('F(x)')
        xlabel('ISI_{msec}')
        
        % Median cumulative distribution
        avgCDFs = cellfun(@median, ttlCDFs);
        plot(avgCDFs,CDFvals,'k','LineWidth',2)
        hold off
    end
    
    close all
end

%% Plot mean burstiness of individual areas

for ar = 1:numel(area_names)
    figure
    errorbar(totalBurstRatio{'mean',area_names{ar}}{:},totalBurstRatio{'std',area_names{ar}}{:},'o')
    title(sprintf('Mean Burtiness of clusters from %s',area_names{ar}))
    xlabel('Session')
    ylabel('Burst ratio')
    ylim([-0.1 0.7])
end

% Plot overall mean with total std
figure
title('Mean Burtiness of Areas')
xticks([1,2,3,4]),xticklabels(area_names)
ylabel('Burst ratio')
xlim([0.5 4.5])
ylim([0 0.5])
hold on
% Function for weighted total standard deviation
std_fct = @(std,weights) sqrt( sum((weights-1).*std.^2, 'omitnan') / (sum(weights, 'omitnan')-numel(weights)) );
x = 1;
for ar = 1:numel(area_names)
    total_mean = sum(totalBurstRatio{'mean',area_names{ar}}{:}.*(totalBurstRatio{'no_units',area_names{ar}}{:}/sum(totalBurstRatio{'no_units',area_names{ar}}{:})), 'omitnan');
    total_std = std_fct(totalBurstRatio{'std',area_names{ar}}{:},totalBurstRatio{'no_units',area_names{ar}}{:});
    std_of_means = std(totalBurstRatio{'mean',area_names{ar}}{:},'omitnan');
    plot(repmat(x,1,numel(totalBurstRatio{'mean',area_names{ar}}{:})), totalBurstRatio{'mean',area_names{ar}}{:},'.','Color','#bfbfbf')
    errorbar(x,total_mean,std_of_means,'o','MarkerFaceColor','k','Color','k')
    x = x+1;
end


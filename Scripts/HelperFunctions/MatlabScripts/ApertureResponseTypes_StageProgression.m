%% Visualize individual responses upon stimulus presentation
% Choose a stage and plot the stage progression
close all; clearvars; clc;
% Make the user pick the two "ResponsePattern" mat files to compare
scriptFullPath = matlab.desktop.editor.getActiveFilename();
load(regexprep(scriptFullPath, 'Scripts.*', 'Scripts\\userDataPath.mat'), 'cohortPath');
startPath = fullfile(cohortPath,'Analysis-Figures\Burstiness-Analysis');

load(fullfile(cohortPath,'animalData.mat'))

burstResponses = cell(1,2);
tonicResponses = cell(1,2);
trialTypeLabels = cell(1,2);
for trialType = 1:2
    [file,fileDir] = uigetfile('*.mat','Choose ResponsePattern mat file for each trial type.',startPath);
    startPath = fileDir;

    % Check if the ResponsePattern file is up to date
    vars = who('-file', fullfile(fileDir,file));
    if ismember('fileVersion', vars)
        % Load fileVersion if it exists
        load(fullfile(fileDir,file),'fileVersion')
        if fileVersion>=2
            load(fullfile(fileDir,file),'burstAPs_Info','tonicAPs_Info');
            burstResponses{trialType} = burstAPs_Info;
            tonicResponses{trialType} = tonicAPs_Info;
            clear('fileVersion')
        else
            fprintf(2,'\nResponsePattern file is not up to date. Re-run burstsUponTrigger.m script...\n')
            return
        end
    else
        fprintf(2,'\nResponsePattern file is not up to date. Re-run burstsUponTrigger.m script...\n')
        return
    end

    prompt = sprintf('Label %i: ',trialType);
    definput = regexp(file, '(?=narrow|wide)([^_]+)', 'match');
    trialTypeLabels{trialType} = cell2mat(inputdlg(prompt,'Trial type labels',[1 50],{strcat(upper(definput{:}(1)),lower(definput{:}(2:end)))}));

    condition = regexp(file, '_([a-zA-Z_]+)(?=narrow|wide)', 'match','once');
    condition = condition(2:end-1);

    stageDescript = cell2mat(regexp(fileparts(fileDir),'\\([^\\]*)$','tokens','once'));
end

% For retrieving the touch-modulated cells
spontWindow = [-600,-400]; % spontaneous window in msec (default: [-600,-400])
respWindow = [0,200]; % response window in msec (default: [0,200])
singleCellReactivityName = ['SingleCellReactivity_', condition, '_withTargetEstimate_', num2str(respWindow/1000,'respond%.3f-%.3f_'), num2str(spontWindow/1000,'baseline%.3f-%.3f.mat')];

% Decide which areas to plot
area_names = {'BC','VPM','POm','ZIv'};
area_colors = {'#377eb8','#4daf4a','#984ea3','#ff7f00'};

f = figure('Name','Waveform type');
set(gcf,'Position',[1000 600 420 300])
ui_field = gobjects(numel(area_names),1);
ui_text = gobjects(numel(area_names)+1,1);

ui_text(1) = uicontrol(f,'Style','text','Units','normalized',...
    'HorizontalAlignment','left','Position',[0.1 0.85 0.8 0.1],...
    'FontSize',10,'String','Choose the desired subpopulation for each area.');

for i = 1:numel(area_names)
    ui_text(i+1) = uicontrol(f,'Style','text','Units','normalized',...
        'HorizontalAlignment','left','Position',[0.1 0.85-0.1*i 0.1 0.05],...
        'FontSize',10,'String',area_names{i});
    
    if isequal(area_names{i},'ZIv')
        initVal = 3;
    else
        initVal = 2;
    end
    ui_field(i) = uicontrol(f,'Style','popupmenu','Units','normalized',...
        'HorizontalAlignment','left','Position',[0.2 0.85-0.1*i 0.2 0.05],...
        'FontSize',10,'String',{'all','RS units','FS units'},'Value', initVal);
end

doneButton = uicontrol(f,'Style','pushbutton','units','normalized',...
    'Position',[0.8 0.1 0.1 0.05],'String','Done',...
    'Callback',{@doneExe,f,ui_field,area_names,...
    burstResponses,tonicResponses});

waitfor(doneButton)

animalNum = cell2mat(cellfun(@(x) str2double(regexp(x,'#(\d+)','tokens','once')), burstResponses{1}.SessionName, 'UniformOutput', false));
stageNum = cell2mat(cellfun(@(x) str2double(regexp(x,'P3\.(\d+\.\d+|\d+)','tokens','once')), burstResponses{1}.SessionName, 'UniformOutput', false));
sessionNum = cell2mat(cellfun(@(x) str2double(regexp(x,'ses\D*(\d+)','tokens','once')), burstResponses{1}.SessionName, 'UniformOutput', false));

assert(isscalar(unique(stageNum)),'So far, only works with one picked stage.')
%% Calculate relevant indices
% Calculate burst-index and burst-preference, e.g.,
% burstIndex [-1;1] = burstProportion_post - burstProportion_pre
% burstPreference [-1;1] = 0.5*(burstIndex_wide - burstIndex_narrow)
burstPreference = 0.5 * (cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)-x(:,1), burstResponses{2}.CellResponses, 'UniformOutput',false)) - ...
    cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)-x(:,1), burstResponses{1}.CellResponses, 'UniformOutput',false)));

% Perform a ranksum test on the results
burstSign = false(numel(burstPreference),1);
burstSign(~isnan(burstPreference)) = cellfun(@(y,z) ranksum(y,z), ...
    cellfun(@(x) x(:,2)-x(:,1), burstResponses{2}.CellResponses(~isnan(burstPreference)), 'UniformOutput',false),...
    cellfun(@(x) x(:,2)-x(:,1), burstResponses{1}.CellResponses(~isnan(burstPreference)), 'UniformOutput',false))<=0.05;

% Perform a ranksum test on tonic responses
tonicNans = isnan(tonicResponses{1}.Significance) | isnan(tonicResponses{2}.Significance);
tonicSign = false(numel(tonicNans),1);
tonicSign(~tonicNans) = cellfun(@(y,z) ranksum(y,z), cellfun(@(x) x(:,2)-x(:,1), tonicResponses{2}.CellResponses(~tonicNans), 'UniformOutput',false),...
    cellfun(@(x) x(:,2)-x(:,1), tonicResponses{1}.CellResponses(~tonicNans), 'UniformOutput',false))<=0.05;

cellSpecs = listdlg('PromptString','Which cells would you like to analyze?', ...
    'Name','Cell types','ListString',{'All cells', 'Only responding with burst bias', 'Only responding with tonic bias',...
    'Only touch-modulated cells'}, ...
    'ListSize',[250,150],'SelectionMode','single');

if cellSpecs==4 % Touch-modulated cells
    fileSelection = unique(burstResponses{1}.SessionName);
    touchModulated = false(numel(burstPreference),1);
    for file = 1:height(fileSelection)
        load(fullfile(fileSelection{file},singleCellReactivityName),'SingleCellResponse')
        % Match the selected units from the tonic and burst responses with
        % the respective ones from the SingleCellResponse file
        for i = find(ismember(burstResponses{1}.SessionName,fileSelection{file}))'
            idx = find(SingleCellResponse.cluster_id==str2double(burstResponses{1}.UnitID{i}));
            touchModulated(i) = SingleCellResponse.respWide(idx) | SingleCellResponse.respNarrow(idx) | SingleCellResponse.respInterm(idx);
        end
    end
end

% For each animal and session calculate the mean burst bias per area
% Output is a cell array with 4 cells (for each area).
% In each of these cells there is a cell for each animal, containing the
% average values for each session.
BurstBiasProgression = cell(1,numel(area_names));
burstCellProgression = cell(1,numel(area_names));
burstCellProgression_type1 = cell(1,numel(area_names));
burstCellProgression_type2 = cell(1,numel(area_names));

% Check how many animals are analyzed 
indiAnimals = unique(animalNum);

dPrimeProgression = cell(1,numel(indiAnimals));
for anim = 1:numel(indiAnimals)
    animalIdx = contains({animalData.cohort(12).animal.animalName},num2str(indiAnimals(anim)));
    dPrimeProgression{anim} = animalData.cohort(12).animal(animalIdx).dvalues_sessions(animalData.cohort(12).animal(animalIdx).stage_sessionCount==unique(stageNum));
end

% Check how many sessions exist for each animal
animalIdx = contains({animalData.cohort(12).animal.animalName},cellstr(num2str(indiAnimals)));
indiSesNum = cellfun(@(x) sum(x==unique(stageNum)), {animalData.cohort(12).animal(animalIdx).stage_sessionCount});

dPrimeConcat = vertcat(dPrimeProgression{:})';

for ar = 1:numel(area_names)
    BurstBiasProgression{ar} = cell(1,numel(indiAnimals));
    burstCellProgression{ar} = cell(1,numel(indiAnimals));
    burstCellProgression_type1{ar} = cell(1,numel(indiAnimals));
    burstCellProgression_type2{ar} = cell(1,numel(indiAnimals));
    % Get the indices of all cells that match the area
    if contains(area_names_waveforms{ar},'RS')
        area_idx = ismember(burstResponses{1}.Area,area_names{ar}) & burstResponses{1}.PutExcitatory;
    elseif contains(area_names_waveforms{ar},'FS')
        area_idx = ismember(burstResponses{1}.Area,area_names{ar}) & ~burstResponses{1}.PutExcitatory;
    else
        area_idx = ismember(burstResponses{1}.Area,area_names{ar});
    end

    for anim = 1:numel(indiAnimals)
        % For each session calculate the mean bias
        BurstBiasProgression{ar}{anim} = NaN(1,indiSesNum(anim));
        burstCellProgression{ar}{anim} = NaN(1,indiSesNum(anim));
        burstCellProgression_type1{ar}{anim} = NaN(1,indiSesNum(anim)); % Response specifically for one trial type
        burstCellProgression_type2{ar}{anim} = NaN(1,indiSesNum(anim)); % Response specifically for one trial type
        for ses = 1:indiSesNum(anim)
            switch cellSpecs
                case 1 % All cells
                    BurstBiasProgression{ar}{anim}(ses) = mean(burstPreference(animalNum==indiAnimals(anim) ...
                    & sessionNum==ses & area_idx),'omitnan');
                case 2 % Only burst biased
                    BurstBiasProgression{ar}{anim}(ses) = mean(burstPreference(animalNum==indiAnimals(anim) ...
                    & sessionNum==ses & area_idx & burstSign),'omitnan');
                case 3 % Only tonic biased
                    BurstBiasProgression{ar}{anim}(ses) = mean(burstPreference(animalNum==indiAnimals(anim) ...
                    & sessionNum==ses & area_idx & tonicSign),'omitnan');
                case 4 % Only touch-modulated
                    BurstBiasProgression{ar}{anim}(ses) = mean(burstPreference(animalNum==indiAnimals(anim) ...
                    & sessionNum==ses & area_idx & touchModulated),'omitnan');
            end
            burstCellProgression{ar}{anim}(ses) = sum(burstSign(animalNum==indiAnimals(anim) ...
                & sessionNum==ses & area_idx))/numel(burstSign(animalNum==indiAnimals(anim) ...
                & sessionNum==ses & area_idx));
            % Negative = type 1 is prefered
            burstCellProgression_type1{ar}{anim}(ses) = sum(burstSign(animalNum==indiAnimals(anim) ...
                & sessionNum==ses & area_idx) & burstPreference(animalNum==indiAnimals(anim) ...
                & sessionNum==ses & area_idx)<0)/numel(burstSign(animalNum==indiAnimals(anim) ...
                & sessionNum==ses & area_idx));
            % Positive = type 2 is prefered
            burstCellProgression_type2{ar}{anim}(ses) = sum(burstSign(animalNum==indiAnimals(anim) ...
                & sessionNum==ses & area_idx) & burstPreference(animalNum==indiAnimals(anim) ...
                & sessionNum==ses & area_idx)>0)/numel(burstSign(animalNum==indiAnimals(anim) ...
                & sessionNum==ses & area_idx));
        end
    end
end

% Calculate values for all areas combined
BurstBiasProgression_allAreas = cell(1,numel(indiAnimals));
burstCellProgression_allAreas = cell(1,numel(indiAnimals));
burstCellProgression_type1_allAreas = cell(1,numel(indiAnimals));
burstCellProgression_type2_allAreas = cell(1,numel(indiAnimals));
for anim = 1:numel(indiAnimals)
    % For each session calculate the mean bias
    for ses = 1:indiSesNum(anim)
        switch cellSpecs
            case 1 % All cells
                BurstBiasProgression_allAreas{anim}(ses) = mean(burstPreference(animalNum==indiAnimals(anim) ...
                    & sessionNum==ses),'omitnan');
            case 2 % Only burst biased
                BurstBiasProgression_allAreas{anim}(ses) = mean(burstPreference(animalNum==indiAnimals(anim) ...
                    & sessionNum==ses & burstSign),'omitnan');
            case 3 % Only tonic biased
                BurstBiasProgression_allAreas{anim}(ses) = mean(burstPreference(animalNum==indiAnimals(anim) ...
                    & sessionNum==ses & tonicSign),'omitnan');
            case 4 % Only touch-modulated
                BurstBiasProgression_allAreas{anim}(ses) = mean(burstPreference(animalNum==indiAnimals(anim) ...
                    & sessionNum==ses & touchModulated),'omitnan');
        end
        burstCellProgression_allAreas{anim}(ses) = sum(burstSign(animalNum==indiAnimals(anim) ...
            & sessionNum==ses))/numel(burstSign(animalNum==indiAnimals(anim) ...
                & sessionNum==ses));
        % Negative = type 1 is prefered
        burstCellProgression_type1_allAreas{anim}(ses) = sum(burstSign(animalNum==indiAnimals(anim) ...
            & sessionNum==ses) & burstPreference(animalNum==indiAnimals(anim) ...
            & sessionNum==ses)<0)/numel(burstSign(animalNum==indiAnimals(anim) ...
                & sessionNum==ses));
        % Positive = type 2 is prefered
        burstCellProgression_type2_allAreas{anim}(ses) = sum(burstSign(animalNum==indiAnimals(anim) ...
            & sessionNum==ses) & burstPreference(animalNum==indiAnimals(anim) ...
            & sessionNum==ses)>0)/numel(burstSign(animalNum==indiAnimals(anim) ...
                & sessionNum==ses));
    end
end

%% Plot the burst bias of all animals, normalizing the session count

for ar = 1:numel(area_names)+1
    if ar==numel(area_names)+1 % All areas
        plotColor = 'k';
        plotTitle = 'allAreas';
    else
        plotColor = area_colors{ar};
        plotTitle = area_names_waveforms{ar};
    end
    % Create sub with normalized x-val (stage progression) in the first
    % row and the corresponding burst bias in the second row.
    if ar==numel(area_names)+1 % All areas
        subResponses = NaN(2,numel([BurstBiasProgression_allAreas{:}]));
        subResponses(2,:) = [BurstBiasProgression_allAreas{:}];
    else
        subResponses = NaN(2,numel([BurstBiasProgression{ar}{:}]));
        subResponses(2,:) = [BurstBiasProgression{ar}{:}];
    end
    switch cellSpecs
        case 1 % All cells
            figure('Name',sprintf('BurstBiasProgression_%s_allCells',plotTitle))
        case 2 % Only burst biased
            figure('Name',sprintf('BurstBiasProgression_%s_onlyBurst',plotTitle))
        case 3 % Only tonic biased
            figure('Name',sprintf('BurstBiasProgression_%s_onlyTonic',plotTitle))
        case 4 % Only touch-modulated
            figure('Name',sprintf('BurstBiasProgression_%s_onlyTouchModul',plotTitle))
    end
    hold on
    
    count = 1;
    for anim = 1:numel(indiAnimals)
        if ar==numel(area_names)+1 % All areas
        subResponses(1,count:count+numel(BurstBiasProgression_allAreas{anim})-1) = linspace(0,1,numel(BurstBiasProgression_allAreas{anim}));
%         plot(linspace(0,1,numel(BurstBiasProgression_allAreas{anim})),BurstBiasProgression_allAreas{anim},'o-','Color','#d9d9d9')
        count = count + numel(BurstBiasProgression_allAreas{anim});
        else
        subResponses(1,count:count+numel(BurstBiasProgression{ar}{anim})-1) = linspace(0,1,numel(BurstBiasProgression{ar}{anim}));
%         plot(linspace(0,1,numel(BurstBiasProgression{ar}{anim})),BurstBiasProgression{ar}{anim},'o-','Color','#d9d9d9')
        count = count + numel(BurstBiasProgression{ar}{anim});
        end
    end
    
    % Sort array by normalized x-vals
    [~,sortIdx] = sort(subResponses(1,:));
    subResponses = subResponses(:,sortIdx);
    % If all animals have the same session count, take those as steps
    if  numel(unique(subResponses(1,:)))==unique(indiSesNum)
        movMeanResponses = NaN(1,unique(indiSesNum));
        movMeanDPrime = NaN(1,unique(indiSesNum));
        movSEMResponses = NaN(1,unique(indiSesNum));
        movSEMDPrime = NaN(1,unique(indiSesNum));
        for i = 0:unique(indiSesNum)-1
            fact = length(subResponses)/unique(indiSesNum);
            idx = 1+(i*fact):(1+i)*fact;
            movMeanResponses(i+1) = mean(subResponses(2,idx),'omitnan');
            movMeanDPrime(i+1) = mean(dPrimeConcat(sortIdx(idx)),'omitnan');
            movSEMResponses(i+1) = std(subResponses(2,idx),'omitnan')/sqrt(sum(idx));
            movSEMDPrime(i+1) = std(dPrimeConcat(sortIdx(idx)),'omitnan')/sqrt(sum(idx));
        end
        xvals = 1:unique(indiSesNum);
        plot([1,unique(indiSesNum)],[0,0],'k--')
        xlabel('Sessions')
    else
        windSize = 0.1; % Size to average over
        windStep = 0.025; % Step for next window (smaller than windSize will result in a moving mean)
        assert(rem((1-windSize),windStep)==0,'windSize and windStep must be chosen in a way, that they smoothly fit into the window of 0 to 1.')

        movMeanResponses = NaN(1,(1-windSize)/windStep+1);
        movMeanDPrime = NaN(1,(1-windSize)/windStep+1);
        movSEMResponses = NaN(1,(1-windSize)/windStep+1);
        movSEMDPrime = NaN(1,(1-windSize)/windStep+1);
        for i = 0:(1-windSize)/windStep
            idx = subResponses(1,:)>=i*windStep & subResponses(1,:)<=windSize+i*windStep;
            movMeanResponses(i+1) = mean(subResponses(2,idx),'omitnan');
            movMeanDPrime(i+1) = mean(dPrimeConcat(sortIdx(idx)),'omitnan');
            movSEMResponses(i+1) = std(subResponses(2,idx),'omitnan')/sqrt(sum(idx));
            movSEMDPrime(i+1) = std(dPrimeConcat(sortIdx(idx)),'omitnan')/sqrt(sum(idx));
        end
        
        xvals = linspace(0,1,(1-windSize)/windStep+1);
        plot([0,1],[0,0],'k--')
        xlabel('Learning progression')
    end

    curve1 = movMeanResponses + movSEMResponses;
    curve2 = movMeanResponses - movSEMResponses;
    plot(xvals, curve1,'-','Color',plotColor);
    plot(xvals, curve2,'-','Color',plotColor);
    x = xvals(~isnan(curve1));
    curve1 = curve1(~isnan(curve1));
    curve2 = curve2(~isnan(curve2));
    fill([x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
        'FaceColor',plotColor,'EdgeColor','none','FaceAlpha',0.2);
    plot(xvals,movMeanResponses,'Color',plotColor,'LineWidth',1.5);

    ylabel(sprintf('%s \\leftarrow Burst bias \\rightarrow %s',trialTypeLabels{:}))

    yyaxis right
    curve1 = movMeanDPrime + movSEMDPrime;
    curve2 = movMeanDPrime - movSEMDPrime;
    plot(xvals, curve1,'-','Color','#cccccc');
    plot(xvals, curve2,'-','Color','#cccccc');
    x = xvals(~isnan(curve1));
    curve1 = curve1(~isnan(curve1));
    curve2 = curve2(~isnan(curve2));
    fill([x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
        'FaceColor','#33adff','EdgeColor','none','FaceAlpha',0.2);
    plot(xvals,movMeanDPrime,'Color','#33adff','LineWidth',1.5);
    ylabel('d'' prime')
    set(gca, 'YColor', '#33adff')


    title(sprintf('Change of burst bias with learning - %s',plotTitle))
    hold off
end

%% Plot the burst responding cells of all animals, normalizing the session count

for ar = 1:numel(area_names)+1
    if ar==numel(area_names)+1 % All areas
        plotColor = 'k';
        plotTitle = 'allAreas';
    else
        plotColor = area_colors{ar};
        plotTitle = area_names_waveforms{ar};
    end
    % Create sub with normalized x-val (stage progression) in the first
    % row and the corresponding burst bias in the second row.
    if ar==numel(area_names)+1 % All areas
        subResponses = NaN(2,numel([burstCellProgression_allAreas{:}]));
        subResponses(2,:) = [burstCellProgression_allAreas{:}];
    else
        subResponses = NaN(2,numel([burstCellProgression{ar}{:}]));
        subResponses(2,:) = [burstCellProgression{ar}{:}];
    end

    figure('Name',sprintf('BurstResponderFractionProgression_%s',plotTitle))
    hold on
    count = 1;
    for anim = 1:numel(indiAnimals)
        if ar==numel(area_names)+1 % All areas
            subResponses(1,count:count+numel(burstCellProgression_allAreas{anim})-1) = linspace(0,1,numel(burstCellProgression_allAreas{anim}));
    %         plot(linspace(0,1,numel(burstCellProgression_allAreas{anim})),burstCellProgression_allAreas{anim},'o-','Color','#d9d9d9')
            count = count + numel(burstCellProgression_allAreas{anim});
        else
            subResponses(1,count:count+numel(burstCellProgression{ar}{anim})-1) = linspace(0,1,numel(burstCellProgression{ar}{anim}));
    %         plot(linspace(0,1,numel(burstCellProgression{ar}{anim})),burstCellProgression{ar}{anim},'o-','Color','#d9d9d9')
            count = count + numel(burstCellProgression{ar}{anim});
        end
    end
    
    % Sort array by normalized x-vals
    [~,sortIdx] = sort(subResponses(1,:));
    subResponses = subResponses(:,sortIdx);
    
    if ar==numel(area_names)+1 % All areas
        type1Sort = [burstCellProgression_type1_allAreas{:}];
        type2Sort = [burstCellProgression_type2_allAreas{:}];
    else
        type1Sort = [burstCellProgression_type1{ar}{:}];
        type2Sort = [burstCellProgression_type2{ar}{:}];
    end
    type1Sort = type1Sort(sortIdx);
    type2Sort = type2Sort(sortIdx);

    % If all animals have the same session count, take those as steps
    if  numel(unique(subResponses(1,:)))==unique(indiSesNum)
        movMeanResponses = NaN(1,unique(indiSesNum));
        movMeanResponses_type1 = NaN(1,unique(indiSesNum));
        movMeanResponses_type2 = NaN(1,unique(indiSesNum));
        movMeanDPrime = NaN(1,unique(indiSesNum));
        movSEMResponses = NaN(1,unique(indiSesNum));
        movSEMDPrime = NaN(1,unique(indiSesNum));
        for i = 0:unique(indiSesNum)-1
            fact = length(subResponses)/unique(indiSesNum);
            idx = 1+(i*fact):(1+i)*fact;
            movMeanResponses(i+1) = mean(subResponses(2,idx),'omitnan');
            movMeanResponses_type1(i+1) = mean(type1Sort(idx),'omitnan');
            movMeanResponses_type2(i+1) = mean(type2Sort(idx),'omitnan');
            movMeanDPrime(i+1) = mean(dPrimeConcat(sortIdx(idx)),'omitnan');
            movSEMResponses(i+1) = std(subResponses(2,idx),'omitnan')/sqrt(sum(idx));
            movSEMDPrime(i+1) = std(dPrimeConcat(sortIdx(idx)),'omitnan')/sqrt(sum(idx));
        end

        xvals = 1:unique(indiSesNum);
        xlabel('Sessions')
    else
        windSize = 0.1; % Size to average over
        windStep = 0.025; % Step for next window (smaller than windSize will result in a moving mean)
        assert(rem((1-windSize),windStep)==0,'windSize and windStep must be chosen in a way, that they smoothly fit into the window of 0 to 1.')

        movMeanResponses = NaN(1,(1-windSize)/windStep+1);
        movMeanResponses_type1 = NaN(1,(1-windSize)/windStep+1);
        movMeanResponses_type2 = NaN(1,(1-windSize)/windStep+1);
        movMeanDPrime = NaN(1,(1-windSize)/windStep+1);
        movSEMResponses = NaN(1,(1-windSize)/windStep+1);
        movSEMDPrime = NaN(1,(1-windSize)/windStep+1);
        for i = 0:(1-windSize)/windStep
            idx = subResponses(1,:)>=i*windStep & subResponses(1,:)<=windSize+i*windStep;
            movMeanResponses(i+1) = mean(subResponses(2,idx),'omitnan');
            movMeanResponses_type1(i+1) = mean(type1Sort(idx),'omitnan');
            movMeanResponses_type2(i+1) = mean(type2Sort(idx),'omitnan');
            movMeanDPrime(i+1) = mean(dPrimeConcat(sortIdx(idx)),'omitnan');
            movSEMResponses(i+1) = std(subResponses(2,idx),'omitnan')/sqrt(sum(idx));
            movSEMDPrime(i+1) = std(dPrimeConcat(sortIdx(idx)),'omitnan')/sqrt(sum(idx));
        end

        xvals = linspace(0,1,(1-windSize)/windStep+1);
        xlabel('Learning progression')
    end

    curve1 = movMeanResponses + movSEMResponses;
    curve2 = movMeanResponses - movSEMResponses;
    plot(xvals, curve1,'-','Color',plotColor);
    plot(xvals, curve2,'-','Color',plotColor);
    x = xvals(~isnan(curve1));
    curve1 = curve1(~isnan(curve1));
    curve2 = curve2(~isnan(curve2));
    fill([x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
        'FaceColor',plotColor,'EdgeColor','none','FaceAlpha',0.2);
    plot(xvals,movMeanResponses,'Color',plotColor,'LineWidth',1.5);
    p1 = plot(xvals,movMeanResponses_type1,'Color','#448D76','LineWidth',1.5);
    p2 = plot(xvals,movMeanResponses_type2,'Color','#CBE432','LineWidth',1.5);
    legend([p1,p2],trialTypeLabels, 'Location','northwest')
    legend('AutoUpdate', 'off')

    ylabel('Fraction')

    yyaxis right
    curve1 = movMeanDPrime + movSEMDPrime;
    curve2 = movMeanDPrime - movSEMDPrime;
    plot(xvals, curve1,'-','Color','#cccccc');
    plot(xvals, curve2,'-','Color','#cccccc');
    x = xvals(~isnan(curve1));
    curve1 = curve1(~isnan(curve1));
    curve2 = curve2(~isnan(curve2));
    fill([x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
        'FaceColor','#33adff','EdgeColor','none','FaceAlpha',0.2);
    plot(xvals,movMeanDPrime,'Color','#33adff','LineWidth',1.5);
    ylabel('d'' prime')
    set(gca(), 'YColor','#33adff')

    title(sprintf('Change of burst responders with learning - %s',plotTitle))
    hold off
end

%% Save figures
% Destination folder for matlab .fig files
destfold = fullfile(cohortPath,'Analysis-Figures','Burstiness-StageProgression');
if exist(destfold,"dir") == 0
    mkdir(destfold)
end

figHandles = handle(sort(double(findall(0, 'type', 'figure'))));
            
fprintf("Saving Figures.\n")

figure_suffix = [condition, '_', strjoin(trialTypeLabels,'&')];
for i = 1:numel(figHandles)
    if ~isempty(figHandles(i).Name)
        destfile = sprintf('%s\\%s_stage%i_%s.fig', destfold, figure_suffix, unique(stageNum), figHandles(i).Name);
        savefig(figHandles(i), destfile);
    end
end
fprintf("\nDone!\n\n")


%% Helper functions
% Function executed by the "doneButton" in subpopulation refinement
function doneExe(~,~,f,ui_field,area_names,...
    burstResponses,tonicResponses)

idx = false(height(burstResponses{1}),1);
waveformTypes = cell(1,numel(area_names));
area_names_waveforms = cell(1,numel(area_names));
for i = 1:numel(area_names)
    switch ui_field(i).Value
        case 1 % 'all'
            waveformTypes{i} = 'all units';
            area_names_waveforms{i} = area_names{i};
            idx = idx | contains(burstResponses{1}.Area,area_names{i});
        case 2 % 'RS units'
            waveformTypes{i} = 'RS units';
            area_names_waveforms{i} = [area_names{i},'-RS'];
            idx = idx | (contains(burstResponses{1}.Area,area_names{i}) & ...
                burstResponses{1}.PutExcitatory);
        case 3 % 'FS units'
            waveformTypes{i} = 'FS units';
            area_names_waveforms{i} = [area_names{i},'-FS'];
            idx = idx | (contains(burstResponses{1}.Area,area_names{i}) & ...
                ~burstResponses{1}.PutExcitatory);
    end
end

for trialType = 1:2
    burstResponses{trialType} = burstResponses{trialType}(idx,:);
    tonicResponses{trialType} = tonicResponses{trialType}(idx,:);
end

% Overwrite variables
assignin('caller','burstResponses',burstResponses)
assignin('caller','tonicResponses',tonicResponses)
assignin('caller','waveformTypes',waveformTypes)
assignin('caller','area_names_waveforms',area_names_waveforms)

close(f)
end
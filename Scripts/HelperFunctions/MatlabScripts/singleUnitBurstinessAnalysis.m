%% Examplary burstiness analyses for single units
% This script takes data from the "ClusterBurstiness*.mat"
% files and displays a raster plot, as well as a comparing box plot.
% It also calculates whether the unit has a significantly different amount
% of bursts in the response window, and saves that information in the
% "ClusterBurstiness*.mat" file.

% Variables to set: fileSelection, condition, two trial types, spontWindow,
% respWindow, timeLapse

close all; clearvars; clc

% Choose sessions to merge together
scriptFullPath = matlab.desktop.editor.getActiveFilename();
load(regexprep(scriptFullPath, 'Scripts.*', 'Scripts\userDataPath.mat'), 'cohortPath');
try
    load(fullfile(cohortPath, 'animalData.mat'))
    load(fullfile(cohortPath, 'allFiles.mat'),'FileInfo')
catch
end

prompt = {'Choose individual sessions',...
    'Stage 1','Stage 2','Stage 3','Stage 4','Stage 5','Stage 6',...
    'Stage 7','Stage 8','Stage 9'};
filePick = listdlg('PromptString',{'Which sessions do you want to analyze?';'Possible to refine stage sessions afterwards.'}, ...
    'ListString',prompt,...
    'ListSize',[300 200],'InitialValue',1,'SelectionMode','single');

if filePick==1
    fileSelection = cellfun(@(x) fullfile(fileparts(x),'intan-signals\automatedCuration'), {FileInfo.folder}, 'UniformOutput', false)';
    answer = listdlg('ListString',fileSelection,...
        'PromptString','Choose sessions to plot.',...
        'ListSize',[600 350]);
    fileSelection = fileSelection(answer,:);
else
    stageNum = str2double(regexp(prompt{filePick},'\d*','match','once'));

    all_sesNames = cell(1,6);
    all_dprimes = cell(1,6);
    for i = 1:6
        idx = find(animalData.cohort(12).animal(i).stage_sessionCount==stageNum);
        ses = animalData.cohort(12).animal(i).session_names(idx);
        dvals = animalData.cohort(12).animal(i).dvalues_sessions(idx);
        all_dprimes{i} = dvals;
        all_sesNames{i} = ses;
    end
    % If a certain stage is not represented for an animal, omit this individual
    all_sesNames = all_sesNames(~cellfun(@isempty,all_sesNames));

    prompt = {'d'' -Inf to -1.65', 'd'' -1.65 to -0.5',...
        'd'' -0.5 to 0.5','d'' 0.5 to 1.65',...
        'd'' 1.65 to Inf','First two sessions',...
        'Last two sessions','All sessions of that stage'};
    performancePick = listdlg('PromptString',{'Which sessions from that stage do you want to analyze?';'Possible to refine stage sessions afterwards.'}, ...
        'ListString',prompt,...
        'ListSize',[300 200],'InitialValue',5,'SelectionMode','single');

    if ismember(performancePick,(1:5))
        all_dprimes = vertcat(all_dprimes{:});
        all_sesNames = vertcat(all_sesNames{:});
        [~,~,group] = histcounts(all_dprimes,'BinEdges',[-Inf,-1.65,-0.5,0.5,1.65,Inf]);
        fileSelection = fullfile(fileparts(fileparts(all_sesNames(group==performancePick))),'intan-signals\automatedCuration');
    elseif performancePick==6
        all_sesNames = cellfun(@(x) x(1:2),all_sesNames,'UniformOutput',false);
        all_sesNames = vertcat(all_sesNames{:});
        fileSelection = fullfile(fileparts(fileparts(all_sesNames)),'intan-signals\automatedCuration');
    elseif performancePick==7
        all_sesNames = cellfun(@(x) x(end-1:end),all_sesNames,'UniformOutput',false);
        all_sesNames = vertcat(all_sesNames{:});
        fileSelection = fullfile(fileparts(fileparts(all_sesNames)),'intan-signals\automatedCuration');
    elseif performancePick==8
        all_sesNames = vertcat(all_sesNames{:});
        fileSelection = fullfile(fileparts(fileparts(all_sesNames)),'intan-signals\automatedCuration');
    end
end

if exist('stageNum','var')
    if contains(prompt{performancePick},'two')
        suffix = lower(strrep(prompt{performancePick},' two ','2'));
    elseif contains(prompt{performancePick},'All')
        suffix = 'allSessions';
    else
        suffix = regexp(prompt{performancePick},'[-+]?(\d+\.?\d*|\.\d+)|[-+]?Inf','match');
        suffix = ['dprime',cell2mat(join(suffix,'_'))];
    end
    definput = sprintf('stage%i-%s',stageNum,suffix);
    sessionDescription = inputdlg('Enter a session description for saving the files:','Output Folder',[1 100],{definput});
else
    sessionDescription = inputdlg('Enter a session description for saving the files:','Output Folder',[1 100]);
end

area_names = {'BC','VPM','POm','ZIv'};
area_colors = {'#377eb8','#4daf4a','#984ea3','#ff7f00'};
% Choose time windows to compare
spontWindow = [-600,-400]; % spontaneous window in msec (default: [-600,-400])
respWindow = [0,200]; % response window in msec (default: [0,200])
timeLapse = [-800,800]; % total time window in msec

% Choose condition to analyze
prompt = {'Reward', 'Punishment', 'Lick', 'onlyFirstLick',...
    'WhiskerContact_left', 'WhiskerContact_right','Middlepoint',...
    'Random', 'WhiskerContact_onlyLeftFirst', 'WhiskerContact_onlyRightFirst'};
condPick = listdlg('PromptString','Pick condition to analyze.', ...
    'ListString',prompt,...
    'ListSize',[300 200],'InitialValue',5,'SelectionMode','single');
condition = prompt{condPick};

% Supress figures, when multiple sessions are analyzed
plotLickTimes = false;
plotUnits = false;
if numel(fileSelection)>1
    set(0,'DefaultFigureVisible','off')
else
    answer = questdlg('Would you like to plot individual units?', ...
        'Plot units', ...
        'Yes','No','No');
    if isequal(answer,'No')
        set(0,'DefaultFigureVisible','off');
    elseif contains(condition,'WhiskerContact')
        plotUnits = true;
        % Decide whether you also want to plot lick times in raster plots
        answer = questdlg('Would you also like to plot lick times in the raster plots?', ...
            'Plot lick times', ...
            'Yes','No','Yes');
        if isequal(answer,'Yes')
            plotLickTimes = true;
        end
    else
        plotUnits = true;
    end
end

% Choose two trial types to analyze (more than that will quickly become
% confusing)
trialType = {'all', 'wide', 'wide&lick', 'wide&no-lick',...
    'narrow', 'narrow&lick', 'narrow&no-lick',...
    'intermediate', 'intermediate&lick', 'intermediate&no-lick'};
[trialPick, tf] = listdlg('ListString',trialType,...
    'PromptString','Which trial type(s) do you want to analyze?',...
    'SelectionMode','multiple','ListSize', [250 250],'InitialValue',[2 5]);
if tf==0
    return
end
assert(numel(trialPick)==2,...
    'singleUniteBurstiness:selectInputNumber',...
    'You have to pick 2 trial types.')
trialType = trialType(trialPick);

%% Choose which cells you would like to analyze
cellSpecs = listdlg('PromptString','Which cells would you like to analyze?', ...
    'Name','Cell types','ListString',{'All cells', 'Only responding with burst bias', 'Only responding with tonic bias',...
    'Only responding with bias in both', 'Only responding with bias in none', 'Only touch-modulated cells'}, ...
    'ListSize',[250,150],'SelectionMode','single');

if ismember(cellSpecs,(2:5))
    % Look if you can access the correct "ResponsePattern*.mat" file
    burstResponses = cell(1,numel(trialType));
    tonicResponses = cell(1,numel(trialType));
    for i = 1:numel(trialType)
        startPath = fullfile(cohortPath, 'Analysis-Figures\Burstiness-Analysis');
        RespPatFile = fullfile(startPath,sessionDescription{:}, ...
            sprintf('ResponsePattern_%s_%s_withTargetEstimate_respond%.3f-%.3f_baseline%.3f-%.3f.mat',condition,trialType{i},respWindow,spontWindow));
        if ~exist(RespPatFile,'file')
            answer = questdlg({sprintf('The ResponsePattern.mat file for the trial type "%s" was not found in the folder "./%s".',trialType{i},sessionDescription{:}); ...
                'Would you like to manually look for it?'}, ...
                'ResponsePattern file', ...
                'Yes','No, abort','Yes');
            if isequal(answer,'Yes')
                [file,fileDir] = uigetfile('*.mat','Choose ResponsePattern mat file for each trial type.',startPath);
                RespPatFile = fullfile(fileDir,file);
            else
                return
            end
        end
        
        % Check if the ResponsePattern file is up to date
        vars = who('-file', RespPatFile);
        if ismember('fileVersion', vars)
            % Load fileVersion if it exists
            load(RespPatFile,'fileVersion')
            if fileVersion>=2
                load(RespPatFile,'burstAPs_Info','tonicAPs_Info');
                burstResponses{i} = burstAPs_Info;
                tonicResponses{i} = tonicAPs_Info;
                clear('fileVersion')
            else
                fprintf(2,'\nResponsePattern file is not up to date. Re-run burstsUponTrigger.m script...\n')
                return
            end
        else
            fprintf(2,'\nResponsePattern file is not up to date. Re-run burstsUponTrigger.m script...\n')
            return            
        end
    end

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
    burstSignTable = [burstResponses{1}(:,1:4),table(burstSign)];

    % Perform a ranksum test on tonic responses
    tonicSign = cellfun(@(y,z) ranksum(y,z), cellfun(@(x) x(:,2)-x(:,1), tonicResponses{2}.CellResponses, 'UniformOutput',false),...
        cellfun(@(x) x(:,2)-x(:,1), tonicResponses{1}.CellResponses, 'UniformOutput',false))<=0.05;
    tonicSignTable = [tonicResponses{1}(:,1:4),table(tonicSign)];
end

switch cellSpecs
    case 1
        cellDescript = 'allCells';
    case 2
        cellDescript = 'onlyBurst';
    case 3
        cellDescript = 'onlyTonic';
    case 4
        cellDescript = 'onlyBoth';
    case 5
        cellDescript = 'onlyNone';
    case 6
        cellDescript = 'onlyTouchModul';
end

%% Loop through each file

% Bin size for continuous firing rates
binSize = 10; % in msec
binEdgesFiring = milliseconds(timeLapse(1):binSize:timeLapse(2));

% Bin size for burst event rate
binSizeBursts = 25; % in msec
binEdgesBursts = milliseconds(timeLapse(1):binSizeBursts:timeLapse(2));

% Check if the selected files have already been analyzed
analyzed = false(height(fileSelection),2);
for file = 1:height(fileSelection)
    currFile = fileSelection{file};
    for i = 1:numel(trialType)
        singleCellFiringRateName = ['SingleCellFiringRates_', condition,num2str(timeLapse/1000,'_timeLapse%.3f-%.3f_'),...
                num2str(respWindow/1000,'respWind%.3f-%.3f_'), num2str(spontWindow/1000,'spontWind%.3f-%.3f_'),...
                num2str(binSize,'%imsecbinSize_'),num2str(binSizeBursts,'%imsecbinSizeBursts_'),trialType{i},'_',cellDescript,'.mat'];
        if exist(fullfile(currFile,singleCellFiringRateName),'file')
            vars = who('-file', fullfile(currFile,singleCellFiringRateName));
            if ismember('fileVersion', vars)
                % Load fileVersion if it exists
                load(fullfile(currFile,singleCellFiringRateName),'fileVersion')
                if fileVersion>=2
                    analyzed(file,i) = true;
                    clear('fileVersion')
                end
            end
        end
    end
end
analyzed = all(analyzed,2);

for file = 1:height(fileSelection)
    currFile = fileSelection{file};
    if analyzed(file) && ~plotUnits
        % If file has been analyzed skip the loop
        fprintf('\nLoading file: "%s" (%i/%i)\n',currFile,file,height(fileSelection))
        continue
    else
        fprintf('\nAnalyzing file: "%s" (%i/%i)\n',currFile,file,height(fileSelection))

        % Get the lick events from the file
        load(fullfile(fileparts(fileparts(currFile)),'videos\HispeedTrials.mat'),'HispeedTrials');
        lickEvents = sort(cell2mat(cellfun(@(x,y) x(y), HispeedTrials.Timestamps(HispeedTrials.Lick==1),num2cell(HispeedTrials.Event_Index(HispeedTrials.Lick==1)), 'UniformOutput', false)));

        % Load the SingleCellResponse.mat file, in order to add the responsiveness
        % of the respective burstiness.
        singleCellReactivityName = ['SingleCellReactivity_', condition, '_withTargetEstimate_', num2str(respWindow/1000,'respond%.3f-%.3f_'), num2str(spontWindow/1000,'baseline%.3f-%.3f.mat')];
        load(fullfile(currFile,singleCellReactivityName),'SingleCellResponse')

        for i = 1:numel(trialType)+1
            % Add variables, if it is missing
            if i == numel(trialType)+1
                newVarNames = {[cell2mat(join(trialType,'VS')),'_Bursts'],[cell2mat(join(trialType,'VS')),'_Tonic']};
            elseif i < numel(trialType)+1
                newVarNames = {[trialType{i},'Burstiness'],[trialType{i},'Tonic']};
            end
            idx = ismember(newVarNames, SingleCellResponse.Properties.VariableNames);
            SingleCellResponse = removevars(SingleCellResponse,newVarNames(idx));
            SingleCellResponse = addvars(SingleCellResponse,false(height(SingleCellResponse),1),false(height(SingleCellResponse),1),'NewVariableNames',newVarNames);
        end


        %% Get file information

        % Get burst data from file
        load(fullfile(currFile,'BurstinessData.mat'),'BurstStarts','TypeOfSpike','NumSpikesInBursts','BurstFrequencies','unitIDs')

        dataInfo = dir(fullfile(currFile,'*all_channels.mat'));
        try
            load(fullfile(dataInfo(1).folder,dataInfo(1).name),'sortedData','fs')
            sortedData = sortedData(all(~cellfun(@isempty, sortedData),2),:);
            clInfo = getClusterInfo(fullfile(currFile,'cluster_info.tsv'));
        catch
            fprintf('\nCould not gather sortedData and clInfo.\n')
            return
        end

        str_idx = regexp(currFile,'#\d*','end');
        animalPath = currFile(1:str_idx);
        % Include only those units, which tetrodes were within the recording area
        load(fullfile(animalPath,'targetHit.mat'),'targetHit')

        clInfo.group(cell2mat(sortedData(:,3))==1) = deal({'good'});
        clInfo.group(cell2mat(sortedData(:,3))==2) = deal({'mua'});
        clInfo.group(cell2mat(sortedData(:,3))==3) = deal({'noise'});

        include = ismember(clInfo.group,{'good','mua'}) & clInfo.isolation_distance > 15 & clInfo.isi_viol < 3 & ismember(clInfo.ch,targetHit);
        sortedData = sortedData(include,:);
        
        % Filter for response type
        switch cellSpecs
            case 2 % onlyBurst
                idx = ismember(burstSignTable.SessionName,fileSelection{file}) & ...
                    ismember(burstSignTable.UnitID, sortedData(:,1)) & ...
                    burstSignTable.burstSign;
            case 3 % onlyTonic
                idx = ismember(tonicSignTable.SessionName,fileSelection{file}) & ...
                    ismember(tonicSignTable.UnitID, sortedData(:,1)) & ...
                    tonicSignTable.tonicSign & ~burstSignTable.burstSign;
            case 4 % onlyBoth
                idx = ismember(tonicSignTable.SessionName,fileSelection{file}) & ...
                    ismember(tonicSignTable.UnitID, sortedData(:,1)) & ...
                    tonicSignTable.tonicSign & burstSignTable.burstSign;
            case 5 % onlyNone
                idx = ismember(tonicSignTable.SessionName,fileSelection{file}) & ...
                    ismember(tonicSignTable.UnitID, sortedData(:,1)) & ...
                    ~tonicSignTable.tonicSign & ~burstSignTable.burstSign;
            case 6
                idx = SingleCellResponse.include;
                touchModulated = SingleCellResponse.respWide(idx) | SingleCellResponse.respNarrow(idx) | SingleCellResponse.respInterm(idx);
                sortedData = sortedData(touchModulated,:);
                BurstStarts = BurstStarts(ismember(unitIDs,sortedData(:,1)));
                TypeOfSpike = TypeOfSpike(ismember(unitIDs,sortedData(:,1)));
                NumSpikesInBursts = NumSpikesInBursts(ismember(unitIDs,sortedData(:,1)));
                BurstFrequencies = BurstFrequencies(ismember(unitIDs,sortedData(:,1)));
        end
        if ismember(cellSpecs,(2:5))
                sortedData = sortedData(ismember(sortedData(:,1),burstSignTable.UnitID(idx)),:);
                BurstStarts = BurstStarts(ismember(unitIDs,sortedData(:,1)));
                TypeOfSpike = TypeOfSpike(ismember(unitIDs,sortedData(:,1)));
                NumSpikesInBursts = NumSpikesInBursts(ismember(unitIDs,sortedData(:,1)));
                BurstFrequencies = BurstFrequencies(ismember(unitIDs,sortedData(:,1)));
                put_ex_idx = burstSignTable.PutExcitatory(idx);
        end

        % Label units with their respective area as subscripts
        sortedData(:,4) = deal({'NaN'});
        for i = 1:height(sortedData)
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

        tempDir = currFile;
        rhdFile = dir(fullfile(tempDir,'*.rhd'));
        iter = 1;
        while isempty(rhdFile) && iter <=3
            tempDir = fileparts(tempDir);
            rhdFile = dir(fullfile(tempDir,'*.rhd'));
            iter = iter + 1;
            if isempty(rhdFile) && iter > 3
                error('No rhd file found.')
            end
        end
        intanDir = rhdFile(1).folder;

        FileInfo = dir(fullfile(currFile,sprintf('StimulusResponse_%s*',condition)));
        load(fullfile(FileInfo(1).folder,FileInfo(1).name),'Wide_Narrow_Intermediate','Lick')


        %% Loop through all units
        firingRates_burstSpikes_temp = cell(1,numel(trialType));
        firingRates_burstSpikes_temp(:) = {NaN(height(sortedData),(timeLapse(2)-timeLapse(1))/binSize)};

        firingRates_tonicSpikes_temp = cell(1,numel(trialType));
        firingRates_tonicSpikes_temp(:) = {NaN(height(sortedData),(timeLapse(2)-timeLapse(1))/binSize)};

        firingRates_allSpikes_temp = cell(1,numel(trialType));
        firingRates_allSpikes_temp(:) = {NaN(height(sortedData),(timeLapse(2)-timeLapse(1))/binSize)};

        burstEvents_temp = cell(1,numel(trialType));
        burstEvents_temp(:) = {NaN(height(sortedData),(timeLapse(2)-timeLapse(1))/binSizeBursts)};

        for unit = 1:height(sortedData)
            unit_idx = find(str2double(sortedData{unit,1})==SingleCellResponse.cluster_id);
            % Retrieve spike times
            spikeTimes_msec = sortedData{unit,2}*1000;

            respondDiff_burst = cell(1,numel(trialType));
            respondDiff_tonic = cell(1,numel(trialType));

            fig = figure;
            tiledlayout(3,4)
            sgtitle(sprintf('Unit #%s | %s',sortedData{unit,1},sortedData{unit,4}))

            for i = 1:numel(trialType)
                if ~isequal(trialType{i},'all')
                    if contains(trialType{i},'wide')
                        logApert = Wide_Narrow_Intermediate==1;
                    elseif contains(trialType{i},'narrow')
                        logApert = Wide_Narrow_Intermediate==2;
                    elseif contains(trialType{i},'intermediate')
                        logApert = Wide_Narrow_Intermediate==3;
                    else
                        logApert = true(length(Wide_Narrow_Intermediate),1);
                    end

                    if contains(trialType{i},'no-lick')
                        logLick = Lick==0;
                    elseif contains(trialType{i},'lick')
                        logLick = Lick==1;
                    else
                        logLick = true(length(Lick),1);
                    end
                    logicalFlag = logApert & logLick;
                end

                % Get the Triggers
                % Set warning to error in order to catch it
                s = warning('error', 'MATLAB:load:variableNotFound');
                condInfo = dir(fullfile(fileparts(currFile),'*analysis.mat'));
                if ~isempty(condInfo)
                    try load(fullfile(condInfo.folder,condInfo.name),'Conditions','update')
                        if update < 5
                            [Conditions, ~] = getConditions(fileparts(intanDir));
                        end
                    catch
                        clear update
                        [Conditions, ~] = getConditions(fileparts(intanDir));
                    end
                else
                    [Conditions, ~] = getConditions(fileparts(intanDir));
                end
                % Reset warning
                warning(s);

                % Adjust Triggers accordingly
                idx = find(cellfun(@(x) isequal(x,condition), {Conditions.name}));
                if ~isequal(trialType{i},'all')
                    for m = 1:numel(idx)
                        Conditions(idx(m)).Triggers = sort(Conditions(idx(m)).Triggers);
                        Conditions(idx(m)).Triggers = Conditions(idx(m)).Triggers(logicalFlag,:);
                    end
                end

                % Get trigger onsets
                triggers = sort(Conditions(idx(1)).Triggers(:,1));
                % Convert trigger points from samples to milliseconds
                triggers = 1000.*triggers./fs;

                withinTimeLapse = false(numel(spikeTimes_msec),1);
                LickWithinTimeLapse = false(numel(lickEvents),1);
                trialNum = numel(triggers);
                trialCount = [];
                lickTrialCount = [];
                spikeGroup = [];
                burstNumSpont = zeros(numel(triggers),1);
                burstNumEvok = zeros(numel(triggers),1);
                tonicNumSpont = zeros(numel(triggers),1);
                tonicNumEvok = zeros(numel(triggers),1);
                burstPropSpont = zeros(numel(triggers),1);
                burstPropEvok = zeros(numel(triggers),1);

                burstsWithinWindow = [];
                doubleSpikes = 0;

                for trig = 1:numel(triggers)
                    % Spikes in time lapse window
                    log_idx = (spikeTimes_msec - triggers(trig)) >= timeLapse(1) & (spikeTimes_msec - triggers(trig)) <= timeLapse(2);
                    lick_log_idx = (lickEvents - triggers(trig)) >= timeLapse(1) & (lickEvents - triggers(trig)) <= timeLapse(2);
                    
                    % With bigger time lapse window sizes, spikes can
                    % occasionally belong to two different trials. Count
                    % those spikes here. Spikes are assigned to the
                    % subsequent trial.
                    if sum(withinTimeLapse(log_idx))
                        doubleSpikes = doubleSpikes + sum(withinTimeLapse(log_idx));
                        trialCount(end-sum(withinTimeLapse(log_idx))+1) = trig; %#ok<SAGROW> 
                        trialCount = [trialCount; repmat(trig,sum(log_idx)-sum(withinTimeLapse(log_idx)),1)]; %#ok<AGROW>
                        spike_type = TypeOfSpike{1,unit}(log_idx);
                        spikeGroup = [spikeGroup(1:end-sum(withinTimeLapse(log_idx))); spike_type'];
                    else
                        trialCount = [trialCount; repmat(trig,sum(log_idx),1)]; %#ok<AGROW>
                        spike_type = TypeOfSpike{1,unit}(log_idx);
                        spikeGroup = [spikeGroup; spike_type']; %#ok<AGROW>
                    end
                    withinTimeLapse(log_idx) = true;

                    if sum(LickWithinTimeLapse(lick_log_idx))
                        lickTrialCount(end-sum(LickWithinTimeLapse(lick_log_idx))+1) = trig; %#ok<SAGROW> 
                        lickTrialCount = [lickTrialCount; repmat(trig,sum(lick_log_idx)-sum(LickWithinTimeLapse(lick_log_idx)),1)]; %#ok<AGROW>
                    else
                        lickTrialCount = [lickTrialCount; repmat(trig,sum(lick_log_idx),1)]; %#ok<AGROW>
                    end
                    LickWithinTimeLapse(lick_log_idx) = true;

                    burst_idx = (BurstStarts{unit} - triggers(trig)) >= timeLapse(1) & (BurstStarts{unit} - triggers(trig)) <= timeLapse(2);
                    burstsWithinWindow = [burstsWithinWindow;milliseconds(BurstStarts{unit}(burst_idx))-milliseconds(triggers(trig))]; %#ok<AGROW>

                    % Get the amount of bursts in the spontaneous window
                    burst_idx = (BurstStarts{unit} - triggers(trig)) >= spontWindow(1) & (BurstStarts{unit} - triggers(trig)) < spontWindow(2);
                    burstNumSpont(trig) = sum(burst_idx);
                    % Get the proportion of burst spikes from all spikes in the spontaneous window
                    log_idx = (spikeTimes_msec - triggers(trig)) >= spontWindow(1) & (spikeTimes_msec - triggers(trig)) < spontWindow(2);
                    burstPropSpont(trig) = max([0 sum(TypeOfSpike{unit}(log_idx))/numel(TypeOfSpike{unit}(log_idx))]);

                    % Get the amount of bursts in the evoked window
                    burst_idx = (BurstStarts{unit} - triggers(trig)) >= respWindow(1) & (BurstStarts{unit} - triggers(trig)) < respWindow(2);
                    burstNumEvok(trig) = sum(burst_idx);
                    % Get the proportion of burst spikes from all spikes in the evoked window
                    log_idx = (spikeTimes_msec - triggers(trig)) >= respWindow(1) & (spikeTimes_msec - triggers(trig)) < respWindow(2);
                    burstPropEvok(trig) = max([0 sum(TypeOfSpike{unit}(log_idx))/numel(TypeOfSpike{unit}(log_idx))]);

                    % Get the proportion of tonic spikes in the spontaneous window
                    log_idx = (spikeTimes_msec - triggers(trig)) >= spontWindow(1) & (spikeTimes_msec - triggers(trig)) < spontWindow(2);
                    tonicNumSpont(trig) = sum(TypeOfSpike{unit}(log_idx)==0);
                    % Get the proportion of tonic spikes in the evoked window
                    log_idx = (spikeTimes_msec - triggers(trig)) >= respWindow(1) & (spikeTimes_msec - triggers(trig)) < respWindow(2);
                    tonicNumEvok(trig) = sum(TypeOfSpike{unit}(log_idx)==0);
                end
                
                if doubleSpikes
                    fprintf(2,'\nUnit #%s: There were %i spike/s that occur in two trials (spikes assigned to latter one). Consider reducing the time window.\n',sortedData{unit,1},doubleSpikes);
                end
                
                if plotLickTimes
                    spikeGroup = categorical(spikeGroup,[0 1],{'Tonic','Burst'});
                    % Concatenate 'Lick' category with the existing array
                    spikeGroup = [spikeGroup; repmat(categorical({'Lick'}), numel(lickEvents(LickWithinTimeLapse)), 1)]; %#ok<AGROW> 

                    timeStamps = [milliseconds(spikeTimes_msec(withinTimeLapse));milliseconds(lickEvents(LickWithinTimeLapse))];
                    trialCount = [trialCount;lickTrialCount]; %#ok<AGROW> 
                else
                    spikeGroup = categorical(spikeGroup,[0 1],{'Tonic','Burst'});
                    timeStamps = milliseconds(spikeTimes_msec(withinTimeLapse));
                end

                % Plot raster for each trial type
                nexttile(i,[2 1])
                s = spikeRasterPlot(timeStamps,...
                    trialCount, 'AllTrialsCount', cellfun(@num2str, num2cell(1:numel(triggers))', 'UniformOutput', false),...
                    'GroupData', spikeGroup,'ColorOrder',[0 0 0;[255 71 26]/255;[26 255 255]/255]);
                xlim(milliseconds(timeLapse))
                s.AlignmentTimes = milliseconds(triggers(unique(trialCount)));
                s.LegendVisible = false;
                s.TitleText = trialType{i};
                yticklabels = {''};
                if i==1
                    ylabel('Trials')
                end

                nexttile(i+8)

                % Extract spikes in the respecticve time window
                spikesAligned = sort(s.SpikeTimeData-s.AlignmentTimes(s.TrialData));
                burstSpikesAligned = sort(s.SpikeTimeData(s.GroupData=='Burst')-s.AlignmentTimes(s.TrialData(s.GroupData=='Burst')));
                tonicSpikesAligned = sort(s.SpikeTimeData(s.GroupData=='Tonic')-s.AlignmentTimes(s.TrialData(s.GroupData=='Tonic')));

                firingRates = (histcounts(spikesAligned,binEdgesFiring)/trialNum)/(binSize/1000);
                firingRates_allSpikes_temp{i}(unit,:) = firingRates;

                firingRates = (histcounts(burstSpikesAligned,binEdgesFiring)/trialNum)/(binSize/1000);
                firingRates_burstSpikes_temp{i}(unit,:) = firingRates;

                firingRates = (histcounts(tonicSpikesAligned,binEdgesFiring)/trialNum)/(binSize/1000);
                firingRates_tonicSpikes_temp{i}(unit,:) = firingRates;

                % Plot mean firing rates of that unit underneath the raster plots
                yyaxis left
                set(gca,'YColor',[0.15 0.15 0.15]);
                plot(binEdgesFiring(1:end-1) + diff(binEdgesFiring)/2,movmean(firingRates_allSpikes_temp{i}(unit,:),2),'k')
                xlim tight
                if i==1
                    ylabel('Firing rate [Hz]','Color','k')
                end
                
                if isempty(burstsWithinWindow)
                    burstEvents_temp{i}(unit,:) = 0;
                else
                    burstEvents_temp{i}(unit,:) = (histcounts(burstsWithinWindow,binEdgesBursts)/trialNum)/(binSizeBursts /1000);
                end

                yyaxis right
                set(gca,'YColor',[255 71 26]/255);
                bar(binEdgesBursts(1:end-1) + diff(binEdgesBursts)/2,burstEvents_temp{i}(unit,:),'FaceColor','None','EdgeColor',[255 71 26]/255,'BarWidth',1)
                if i==numel(trialType)
                    ylabel('Burst event rate [Hz]','Color',[255 71 26]/255)
                end

                xlabel('Time [sec]')
                ax = gca;
                ax.XTickLabel = cellfun(@num2str, num2cell(timeLapse(1)/1000:0.2:timeLapse(2)/1000)', 'UniformOutput', false);

                % Dot dash plot with mean of baseline and evoked window burst proportions
                % burstIndex [-1;1] = burstProportion_post - burstProportion_pre
                nexttile((i-1)*4+3)
                title(trialType{i})
                hold on

                respondDiff_burst{i} = burstPropEvok-burstPropSpont;
                % Thicken the lines, of stronger represented data pairs
                uniquePairs = unique([burstPropSpont,burstPropEvok],'rows');
                for pair = 1:height(uniquePairs)
                    numPairs = sum(ismember([burstPropSpont,burstPropEvok],uniquePairs(pair,:),'rows'));
                    % Adjust line thickness from 0.5 (default) to 2.5, depending on
                    % how often the data pair is represented.
                    thickness = 0.5+round(2*numPairs/numel(triggers),1);
                    line([1,2],uniquePairs(pair,:),'LineWidth',thickness, 'Color', 0.7*[1 1 1])
                end
                % Plot mean burst number
                line([0.9,1.1],[mean(burstPropSpont,'omitnan') mean(burstPropSpont,'omitnan')],'LineWidth',1.5, 'Color', [0 0 0])
                line([1.9,2.1],[mean(burstPropEvok,'omitnan') mean(burstPropEvok,'omitnan')],'LineWidth',1.5, 'Color', [0 0 0])

                xlim([0.75 2.25])
                ylim([-0.1 1.1])
                ylabel('Burst spike prop')
                xticks([1,2])
                xticklabels({'spont','evoked'})

                if ~isempty(burstPropSpont) && ~isempty(burstPropEvok)
                    p = signrank(burstPropSpont,burstPropEvok);

                    yt = get(gca, 'YTick');
                    xt = get(gca, 'XTick');
                    if p < 0.001
                        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '***','HorizontalAlignment','center')
                        SingleCellResponse{unit_idx,[trialType{i},'Burstiness']} = true;
                    elseif p < 0.01
                        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '**','HorizontalAlignment','center')
                        SingleCellResponse{unit_idx,[trialType{i},'Burstiness']} = true;
                    elseif p < 0.05
                        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '*','HorizontalAlignment','center')
                        SingleCellResponse{unit_idx,[trialType{i},'Burstiness']} = true;
                    else
                        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.9, 'n.s.','HorizontalAlignment','center')
                    end
                end

                % Dot dash plot with mean of baseline and evoked window tonic proportions
                % tonicIndex [0;1] = tonicSpikes_post / (tonicSpikes_post + tonicSpikes_pre)
                nexttile((i-1)*4+4)
                title(trialType{i})
                hold on

                respondDiff_tonic{i} = tonicNumEvok./(tonicNumEvok+tonicNumSpont);
                % Thicken the lines, of stronger represented data pairs
                uniquePairs = unique([tonicNumSpont,tonicNumEvok],'rows');
                for pair = 1:height(uniquePairs)
                    numPairs = sum(ismember([tonicNumSpont,tonicNumEvok],uniquePairs(pair,:),'rows'));
                    % Adjust line thickness from 0.5 (default) to 2.5, depending on
                    % how often the data pair is represented.
                    thickness = 0.5+round(2*numPairs/numel(triggers),1);
                    line([1,2],uniquePairs(pair,:),'LineWidth',thickness, 'Color', 0.7*[1 1 1])
                end
                % Plot mean tonic number
                line([0.9,1.1],[mean(tonicNumSpont,'omitnan') mean(tonicNumSpont,'omitnan')],'LineWidth',1.5, 'Color', [0 0 0])
                line([1.9,2.1],[mean(tonicNumEvok,'omitnan') mean(tonicNumEvok,'omitnan')],'LineWidth',1.5, 'Color', [0 0 0])

                xlim([0.75 2.25])
                ylabel('Tonic spikes')
                xticks([1,2])
                xticklabels({'spont','evoked'})

                if ~isempty(tonicNumSpont) && ~isempty(tonicNumEvok)
                    p = signrank(tonicNumSpont,tonicNumEvok);

                    yt = get(gca, 'YTick');
                    xt = get(gca, 'XTick');
                    if p < 0.001
                        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '***','HorizontalAlignment','center')
                        SingleCellResponse{unit_idx,[trialType{i},'Tonic']} = true;
                    elseif p < 0.01
                        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '**','HorizontalAlignment','center')
                        SingleCellResponse{unit_idx,[trialType{i},'Tonic']} = true;
                    elseif p < 0.05
                        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '*','HorizontalAlignment','center')
                        SingleCellResponse{unit_idx,[trialType{i},'Tonic']} = true;
                    else
                        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.9, 'n.s.','HorizontalAlignment','center')
                    end
                end
            end

            % Compare response magnitudes for each trial type for burst spikes
            % burstIndex [-1;1] = burstProportion_post - burstProportion_pre
            % burstPreference [-1;1] = 0.5*(burstIndex_wide - burstIndex_narrow)
            nexttile(11)
            title(join(trialType,' vs. '))
            hold on

            xlim([0.75 2.25])
            ylim padded
            ylabel('Burst diff')
            xticks([1,2])
            xticklabels(trialType)

            if ~any(cellfun(@isempty, respondDiff_burst))
                for i = 1:2
                    % Thicken the markers of stronger represented data
                    uniqueDiff = unique(respondDiff_burst{i});
                    for pair = 1:height(uniqueDiff)
                        numPairs = sum(respondDiff_burst{i}==uniqueDiff(pair));
                        % Adjust marker thickness from 10 (default) to 30, depending on
                        % how often the data pair is represented.
                        thickness = 10+round(20*numPairs/numel(respondDiff_burst{i}));
                        plot(i,uniqueDiff(pair),'.', 'MarkerSize', thickness, 'Color', 0.7*[1 1 1]);
                    end
                    % Plot mean burst number
                    line([i-0.1 i+0.1],[mean(respondDiff_burst{i},'omitnan') mean(respondDiff_burst{i},'omitnan')],'LineWidth',1.5, 'Color', [0 0 0])
                end
                p = ranksum(respondDiff_burst{1},respondDiff_burst{2});
                yt = get(gca, 'YTick');
                xt = get(gca, 'XTick');
                if p < 0.001
                    plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '***','HorizontalAlignment','center')
                    SingleCellResponse{unit_idx,[cell2mat(join(trialType,'VS')),'_Bursts']} = true;
                elseif p < 0.01
                    plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '**','HorizontalAlignment','center')
                    SingleCellResponse{unit_idx,[cell2mat(join(trialType,'VS')),'_Bursts']} = true;
                elseif p < 0.05
                    plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '*','HorizontalAlignment','center')
                    SingleCellResponse{unit_idx,[cell2mat(join(trialType,'VS')),'_Bursts']} = true;
                else
                    plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.9, 'n.s.','HorizontalAlignment','center')
                end
            end

            % Compare response magnitudes for each trial type for tonic spikes
            % tonicIndex [0;1] = tonicSpikes_post / (tonicSpikes_post + tonicSpikes_pre)
            % tonicPreference [-1;1] = tonicIndex_wide - tonicIndex_narrow
            nexttile(12)
            title(join(trialType,' vs. '))
            hold on

            xlim([0.75 2.25])
            ylim padded
            ylabel('Tonic diff')
            xticks([1,2])
            xticklabels(trialType)

            if ~any(cellfun(@isempty, respondDiff_tonic)) && ~any(cellfun(@(x) all(isnan(x)), respondDiff_tonic))
                for i = 1:2
                    % Thicken the markers of stronger represented data
                    uniqueDiff = unique(respondDiff_tonic{i});
                    for pair = 1:height(uniqueDiff)
                        numPairs = sum(respondDiff_tonic{i}==uniqueDiff(pair));
                        % Adjust marker thickness from 10 (default) to 30, depending on
                        % how often the data pair is represented.
                        thickness = 10+round(20*numPairs/numel(respondDiff_tonic{i}));
                        plot(i,uniqueDiff(pair),'.', 'MarkerSize', thickness, 'Color', 0.7*[1 1 1]);
                    end
                    % Plot mean tonic number
                    line([i-0.1 i+0.1],[mean(respondDiff_tonic{i},'omitnan') mean(respondDiff_tonic{i},'omitnan')],'LineWidth',1.5, 'Color', [0 0 0])
                end
                p = ranksum(respondDiff_tonic{1},respondDiff_tonic{2});
                yt = get(gca, 'YTick');
                xt = get(gca, 'XTick');
                if p < 0.001
                    plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '***','HorizontalAlignment','center')
                    SingleCellResponse{unit_idx,[cell2mat(join(trialType,'VS')),'_Tonic']} = true;
                elseif p < 0.01
                    plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '**','HorizontalAlignment','center')
                    SingleCellResponse{unit_idx,[cell2mat(join(trialType,'VS')),'_Tonic']} = true;
                elseif p < 0.05
                    plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '*','HorizontalAlignment','center')
                    SingleCellResponse{unit_idx,[cell2mat(join(trialType,'VS')),'_Tonic']} = true;
                else
                    plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.9, 'n.s.','HorizontalAlignment','center')
                end
            end
            % if p > 0.05  % To plot only tonic responding
            %     close(fig)
            % end
            % close all % close figure, in order to save space
        end
        
        if cellSpecs==1 % only save when all cells are analyzed
            % Save SingleCellResponse.mat file with added burstiness data
            save(fullfile(currFile,singleCellReactivityName),'SingleCellResponse')
        end

        WaveformInfo = dir(fullfile(dataInfo(1).folder, '*waveforms_all.mat'));
        try load(fullfile(WaveformInfo.folder,WaveformInfo.name),'clWaveforms')
        catch
        end
        
        if cellSpecs==1
            if all(cellfun(@(x,y) isequal(x,y), sortedData(:,1),clWaveforms(:,1)))
                [~, put_ex_idx] = getWaveformDistribution(clWaveforms,sortedData(:,4));
            elseif all(cellfun(@(x) any(ismember(sortedData(:,1),x)), clWaveforms(:,1)))
                % Order the units according to the ClusterBurstiness table
                unitOrder = nan(height(clWaveforms),1);
                for ii = 1:height(sortedData)
                    unitOrder(ii) = find(ismember(clWaveforms(:,1),sortedData{ii,1}));
                end
                clWaveforms = clWaveforms(unitOrder,:);

                [~, put_ex_idx] = getWaveformDistribution(clWaveforms,sortedData(:,4));
            else
                error('burstUponTrigger:waveformUnitsError','Units in the clWaveforms variable do not match the ones from the sortedData.')
            end
        elseif cellSpecs==6 % Touch-modulated cells
            idx = cellfun(@(x) any(ismember(sortedData(:,1),x)), clWaveforms(:,1));
            clWaveforms = clWaveforms(idx,:);
            % Order the units according to the ClusterBurstiness table
            unitOrder = nan(height(clWaveforms),1);
            for ii = 1:height(sortedData)
                unitOrder(ii) = find(ismember(clWaveforms(:,1),sortedData{ii,1}));
            end
            clWaveforms = clWaveforms(unitOrder,:);
            
            [~, put_ex_idx] = getWaveformDistribution(clWaveforms,sortedData(:,4));
        end

        unitInfo = table(sortedData(:,1),sortedData(:,4),put_ex_idx,...
            'VariableNames',{'ID','Area','putExcit'});

        for i = 1:numel(trialType)
            singleCellFiringRateName = ['SingleCellFiringRates_', condition,num2str(timeLapse/1000,'_timeLapse%.3f-%.3f_'),...
                num2str(respWindow/1000,'respWind%.3f-%.3f_'), num2str(spontWindow/1000,'spontWind%.3f-%.3f_'),...
                num2str(binSize,'%imsecbinSize_'),num2str(binSizeBursts,'%imsecbinSizeBursts_'),trialType{i},'_',cellDescript,'.mat'];

            firingRates_burstSpikes = firingRates_burstSpikes_temp{i};
            firingRates_tonicSpikes = firingRates_tonicSpikes_temp{i};
            firingRates_allSpikes = firingRates_allSpikes_temp{i};
            burstEvents = burstEvents_temp{i};
            % Save singleCellFiringRate.mat file with firing rate information
            % Change this to track changes and automatically re-run the script, if files are not up to date
            fileVersion = 2; 
            save(fullfile(currFile,singleCellFiringRateName),...
                'firingRates_burstSpikes','firingRates_tonicSpikes','firingRates_allSpikes','burstEvents','unitInfo','fileVersion')
        end
    end
end

%% Load all relevant data and concatenate them

firingRates_burstSpikes_allFiles = cell(1,height(fileSelection));
firingRates_burstSpikes_allFiles(:) = {cell(1,numel(trialType))};

firingRates_tonicSpikes_allFiles = cell(1,height(fileSelection));
firingRates_tonicSpikes_allFiles(:) = {cell(1,numel(trialType))};

firingRates_allSpikes_allFiles = cell(1,height(fileSelection));
firingRates_allSpikes_allFiles(:) = {cell(1,numel(trialType))};

burstEvents_allFiles = cell(1,height(fileSelection));
burstEvents_allFiles(:) = {cell(1,numel(trialType))};

UnitIDs_allFiles = {};
UnitAreas_allFiles = {};
putExcit_allFiles = logical([]);

for file = 1:height(fileSelection)
    currFile = fileSelection{file};
    % Load respective variables for each trial type
    for i = 1:numel(trialType)
        singleCellFiringRateName = ['SingleCellFiringRates_', condition,num2str(timeLapse/1000,'_timeLapse%.3f-%.3f_'),...
                num2str(respWindow/1000,'respWind%.3f-%.3f_'), num2str(spontWindow/1000,'spontWind%.3f-%.3f_'),...
                num2str(binSize,'%imsecbinSize_'),num2str(binSizeBursts,'%imsecbinSizeBursts_'),trialType{i},'_',cellDescript,'.mat'];
        load(fullfile(currFile,singleCellFiringRateName),...
            'firingRates_burstSpikes','firingRates_tonicSpikes','firingRates_allSpikes','burstEvents','unitInfo')
        firingRates_burstSpikes_allFiles{file}{i} = firingRates_burstSpikes;
        firingRates_tonicSpikes_allFiles{file}{i} = firingRates_tonicSpikes;
        firingRates_allSpikes_allFiles{file}{i} = firingRates_allSpikes;
        burstEvents_allFiles{file}{i} = burstEvents;
    end
    UnitIDs_allFiles = [UnitIDs_allFiles;unitInfo.ID]; %#ok<AGROW>
    UnitAreas_allFiles = [UnitAreas_allFiles;unitInfo.Area]; %#ok<AGROW>
    putExcit_allFiles = [putExcit_allFiles;unitInfo.putExcit]; %#ok<AGROW>
end

%% Choose which unit waveforms to plot
set(0,'DefaultFigureVisible','on');

% Choose which waveforms to analyze
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
    'Callback',{@doneExeWave,f,ui_field,area_names});

waitfor(doneButton)

%% Plot average area activity
firingRatesAll_allSpikes = cell(1,numel(trialType));
firingRatesAll_burstSpikes = cell(1,numel(trialType));
firingRatesAll_tonicSpikes = cell(1,numel(trialType));

% Array with only the start times of each burst
burstActivityAll = cell(1,numel(trialType));

% Loop through the main cell array XX
figBurstSpikesProp = gobjects(1, numel(trialType));
for i = 1:numel(trialType)
    % Concatenate the subcells using cell2mat
    firingRatesAll_allSpikes{i} = cell2mat(cellfun(@(x) cell2mat(x(i)), firingRates_allSpikes_allFiles, 'UniformOutput', false)');
    firingRatesAll_burstSpikes{i} = cell2mat(cellfun(@(x) cell2mat(x(i)), firingRates_burstSpikes_allFiles, 'UniformOutput', false)');
    firingRatesAll_tonicSpikes{i} = cell2mat(cellfun(@(x) cell2mat(x(i)), firingRates_tonicSpikes_allFiles, 'UniformOutput', false)');

    burstActivityAll{i} = cell2mat(cellfun(@(x) cell2mat(x(i)), burstEvents_allFiles, 'UniformOutput', false)');
    
    % Preallocate burst proportion figures for plotting
    figBurstSpikesProp(i) = figure('Name',sprintf('BurstSpikesProp_allAreas_%s',trialType{i}));
    if binSize == binSizeBursts % Only generate figure, if the burst events are binned in the same manner
        figBurstEventProp(i) = figure('Name',sprintf('BurstEventProp_allAreas_%s',trialType{i})); %#ok<SAGROW>
    end
end

% Comparison figure for respective trial type burst firing rates
% Significance tests will be performed every individual data point
figTrialComp_indi = figure('Name','BurstRateComp_allAreas_individualPoints');

% Significance tests will be performed on time windows rather than
% individual data points
windowLength = 50; % in msec
localBonferroni = 50; % local correction for multiple testing, set to 'off' to correct over the whole time lapse
windowBins = diff(timeLapse)/windowLength;
if isnumeric(localBonferroni)
    bonferroniString = sprintf('local bonferroni %imsec, i.e., %.1f bins',localBonferroni,localBonferroni/windowLength);
elseif isequal(localBonferroni, 'off')
    bonferroniString = sprintf('bonferroni over whole time window, i.e., %i bins',windowBins);
    localBonferroni = diff(timeLapse);
end
figTrialBurstComp_wind = figure('Name',sprintf('BurstRateComp_allAreas_%imsecWindow',windowLength));
figTrialTonicComp_wind = figure('Name',sprintf('TonicRateComp_allAreas_%imsecWindow',windowLength));
figTrialAllComp_wind = figure('Name',sprintf('FiringRateComp_allAreas_%imsecWindow',windowLength));

% Figures for comparing the areas for each trial type
axAreaComp_abs = gobjects(1,numel(trialType));
axAreaComp_norm = gobjects(1,numel(trialType));
axAreaComp_smooth = gobjects(1,numel(trialType));

for i = 1:numel(trialType)
    fig = figure('Name',sprintf('AreaComp_absolute_%s',trialType{i}));
    axAreaComp_abs(i) = axes(fig); %#ok<LAXES> 
    title(sprintf('Area burst activtiy - %s',trialType{i}))
    ylabel('Burst spike firing rate [Hz]')
    
    fig = figure('Name',sprintf('AreaComp_norm_%s',trialType{i}));
    axAreaComp_norm(i) = axes(fig); %#ok<LAXES> 
    title(sprintf('Normalized area burst activtiy - %s',trialType{i}))
    ylabel('Norm burst spike firing rate')

    fig = figure('Name',sprintf('AreaComp_smooth_%s',trialType{i}));
    smoothWind = 8; % Number of bin sizes to smooth over
    fig.UserData = struct('numSmoothWind', smoothWind);
    axAreaComp_smooth(i) = axes(fig); %#ok<LAXES> 
    title(sprintf('Normalized smoothed area burst activtiy - %s',trialType{i}))
    ylabel('Norm burst spike firing rate')
end

% Preallocate burst plots and cell arrays
p_AreaComp_abs = gobjects(numel(trialType),numel(area_names));
p_AreaComp_norm = gobjects(numel(trialType),numel(area_names));
p_AreaComp_smooth = gobjects(numel(trialType),numel(area_names));
for ar = 1:numel(area_names)
    % Plot burst events (i.e., starting time of each burst) of each unit as colored PSTH
    % and underneath the population mean firing rate (all spikes) with an
    % overlay of the respective mean burst event rates
    if isequal(area_waveforms{ar},'All')
        fig = figure('Name',sprintf('%s_BurstEventRate',area_names{ar}));
    else
        fig = figure('Name',sprintf('%s_BurstEventRate_%sunits',area_names{ar},area_waveforms{ar}));
    end
    fig.UserData = struct('normMethodBursts', 'z-score');

    tiledlayout(3,2)
    sgtitle(sprintf('%s burst normalized activity',area_names{ar}))
    fr_all_area = cell(1,numel(trialType));
    p_TrialCompAll = gobjects(1,numel(trialType));    
    for i = 1:numel(trialType)
        % Subplot the mean firing rate for each trial type and area
        if i==1
            lineStyle = '-';
        else
            lineStyle = '--';
        end
        
        nexttile(i,[2 1])
        
        % Filter respective units
        if isequal(area_waveforms{ar},'RS')
            idx = ismember(UnitAreas_allFiles,area_names{ar}) & putExcit_allFiles;
        elseif isequal(area_waveforms{ar},'FS')
            idx = ismember(UnitAreas_allFiles,area_names{ar}) & ~putExcit_allFiles;
        else
            idx = ismember(UnitAreas_allFiles,area_names{ar});
        end

        fr_all_area{i} = firingRatesAll_allSpikes{i}(idx,:);

        burst_area = burstActivityAll{i}(idx,:);
        burst_area = zscore(burst_area,0,2);
        %         burst_area = (burst_area-min(burst_area,[],2))./(max(burst_area,[],2)-min(burst_area,[],2));

        imagesc(burst_area);
        title(trialType{i})
        xticklabels([])
        ylabel('Units')
        set(gca, 'YTick', 1:height(burst_area), 'YTickLabel', UnitIDs_allFiles(idx))
        % colormap("gray")
        colormap(magma)

        nexttile(i+4)
        hold on
        yyaxis left
        set(gca,'YColor',[0.15 0.15 0.15]);
        if i==1
            ylabel('Mean firing rate [Hz]','Color','k')
        end

        xvals = binEdgesFiring(1:end-1) + diff(binEdgesFiring)/2;
        % Plot mean with standard error (SEM)
        curve1 = mean(fr_all_area{i},1,'omitnan') + std(fr_all_area{i},[],1,'omitnan')/sqrt(width(fr_all_area{i}));
        curve2 = mean(fr_all_area{i},1,'omitnan') - std(fr_all_area{i},[],1,'omitnan')/sqrt(width(fr_all_area{i}));
        plot(xvals, curve1,'-k');
        plot(xvals, curve2,'-k');
        
        idx = find(~isnan(curve1));
        fill([xvals(idx) fliplr(xvals(idx))], [curve1(idx) fliplr(curve2(idx))],[0 0 .85],...
            'FaceColor','#000000','EdgeColor','none','FaceAlpha',0.1);
        plot(xvals,mean(fr_all_area{i},1,'omitnan'),'-k','LineWidth',1.5)
        xlim tight

        yyaxis right
        set(gca,'YColor',[255 71 26]/255);
        if i==numel(trialType)
            ylabel('Burst rate [z-score]','Color',[255 71 26]/255)
        end
        bar(binEdgesBursts(1:end-1) + diff(binEdgesBursts)/2,mean(burst_area,1,'omitnan'),'FaceColor','None','EdgeColor',[255 71 26]/255,'BarWidth',1)

        % Plot into comparison figure
        subAxAll = subplot(2, 2, ar, 'Parent', figTrialAllComp_wind);
        hold(subAxAll,'on')
        
        plot(subAxAll,xvals, curve1,lineStyle,'Color',area_colors{ar});
        plot(subAxAll,xvals, curve2,lineStyle,'Color',area_colors{ar});

        fill(subAxAll,[xvals(idx) fliplr(xvals(idx))], [curve1(idx) fliplr(curve2(idx))],[0 0 .85],...
            'FaceColor',area_colors{ar},'EdgeColor','none','FaceAlpha',0.1);
        p_TrialCompAll(i) = plot(subAxAll,xvals,mean(fr_all_area{i},1,'omitnan'),lineStyle,'Color',area_colors{ar},'LineWidth',1.5);

        title(subAxAll,sprintf('%s',area_names{ar}))
        ylabel(subAxAll,'Mean firing rate [Hz]')

        if ar==1 && i==numel(trialType)
            axLegAll = legend(subAxAll,[p_TrialCompAll(:)],trialType,'Location','northwest','AutoUpdate','off');
        end
        hold(subAxAll,'off')
    end

    % Plot overall firing rates of each unit as colored PSTH and underneath
    % the population mean firing rate for tonic and burst spikes

    if isequal(area_waveforms{ar},'All')
        fig = figure('Name',sprintf('%s_BurstTonicFR',area_names{ar}));
    else
        fig = figure('Name',sprintf('%s_BurstTonicFR_%sunits',area_names{ar},area_waveforms{ar}));
    end
    fig.UserData = struct('normMethodPSTH', 'MinMax');

    tiledlayout(3,2)
    sgtitle(sprintf('%s burst and tonic firing rate',area_names{ar}))
    
    % Preallocate burst plots and cell arrays
    p_TrialCompBurst = gobjects(1,numel(trialType));
    p_TrialCompTonic = gobjects(1,numel(trialType));
    fr_burst_area = cell(1,numel(trialType));
    fr_tonic_area = cell(1,numel(trialType));
    for i = 1:numel(trialType)
        nexttile(i,[2 1])

        % Filter respective units
        if isequal(area_waveforms{ar},'RS')
            idx = ismember(UnitAreas_allFiles,area_names{ar}) & putExcit_allFiles;
        elseif isequal(area_waveforms{ar},'FS')
            idx = ismember(UnitAreas_allFiles,area_names{ar}) & ~putExcit_allFiles;
        else
            idx = ismember(UnitAreas_allFiles,area_names{ar});
        end

        % fr_all_area_norm = zscore(fr_all_area{i},0,2);
        fr_all_area_norm = (fr_all_area{i}-min(fr_all_area{i},[],2))./(max(fr_all_area{i},[],2)-min(fr_all_area{i},[],2));
        
        % Sort the responses in descending order regarding their reponse magnitude
        spontBins = find(binEdgesFiring >= milliseconds(spontWindow(1)) & binEdgesFiring < milliseconds(spontWindow(2)));
        respBins = find(binEdgesFiring >= milliseconds(respWindow(1)) & binEdgesFiring < milliseconds(respWindow(2)));
        respMagnitude = arrayfun(@(x) mean(fr_all_area_norm(x,respBins))-mean(fr_all_area_norm(x,spontBins)),1:height(fr_all_area_norm));
        [~, sortIdx] = sort(respMagnitude);
        % Print amount of depressed and excited units
        fprintf('No. of excited units in %s (condition: %s): %i (%.2f%%)\n', ...
            area_names{ar}, trialType{i},sum(respMagnitude>=0), 100*sum(respMagnitude>=0)/numel(respMagnitude))
        fprintf('No. of depressed units in %s (condition: %s): %i (%.2f%%)\n\n', ...
            area_names{ar}, trialType{i},sum(respMagnitude<0), 100*sum(respMagnitude<0)/numel(respMagnitude))
        

        fr_burst_area{i} = firingRatesAll_burstSpikes{i}(idx,:);
        fr_tonic_area{i} = firingRatesAll_tonicSpikes{i}(idx,:);

        imagesc(fr_all_area_norm(sortIdx,:));
        title(trialType{i})
        xticklabels([])
        ylabel('Units')
        unitNames = UnitIDs_allFiles(idx);
        set(gca, 'YTick', 1:height(burst_area), 'YTickLabel', unitNames(sortIdx))
        colormap(magma)

        nexttile(i+4)
        hold on
        yyaxis left
        set(gca,'YColor',[0.15 0.15 0.15]);
        if i==1
            ylabel('Mean tonic firing rate [Hz]','Color','k')
        end

        xvals = binEdgesFiring(1:end-1) + diff(binEdgesFiring)/2;
        % Plot mean with standard error (SEM)
        curve1 = mean(fr_tonic_area{i},1,'omitnan') + std(fr_tonic_area{i},[],1,'omitnan')/sqrt(width(fr_tonic_area{i}));
        curve2 = mean(fr_tonic_area{i},1,'omitnan') - std(fr_tonic_area{i},[],1,'omitnan')/sqrt(width(fr_tonic_area{i}));
        plot(xvals, curve1,'-k');
        plot(xvals, curve2,'-k');
        x = xvals(~isnan(curve1));
        curve1 = curve1(~isnan(curve1));
        curve2 = curve2(~isnan(curve2));
        fill([x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
            'FaceColor','#000000','EdgeColor','none','FaceAlpha',0.1);
        plot(xvals,mean(fr_tonic_area{i},1,'omitnan'),'-k','LineWidth',1.5)
        xlim tight

        yyaxis right
        set(gca,'YColor',[255 71 26]/255);
        if i==numel(trialType)
            ylabel('Mean burst firing rate [Hz]','Color',[255 71 26]/255)
        end        

        % Plot mean with standard error (SEM)
        curve1 = mean(fr_burst_area{i},1,'omitnan') + std(fr_burst_area{i},[],1,'omitnan')/sqrt(width(fr_burst_area{i}));
        curve2 = mean(fr_burst_area{i},1,'omitnan') - std(fr_burst_area{i},[],1,'omitnan')/sqrt(width(fr_burst_area{i}));
        plot(xvals, curve1,'-','Color',[255 71 26]/255);
        plot(xvals, curve2,'-','Color',[255 71 26]/255);
        x = xvals(~isnan(curve1));
        curve1 = curve1(~isnan(curve1));
        curve2 = curve2(~isnan(curve2));
        fill([x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
            'FaceColor',[255 71 26]/255,'EdgeColor','none','FaceAlpha',0.1);
        plot(xvals,mean(fr_burst_area{i},1,'omitnan'),'-','Color',[255 71 26]/255,'LineWidth',1.5)
        xlim tight
        
        % Subplot the burst firing rate for each trial type and area
        if i==1
            lineStyle = '-';
        else
            lineStyle = '--';
        end
        subAxBurst = subplot(2, 2, ar, 'Parent', figTrialComp_indi);
        subAxTonic = subplot(2, 2, ar, 'Parent', figTrialTonicComp_wind);
        hold(subAxBurst,'on')
        hold(subAxTonic,'on')
        hold(axAreaComp_abs(i),'on')
        hold(axAreaComp_norm(i),'on')
        hold(axAreaComp_smooth(i),'on')

        % Plot burst spike rates
        curve1 = mean(fr_burst_area{i},1,'omitnan') + std(fr_burst_area{i},[],1,'omitnan')/sqrt(width(fr_burst_area{i}));
        curve2 = mean(fr_burst_area{i},1,'omitnan') - std(fr_burst_area{i},[],1,'omitnan')/sqrt(width(fr_burst_area{i}));
        
        plot(subAxBurst,xvals, curve1,lineStyle,'Color',area_colors{ar});
        plot(axAreaComp_abs(i),xvals, curve1,lineStyle,'Color',area_colors{ar});
        plot(subAxBurst,xvals, curve2,lineStyle,'Color',area_colors{ar});
        plot(axAreaComp_abs(i),xvals, curve2,lineStyle,'Color',area_colors{ar});
        
        x = xvals(~isnan(curve1));
        curve1 = curve1(~isnan(curve1));
        curve2 = curve2(~isnan(curve2));
        fill(subAxBurst,[x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
            'FaceColor',area_colors{ar},'EdgeColor','none','FaceAlpha',0.1);
        fill(axAreaComp_abs(i),[x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
            'FaceColor',area_colors{ar},'EdgeColor','none','FaceAlpha',0.1);

        p_TrialCompBurst(i) = plot(subAxBurst,xvals,mean(fr_burst_area{i},1,'omitnan'),lineStyle,'Color',area_colors{ar},'LineWidth',1.5);
        p_AreaComp_abs(i,ar) = plot(axAreaComp_abs(i),xvals,mean(fr_burst_area{i},1,'omitnan'),lineStyle,'Color',area_colors{ar},'LineWidth',1.5);
        
        % Min-max normalization of burst firing rates
        burstMean = mean(fr_burst_area{i},1,'omitnan');
        normVals = (burstMean-min(burstMean))/(max(burstMean)-min(burstMean));
        p_AreaComp_norm(i,ar) = plot(axAreaComp_norm(i),xvals,normVals,lineStyle,'Color',area_colors{ar},'LineWidth',1.5);

        % Min-max normalization of smoothed burst firing rates
        normVals = (movmean(burstMean,smoothWind)-min(movmean(burstMean,smoothWind)))/...
            (max(movmean(burstMean,smoothWind))-min(movmean(burstMean,smoothWind)));
        p_AreaComp_smooth(i,ar) = plot(axAreaComp_smooth(i),xvals,normVals,lineStyle,'Color',area_colors{ar},'LineWidth',1.5);

        title(subAxBurst,sprintf('%s',area_names{ar}))
        ylabel(subAxBurst,'Burst spike firing rate [Hz]')

        % Plot tonic spike rates
        curve1 = mean(fr_tonic_area{i},1,'omitnan') + std(fr_tonic_area{i},[],1,'omitnan')/sqrt(width(fr_tonic_area{i}));
        curve2 = mean(fr_tonic_area{i},1,'omitnan') - std(fr_tonic_area{i},[],1,'omitnan')/sqrt(width(fr_tonic_area{i}));
        
        plot(subAxTonic,xvals, curve1,lineStyle,'Color',area_colors{ar});
        plot(subAxTonic,xvals, curve2,lineStyle,'Color',area_colors{ar});
        
        x = xvals(~isnan(curve1));
        curve1 = curve1(~isnan(curve1));
        curve2 = curve2(~isnan(curve2));
        fill(subAxTonic,[x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
            'FaceColor',area_colors{ar},'EdgeColor','none','FaceAlpha',0.1);

        p_TrialCompTonic(i) = plot(subAxTonic,xvals,mean(fr_tonic_area{i},1,'omitnan'),lineStyle,'Color',area_colors{ar},'LineWidth',1.5);

        title(subAxTonic,sprintf('%s',area_names{ar}))
        ylabel(subAxTonic,'Tonic spike firing rate [Hz]')

        if ar==1 && i==numel(trialType)
            axLegBurst = legend(subAxBurst,[p_TrialCompBurst(:)],trialType,'Location','northwest','AutoUpdate','off');
            axLegTonic = legend(subAxTonic,[p_TrialCompTonic(:)],trialType,'Location','northwest','AutoUpdate','off');
        elseif ar==numel(area_names)
            legend(axAreaComp_abs(i),[p_AreaComp_abs(i,:)],area_names,'Location','northwest','AutoUpdate','off')
            legend(axAreaComp_norm(i),[p_AreaComp_norm(i,:)],area_names,'Location','northwest','AutoUpdate','off')
            legend(axAreaComp_smooth(i),[p_AreaComp_smooth(i,:)],area_names,'Location','northwest','AutoUpdate','off')
        end
        hold(subAxBurst,'off')
        hold(subAxTonic,'off')
        hold(axAreaComp_abs(i),'off')
        hold(axAreaComp_norm(i),'off')
        hold(axAreaComp_smooth(i),'off')

        % Subplot the burst proportions for each trial type and area
        subAxBurstSpikeProp(i) = subplot(2, 2, ar, 'Parent', figBurstSpikesProp(i)); %#ok<SAGROW>
        hold(subAxBurstSpikeProp(i),'on')

        % Plot burst spike proportions
        burstSpikesPropVals = fr_burst_area{i}./fr_all_area{i};
        curve1 = mean(burstSpikesPropVals,1,'omitnan') + std(burstSpikesPropVals,[],1,'omitnan')/sqrt(width(burstSpikesPropVals));
        curve2 = mean(burstSpikesPropVals,1,'omitnan') - std(burstSpikesPropVals,[],1,'omitnan')/sqrt(width(burstSpikesPropVals));
        
        plot(subAxBurstSpikeProp(i),xvals, curve1,lineStyle,'Color',area_colors{ar});
        plot(subAxBurstSpikeProp(i),xvals, curve2,lineStyle,'Color',area_colors{ar});
        
        x = xvals(~isnan(curve1));
        curve1 = curve1(~isnan(curve1));
        curve2 = curve2(~isnan(curve2));
        fill(subAxBurstSpikeProp(i),[x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
            'FaceColor',area_colors{ar},'EdgeColor','none','FaceAlpha',0.1);
        fill(subAxBurstSpikeProp(i),[x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
            'FaceColor',area_colors{ar},'EdgeColor','none','FaceAlpha',0.1);

        plot(subAxBurstSpikeProp(i),xvals,mean(burstSpikesPropVals,1,'omitnan'),'-','Color',area_colors{ar},'LineWidth',1.5);

        title(subAxBurstSpikeProp(i),sprintf('%s',area_names{ar}))
        ylabel(subAxBurstSpikeProp(i),'Burst spike proportion')
        hold(subAxBurstSpikeProp(i),'off')

        % If burst events were sampled at the same time frame, plot burst
        % event proportion
        if binSize == binSizeBursts
            subAxBurstEventProp(i) = subplot(2, 2, ar, 'Parent', figBurstEventProp(i)); %#ok<SAGROW>
            hold(subAxBurstEventProp(i),'on')
            
            burst_area = burstActivityAll{i}(idx,:);
            burstEventPropVals = burst_area./(burst_area + fr_tonic_area{i});

            curve1 = mean(burstEventPropVals,1,'omitnan') + std(burstEventPropVals,[],1,'omitnan')/sqrt(width(burstEventPropVals));
            curve2 = mean(burstEventPropVals,1,'omitnan') - std(burstEventPropVals,[],1,'omitnan')/sqrt(width(burstEventPropVals));

            plot(subAxBurstEventProp(i),xvals, curve1,lineStyle,'Color',area_colors{ar});
            plot(subAxBurstEventProp(i),xvals, curve2,lineStyle,'Color',area_colors{ar});

            x = xvals(~isnan(curve1));
            curve1 = curve1(~isnan(curve1));
            curve2 = curve2(~isnan(curve2));
            fill(subAxBurstEventProp(i),[x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
                'FaceColor',area_colors{ar},'EdgeColor','none','FaceAlpha',0.1);
            fill(subAxBurstEventProp(i),[x fliplr(x)], [curve1 fliplr(curve2)],[0 0 .85],...
                'FaceColor',area_colors{ar},'EdgeColor','none','FaceAlpha',0.1);

            plot(subAxBurstEventProp(i),xvals,mean(burstEventPropVals,1,'omitnan'),'-','Color',area_colors{ar},'LineWidth',1.5);

            title(subAxBurstEventProp(i),sprintf('%s',area_names{ar}))
            ylabel(subAxBurstEventProp(i),'Burst event proportion')
            hold(subAxBurstEventProp(i),'off')
        end
    end
    
    % Compare the different trial types burst firing rates (only works for
    % two trial types)
    if ar==1 
        subAxBurst_wind = copyobj([axLegBurst,subAxBurst],figTrialBurstComp_wind);
        subAxBurst_wind = subAxBurst_wind(2);
    else
        subAxBurst_wind = copyobj(subAxBurst,figTrialBurstComp_wind);
    end
    
    if numel(trialType)==2
        % Compare each time point for significance with bonferroni correction
        [H_indi, p_indi] = ttest2(fr_burst_area{1}, fr_burst_area{2},'Alpha',0.05/width(fr_burst_area{1}),'Vartype','unequal');
        figTrialComp_indi.UserData = struct('statTest', 'ttest2',...
            'Vartype','unequal, i.e. Welsch t-test',...
            'alphaCorrection','bonferroni');
%         [p_indi, H_indi] = arrayfun(@(x) signrank(fr_burst_area{1}(:,x), fr_burst_area{2}(:,x),'Alpha',0.05/width(fr_burst_area{1})), 1:width(fr_burst_area{1}));
        % False discovery rate correction after Hochberg
%         [p_indi, H_indi] = mafdr(p_indi);

        % Convert subsequent significant data points into continuous bars
        barBegins = strfind(H_indi, [0 1]);
        if H_indi(1)
            barBegins = [1, barBegins+1]; 
        end
        barEnds = strfind(H_indi, [1 0]);
        if H_indi(end)
            barEnds = [barEnds, numel(H_indi)]; %#ok<AGROW> 
        end
        % Plot bars underneath the respective plots
        if ~isempty(barBegins)
            % Plot the significance bar 20% lower than the line plots
            lines = findobj(subAxBurst, 'Type', 'line');
            yvals = [lines(:).YData];
            barHeight = min(yvals)-(max(yvals)-min(yvals))*0.2;
            
            hold(subAxBurst,'on')
            for i = 1:numel(barBegins)
                plot(subAxBurst,[xvals(barBegins(i)),xvals(barEnds(i))],...
                   [barHeight, barHeight], 'Color',area_colors{ar},'LineWidth',6)
            end
            % 10% margin at top and bottom
            ylim(subAxBurst,[barHeight-(max(yvals)-min(yvals))*0.1 max(yvals)+(max(yvals)-min(yvals))*0.1])
            
            hold(subAxBurst,'off')
        end

        % Compare time windows for significance with bonferroni correction
        xvalsWind = milliseconds(timeLapse(1):windowLength:timeLapse(2)-windowLength);
        
        % For plotting burst rate statistical results
        H_wind = NaN(1,windowBins);
        p_wind = NaN(1,windowBins);
        for wind = 1:windowBins
            fact = width(fr_burst_area{1})/windowBins;
            startPoint = 1+(wind-1)*fact;
%             [h_temp, p_temp] = ttest(mean(fr_burst_area{1}(:,startPoint:startPoint+fact-1),'omitnan'),...
%                 mean(fr_burst_area{2}(:,startPoint:startPoint+fact-1),'omitnan'),'Alpha',0.05/windowBins);
            % Account for difference in variability
            [h_temp, p_temp] = ttest2(mean(fr_burst_area{1}(:,startPoint:startPoint+fact-1),2,'omitnan'),...
                mean(fr_burst_area{2}(:,startPoint:startPoint+fact-1),2,'omitnan'),'Alpha',0.05/(localBonferroni/windowLength),'Vartype','unequal');
            H_wind(wind) = h_temp;
            p_wind(wind) = p_temp;
        end
        figTrialBurstComp_wind.UserData = struct('statTest', 'ttest2',...
            'Vartype','unequal, i.e. Welsch t-test',...
            'alphaCorrection',bonferroniString);
        % Convert subsequent significant data points into continuous bars
        barBegins = strfind(H_wind, [0 1]);
        if H_wind(1)
            barBegins = [1, barBegins+1]; 
        end
        barEnds = strfind(H_wind, [1 0]);
        if H_wind(end)
            barEnds = [barEnds, numel(H_wind)]; %#ok<AGROW> 
        end
        % Plot bars underneath the respective plots
        if ~isempty(barBegins)
            % Plot the significance bar 20% lower than the line plots
            lines = findobj(subAxBurst_wind, 'Type', 'line');
            yvals = [lines(:).YData];
            barHeight = min(yvals)-(max(yvals)-min(yvals))*0.2;
            
            hold(subAxBurst_wind,'on')
            for i = 1:numel(barBegins)
                plot(subAxBurst_wind,[xvalsWind(barBegins(i)),xvalsWind(barEnds(i))+milliseconds(windowLength)],...
                   [barHeight, barHeight], 'Color',area_colors{ar},'LineWidth',6)
            end
            % 10% margin at top and bottom
            ylim(subAxBurst_wind,[barHeight-(max(yvals)-min(yvals))*0.1 max(yvals)+(max(yvals)-min(yvals))*0.1])
            
            hold(subAxBurst_wind,'off')
        end

        % For plotting tonic rate statistical results
        H_wind = NaN(1,windowBins);
        p_wind = NaN(1,windowBins);
        for wind = 1:windowBins
            fact = width(fr_tonic_area{1})/windowBins;
            startPoint = 1+(wind-1)*fact;
%             [h_temp, p_temp] = ttest(mean(fr_tonic_area{1}(:,startPoint:startPoint+fact-1),'omitnan'),...
%                 mean(fr_tonic_area{2}(:,startPoint:startPoint+fact-1),'omitnan'),'Alpha',0.05/windowBins);
            % Account for difference in variability
            [h_temp, p_temp] = ttest2(mean(fr_tonic_area{1}(:,startPoint:startPoint+fact-1),2,'omitnan'),...
                mean(fr_tonic_area{2}(:,startPoint:startPoint+fact-1),2,'omitnan'),'Alpha',0.05/(localBonferroni/windowLength),'Vartype','unequal');
            H_wind(wind) = h_temp;
            p_wind(wind) = p_temp;
        end
        figTrialTonicComp_wind.UserData = struct('statTest', 'ttest2',...
            'Vartype','unequal, i.e. Welsch t-test',...
            'alphaCorrection',bonferroniString);
        % Convert subsequent significant data points into continuous bars
        barBegins = strfind(H_wind, [0 1]);
        if H_wind(1)
            barBegins = [1, barBegins+1]; 
        end
        barEnds = strfind(H_wind, [1 0]);
        if H_wind(end)
            barEnds = [barEnds, numel(H_wind)]; %#ok<AGROW> 
        end
        % Plot bars underneath the respective plots
        if ~isempty(barBegins)
            % Plot the significance bar 20% lower than the line plots
            lines = findobj(subAxTonic, 'Type', 'line');
            yvals = [lines(:).YData];
            barHeight = min(yvals)-(max(yvals)-min(yvals))*0.2;
            
            hold(subAxTonic,'on')
            for i = 1:numel(barBegins)
                plot(subAxTonic,[xvalsWind(barBegins(i)),xvalsWind(barEnds(i))+milliseconds(windowLength)],...
                   [barHeight, barHeight], 'Color',area_colors{ar},'LineWidth',6)
            end
            % 10% margin at top and bottom
            ylim(subAxTonic,[barHeight-(max(yvals)-min(yvals))*0.1 max(yvals)+(max(yvals)-min(yvals))*0.1])
            
            hold(subAxTonic,'off')
        end

        % For plotting overall mean firing rate statistical results
        H_wind = NaN(1,windowBins);
        p_wind = NaN(1,windowBins);
        for wind = 1:windowBins
            fact = width(fr_all_area{1})/windowBins;
            startPoint = 1+(wind-1)*fact;
%             [h_temp, p_temp] = ttest(mean(fr_all_area{1}(:,startPoint:startPoint+fact-1),'omitnan'),...
%                 mean(fr_all_area{2}(:,startPoint:startPoint+fact-1),'omitnan'),'Alpha',0.05/windowBins);
            % Account for difference in variability
            [h_temp, p_temp] = ttest2(mean(fr_all_area{1}(:,startPoint:startPoint+fact-1),2,'omitnan'),...
                mean(fr_all_area{2}(:,startPoint:startPoint+fact-1),2,'omitnan'),'Alpha',0.05/(localBonferroni/windowLength),'Vartype','unequal');
            H_wind(wind) = h_temp;
            p_wind(wind) = p_temp;
        end
        figTrialAllComp_wind.UserData = struct('statTest', 'ttest2',...
            'Vartype','unequal, i.e. Welsch t-test',...
            'alphaCorrection',bonferroniString);
        % Convert subsequent significant data points into continuous bars
        barBegins = strfind(H_wind, [0 1]);
        if H_wind(1)
            barBegins = [1, barBegins+1]; 
        end
        barEnds = strfind(H_wind, [1 0]);
        if H_wind(end)
            barEnds = [barEnds, numel(H_wind)]; %#ok<AGROW> 
        end
        % Plot bars underneath the respective plots
        if ~isempty(barBegins)
            % Plot the significance bar 20% lower than the line plots
            lines = findobj(subAxAll, 'Type', 'line');
            yvals = [lines(:).YData];
            barHeight = min(yvals)-(max(yvals)-min(yvals))*0.2;
            
            hold(subAxAll,'on')
            for i = 1:numel(barBegins)
                plot(subAxAll,[xvalsWind(barBegins(i)),xvalsWind(barEnds(i))+milliseconds(windowLength)],...
                   [barHeight, barHeight], 'Color',area_colors{ar},'LineWidth',6)
            end
            % 10% margin at top and bottom
            ylim(subAxAll,[barHeight-(max(yvals)-min(yvals))*0.1 max(yvals)+(max(yvals)-min(yvals))*0.1])
            
            hold(subAxAll,'off')
        end
    end
end

%% Save figures
% Destination folder for matlab .fig files
destfold = fullfile(cohortPath,'Analysis-Figures','Population-Burst-PSTH',sessionDescription{:});
if exist(destfold,"dir") == 0
    mkdir(destfold)
end

figHandles = handle(sort(double(findall(0, 'type', 'figure'))));

fprintf("\nSaving Figures.\n")

for i = 1:numel(figHandles)
    % Only save the population figures, not the individual units
    if contains(figHandles(i).Name,'BurstEventRate')
        destfile = sprintf('%s\\%s_%s.fig', destfold, figHandles(i).Name,...
            sprintf('%s_%imsecBinSizeFR_%imsecBinSizeBursts_timeLapse%i-%i_%s',strjoin(trialType,'&'), ...
            binSize,binSizeBursts,timeLapse(1),timeLapse(2),cellDescript));
        savefig(figHandles(i), destfile);
    elseif contains(figHandles(i).Name,{'AreaComp','BurstSpikesProp'})
        destfile = sprintf('%s\\%s_%s.fig', destfold, figHandles(i).Name,...
            sprintf('%imsecBinSizeFR_timeLapse%i-%i_%s',binSize,timeLapse(1),timeLapse(2),cellDescript));
        savefig(figHandles(i), destfile);
    elseif contains(figHandles(i).Name,{'BurstEventProp'})
        destfile = sprintf('%s\\%s_%s.fig', destfold, figHandles(i).Name,...
            sprintf('%imsecBinSizeFR_%imsecBinSizeBursts_timeLapse%i-%i_%s', ...
            binSize,binSizeBursts,timeLapse(1),timeLapse(2),cellDescript));
        savefig(figHandles(i), destfile);
    elseif contains(figHandles(i).Name,{'BurstTonicFR','BurstRateComp','TonicRateComp','FiringRateComp'})
        destfile = sprintf('%s\\%s_%s.fig', destfold, figHandles(i).Name,...
            sprintf('%s_%imsecBinSizeFR_timeLapse%i-%i_%s',strjoin(trialType,'&'), ...
            binSize,timeLapse(1),timeLapse(2),cellDescript));
        savefig(figHandles(i), destfile);
    end
end

fprintf("\nDone!\n\n")

%% Helper functions

function doneExeWave(~,~,f,ui_field,area_names)

area_waveforms = cell(1,numel(area_names));
for ar = 1:numel(area_names)
    switch ui_field(ar).Value
        case 2
            area_waveforms{ar} = 'RS';
        case 3
            area_waveforms{ar} = 'FS';
        otherwise
            area_waveforms{ar} = 'All';
    end
end

% Overwrite variables
assignin('caller','area_waveforms',area_waveforms)

close(f)
end

function [put_in_idx, put_ex_idx] = getWaveformDistribution(clWaveforms,cellArea)
trough2peak_dur = nan(1,size(clWaveforms,1));
for unit = 1:size(clWaveforms,1)
    mean_waveform = mean(clWaveforms{unit,2},2);
    [~,troughs_loc,~,~] = findpeaks(-mean_waveform,'SortStr','descend','NPeaks',1);
    [~,peak_loc,~,~] = findpeaks(mean_waveform);
    if ~isempty(find(peak_loc > troughs_loc,1))
        peak_loc = peak_loc(find(peak_loc > troughs_loc,1));
        trough2peak_dur(unit) = abs(diff([peak_loc,troughs_loc]));
    else
        try
            peak_loc = peak_loc(find(peak_loc < troughs_loc,1,'last'));
            trough2peak_dur(unit) = abs(diff([peak_loc,troughs_loc]));
        catch
            trough2peak_dur(unit) = nan;
        end
    end
end

% Widths are given in data points. Calculate in usecs
framerate = 30000; % in Hz
trough2peak_usec = (trough2peak_dur/framerate)*1000000;

put_in_idx = false(height(clWaveforms),1);
put_ex_idx = false(height(clWaveforms),1);

% Choose trough to peak cutoff of 300 µs for thalamic nuclei
idx = find(cellfun(@(x) ismember(x, {'VPM','POm'}),cellArea));
put_in_idx(idx) = trough2peak_usec(idx) < 300;
put_ex_idx(idx) = trough2peak_usec(idx) >= 300;

% Choose trough to peak cutoff of 350 µs for cortex and ZI
idx = find(cellfun(@(x) ismember(x, {'BC','ZIv'}),cellArea));
put_in_idx(idx) = trough2peak_usec(idx) < 350;
put_ex_idx(idx) = trough2peak_usec(idx) >= 350;

end
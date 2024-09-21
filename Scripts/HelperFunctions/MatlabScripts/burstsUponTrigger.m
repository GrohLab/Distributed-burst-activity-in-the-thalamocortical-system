%% Bursts before and after trigger onset
% This script calculates the change in burstiness for a given trigger and
% the respective recorded brain areas.
% The "getBurstiness.m" script already calculates important information
% about burst on- or offsets. Load the respective "BurstinessData.mat"
% files and the triggers, in order to calculate the trigger response.

% Variables to set: fileSelection, chCond, adjustConditions, spontWindow,
% respWindow, trigOffset

close all; clearvars; clc

% Change this to track changes and automatically re-run the script, if files are not up to date
fileVersion = 2; 

% Access correct individual
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

filesMissing = cellfun(@isfile, fullfile(fileSelection, 'BurstinessData.mat'));
assert(all(filesMissing),'burstsUponTrigger:notFullyAnalyzed',...
    ['Not all files have been burst analyzed ("getBurstiness.m"). BurstinessData.mat missing in:',strrep(sprintf('\n- %s', fileSelection{~filesMissing}),'\','\\')])

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
% Strip the string from any special characters in the beginning or end
sessionDescription = regexprep(sessionDescription,'[^a-zA-Z0-9]*$','');
sessionDescription = regexprep(sessionDescription,'^[^a-zA-Z0-9]*','');
sessionDescription = strrep(sessionDescription,' ','');

% Define brain areas
area_names = {'BC','VPM','POm','ZIv'};

% Choose condition to analyze
prompt = {'Reward', 'Punishment', 'Lick', 'onlyFirstLick',...
    'WhiskerContact_left', 'WhiskerContact_right','Middlepoint',...
    'Random', 'WhiskerContact_onlyLeftFirst', 'WhiskerContact_onlyRightFirst'};
condPick = listdlg('PromptString','Pick condition to analyze.', ...
    'ListString',prompt,...
    'ListSize',[300 200],'InitialValue',5,'SelectionMode','single');
condition = prompt{condPick};

% Choose time windows to compare
spontWindow = [-1000,-800]; % spontaneous window in msec (default: [-600,-400])
respWindow = [0,200]; % response window in msec (default: [0,-200])
trigOffset = false; % for response on 'offset', set to true
targetEstimate = true;

assert(round(diff(respWindow),3)==round(diff(spontWindow),3),...
    'Salience_function:WindowsNotSameLength',...
    'responseWindow and spontaneousWindow must be of the same length.')
windowLength = round(diff(spontWindow),3);

%% Loop through files to create ClusterBurstiness variable
% For every occuring event retrieve the number of bursts and the total
% number of spikes.
close all;

% For only looking at specific trial types. Variable must be a cell array,
% with only one or two of the following values:
% 'wide', 'narrow', 'intermediate', 'lick', 'no-lick'
% To analyze all trial types, set variable to 'all'
prompt = {'all', 'wide', 'narrow', 'intermediate', 'lick', 'no-lick'};
trialPick = listdlg('PromptString',{'Pick trial type to analyze.'; 'An aperture setting can be combined with an outcome.'}, ...
    'ListString',prompt,...
    'ListSize',[300 200],'InitialValue',2,'SelectionMode','multiple');
if numel(trialPick) > 2
    error('burstUponTrigger:trialPickError','No more than 2 trial types can be picked.')
elseif numel(trialPick) > 1
    assert(ismember(trialPick(1),(2:4)) && ismember(trialPick(2),(5:6)),...
        'burstUponTrigger:trialPickError',...
        'When picking two trial types, one must be an aperture setting, and the other an outcome type.')
end
trialType = prompt(trialPick);

% Cutoff for putative inhibitory cells in µsec
prompt = {'all', 'putative inhibitory', 'putative excitatory'};
waveformPick = listdlg('PromptString','Which units do you want to analyze?', ...
    'ListString',prompt,...
    'ListSize',[200 150]);
waveformType = prompt(waveformPick);

if targetEstimate
    outputFileName = ['ClusterBurstiness_', condition, '_', cell2mat(join(trialType,'&')),...
        '_withTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f.mat')];
    responseFileName = ['ResponsePattern_', condition, '_', cell2mat(join(trialType,'&')),...
        '_withTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f.mat')];
else
    outputFileName = ['ClusterBurstiness_', condition, '_', cell2mat(join(trialType,'&')),...
        '_withoutTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f.mat')]; %#ok<UNRCH>
    responseFileName = ['ResponsePattern_', condition, '_', cell2mat(join(trialType,'&')),...
        '_withoutTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f.mat')];
end

for ses = 1:numel(fileSelection)

    % Set warning to error in order to catch it
    s = warning('error', 'MATLAB:load:variableNotFound');

    try load(fullfile(fileSelection{ses},outputFileName),'update','ClusterBurstiness','trialNum')
        % Version of file (aka 'update') must be the most recent
        assert(update >= 3)
    catch
        clear update
        load(fullfile(fileSelection{ses},'BurstinessData.mat'),'BurstStarts','TypeOfSpike','NumSpikesInBursts','BurstFrequencies','unitIDs')

        dataInfo = dir(fullfile(fileSelection{ses},'*all_channels.mat'));
        try
            load(fullfile(dataInfo(1).folder,dataInfo(1).name),'sortedData','fs')
            sortedData = sortedData(all(~cellfun(@isempty, sortedData),2),:);
            clInfo = getClusterInfo(fullfile(fileSelection{ses},'cluster_info.tsv'));
        catch
            fprintf('\nCould not gather sortedData and clInfo.\n')
            return
        end

        str_idx = regexp(fileSelection{ses},'#\d*','end');
        animalPath = fileSelection{ses}(1:str_idx);
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

        tempDir = fileSelection{ses};
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

        % Get the Triggers
        % Set warning to error in order to catch it
        s = warning('error', 'MATLAB:load:variableNotFound');
        condInfo = dir(fullfile(fileparts(fileSelection{ses}),'*analysis.mat'));
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

        % trial types must have been analyzed for the respective stimulus
        % condition
        FileInfo = dir(fullfile(fileSelection{ses},sprintf('StimulusResponse_%s*',condition)));
        if ~any(cellfun(@isempty, trialType)) && ~isempty(FileInfo) && ~any(contains(trialType,'all'))
            load(fullfile(FileInfo(1).folder,FileInfo(1).name),'Wide_Narrow_Intermediate','Lick')
            if any(ismember(trialType,'wide'))
                logApert = Wide_Narrow_Intermediate==1;
            elseif any(ismember(trialType,'narrow'))
                logApert = Wide_Narrow_Intermediate==2;
            elseif any(ismember(trialType,'intermediate'))
                logApert = Wide_Narrow_Intermediate==3;
            else
                logApert = true(length(Wide_Narrow_Intermediate),1);
            end

            if any(ismember(trialType,'lick'))
                logLick = Lick==1;
            elseif any(ismember(trialType,'no-lick'))
                logLick = Lick==0;
            else
                logLick = true(length(Lick),1);
            end
            logicalFlag = logApert & logLick;
        end

        % Adjust Triggers accordingly
        idx = find(cellfun(@(x) isequal(x,condition), {Conditions.name}));
        if ~any(cellfun(@isempty, trialType)) && ~any(contains(trialType,'all'))
            for i = 1:numel(idx)
                Conditions(idx(i)).Triggers = sort(Conditions(idx(i)).Triggers);
                Conditions(idx(i)).Triggers = Conditions(idx(i)).Triggers(logicalFlag,:);
            end
        end

        if trigOffset
            triggers = sort(Conditions(idx(1)).Triggers(:,2)); %#ok<UNRCH>
        else
            triggers = sort(Conditions(idx(1)).Triggers(:,1));
        end
        % Convert trigger points from samples to milliseconds
        triggers = 1000.*triggers./fs;

        if ~isequal(sortedData(:,1),unitIDs)
            error('burstsUponTrigger:unitMismatch',...
                'Units IDs of analyzed BurstinessData.mat file and sortedData do not match.')
        end

        % Get burst metrics for every unit seperately
        ClusterBurstiness = table('Size',[height(sortedData), 16],'VariableTypes',...
            [repmat({'cell'},1,4), repmat({'double'},1,12)],'VariableNames',...
            {'unit','area','spikeTypePerTrialPre','spikeTypePerTrialPost',...
            'spikesPre','burstsPre','burstSpikesPre','nonBurstSpikesPre',...
            'avgSpikesPerBurstPre','avgBurstFrequencyPre',...
            'spikesPost','burstsPost','burstSpikesPost','nonBurstSpikesPost',...
            'avgSpikesPerBurstPost','avgBurstFrequencyPost'});
        ClusterBurstiness.unit = sortedData(:,1);
        ClusterBurstiness.area = sortedData(:,4);
        for unit = 1:height(sortedData)
            spiketimes_msec = sortedData{unit,2}.*1000;
            ClusterBurstiness.spikeTypePerTrialPre{unit} = cell(numel(triggers),1);
            ClusterBurstiness.spikeTypePerTrialPost{unit} = cell(numel(triggers),1);
            for trig = 1:numel(triggers)
                % Spikes in spontaneous window
                log_idx = (spiketimes_msec - triggers(trig)) >= spontWindow(1) & (spiketimes_msec - triggers(trig)) < spontWindow(2);
                burst_idx = (BurstStarts{unit} - triggers(trig)) >= spontWindow(1) & (BurstStarts{unit} - triggers(trig)) < spontWindow(2);
                spike_type = TypeOfSpike{1,unit}(log_idx);

                ClusterBurstiness.spikeTypePerTrialPre{unit}{trig} = spike_type;
                ClusterBurstiness.spikesPre(unit) = ClusterBurstiness.spikesPre(unit)+numel(spike_type);
                ClusterBurstiness.nonBurstSpikesPre(unit) = ClusterBurstiness.nonBurstSpikesPre(unit)+sum(~spike_type);
                ClusterBurstiness.burstSpikesPre(unit) = ClusterBurstiness.burstSpikesPre(unit)+sum(spike_type);
                ClusterBurstiness.burstsPre(unit) = ClusterBurstiness.burstsPre(unit)+sum(burst_idx);
                if any(burst_idx)
                    ClusterBurstiness.avgSpikesPerBurstPre(unit) = ClusterBurstiness.avgSpikesPerBurstPre(unit) + sum(NumSpikesInBursts{unit}(burst_idx));
                    ClusterBurstiness.avgBurstFrequencyPre(unit) = ClusterBurstiness.avgBurstFrequencyPre(unit) + sum(BurstFrequencies{unit}(burst_idx));
                end

                % Spikes in response window
                log_idx = (spiketimes_msec - triggers(trig)) >= respWindow(1) & (spiketimes_msec - triggers(trig)) < respWindow(2);
                burst_idx = (BurstStarts{unit} - triggers(trig)) >= respWindow(1) & (BurstStarts{unit} - triggers(trig)) < respWindow(2);
                spike_type = TypeOfSpike{1,unit}(log_idx);

                ClusterBurstiness.spikeTypePerTrialPost{unit}{trig} = spike_type;
                ClusterBurstiness.spikesPost(unit) = ClusterBurstiness.spikesPost(unit)+numel(spike_type);
                ClusterBurstiness.nonBurstSpikesPost(unit) = ClusterBurstiness.nonBurstSpikesPost(unit)+sum(~spike_type);
                ClusterBurstiness.burstSpikesPost(unit) = ClusterBurstiness.burstSpikesPost(unit)+sum(spike_type);
                ClusterBurstiness.burstsPost(unit) = ClusterBurstiness.burstsPost(unit)+sum(burst_idx);
                if any(burst_idx)
                    ClusterBurstiness.avgSpikesPerBurstPost(unit) = ClusterBurstiness.avgSpikesPerBurstPost(unit) + sum(NumSpikesInBursts{unit}(burst_idx));
                    ClusterBurstiness.avgBurstFrequencyPost(unit) = ClusterBurstiness.avgBurstFrequencyPost(unit) + sum(BurstFrequencies{unit}(burst_idx));
                end
            end
            % Devide by number of bursts
            ClusterBurstiness.avgSpikesPerBurstPre(unit) = ClusterBurstiness.avgSpikesPerBurstPre(unit)/ClusterBurstiness.burstsPre(unit);
            ClusterBurstiness.avgBurstFrequencyPre(unit) = ClusterBurstiness.avgBurstFrequencyPre(unit)/ClusterBurstiness.burstsPre(unit);
            ClusterBurstiness.avgSpikesPerBurstPost(unit) = ClusterBurstiness.avgSpikesPerBurstPost(unit)/ClusterBurstiness.burstsPost(unit);
            ClusterBurstiness.avgBurstFrequencyPost(unit) = ClusterBurstiness.avgBurstFrequencyPost(unit)/ClusterBurstiness.burstsPost(unit);
        end

        trialNum = numel(triggers);
        % Current version of the script in order to keep track of changes
        update = 3;
        save(fullfile(fileSelection{ses},outputFileName),'update','ClusterBurstiness','trialNum')
    end
    % Reset warning
    warning(s);

end

%% Concatenate information of selected files

burstRatio_pre = cell(1,numel(fileSelection));
burstRatio_post = cell(1,numel(fileSelection));

AP_belonging_to_bursts_pre = cell(1,numel(fileSelection));
AP_belonging_to_bursts_post = cell(1,numel(fileSelection));
burstAPs_signif = cell(1,numel(fileSelection));
tonicAPs_signif = cell(1,numel(fileSelection));

non_burst_spikesPerTrial_pre = cell(1,numel(fileSelection));
non_burst_spikesPerTrial_post = cell(1,numel(fileSelection));

burstsPerTrial_pre = cell(1,numel(fileSelection));
burstsPerTrial_post = cell(1,numel(fileSelection));

spikesPerBurst_pre = cell(1,numel(fileSelection));
spikesPerBurst_post = cell(1,numel(fileSelection));

burstFrequency_pre = cell(1,numel(fileSelection));
burstFrequency_post = cell(1,numel(fileSelection));

putEx_allFiles = cell(1,numel(fileSelection));
trialNum_allFiles = cell(1,numel(fileSelection));
areas_allFiles = cell(1,numel(fileSelection));
units_allFiles = cell(1,numel(fileSelection));
sessions_allFiles = cell(1,numel(fileSelection));

% singleCellFileName = ['SingleCellReactivity_', condition, '_withTargetEstimate_', num2str(respWindow/1000,'respond%.3f-%.3f_'), num2str(spontWindow/1000,'baseline%.3f-%.3f.mat')];

for ses = 1:numel(fileSelection)
    load(fullfile(fileSelection{ses},outputFileName),'ClusterBurstiness','trialNum')
    %     load(fullfile(fileSelection{ses},singleCellFileName),'SingleCellResponse')

    assert(~isempty(dir(fullfile(fileSelection{ses}, '*waveforms_all.mat'))), ...
        'burstUponTrigger:waveformFileEror',sprintf('No waveforms_all.mat file for dir "%s"',strrep(fileSelection{ses},'\','\\')))

    WaveformInfo = dir(fullfile(fileSelection{ses}, '*waveforms_all.mat'));
    load(fullfile(WaveformInfo.folder,WaveformInfo.name),'clWaveforms')
    if all(cellfun(@(x,y) isequal(x,y), ClusterBurstiness.unit,clWaveforms(:,1)))
        [put_in_idx, put_ex_idx] = getWaveformDistribution(clWaveforms,ClusterBurstiness.area);
    elseif all(cellfun(@(x) any(ismember(ClusterBurstiness.unit,x)), clWaveforms(:,1)))
        % Order the units according to the ClusterBurstiness table
        unitOrder = nan(height(clWaveforms),1);
        for i = 1:height(ClusterBurstiness)
            unitOrder(i) = find(ismember(clWaveforms(:,1),ClusterBurstiness.unit{i}));
        end
        clWaveforms = clWaveforms(unitOrder,:);

        [put_in_idx, put_ex_idx] = getWaveformDistribution(clWaveforms,ClusterBurstiness.area);
    else
        error('burstUponTrigger:waveformUnitsError','Units in the clWaveforms variable do not match the ones from the ClusterBurstiness table.')
    end
    putEx_allFiles{ses} = put_ex_idx;

    % Filter ClusterBurstiness table for waveform type
    if ismember(waveformType,'putative inhibitory')
        ClusterBurstiness = ClusterBurstiness(put_in_idx,:);
        putEx_allFiles{ses} = putEx_allFiles{ses}(putEx_allFiles{ses}==0);
    elseif ismember(waveformType, 'putative excitatory')
        ClusterBurstiness = ClusterBurstiness(put_ex_idx,:);
        putEx_allFiles{ses} = putEx_allFiles{ses}(putEx_allFiles{ses}==1);
    end

    burstRatio_pre{ses} = ClusterBurstiness.burstsPre./(ClusterBurstiness.burstsPre+ClusterBurstiness.nonBurstSpikesPre);
    burstRatio_post{ses} = ClusterBurstiness.burstsPost./(ClusterBurstiness.burstsPost+ClusterBurstiness.nonBurstSpikesPost);

    AP_belonging_to_bursts_pre{ses} = ClusterBurstiness.burstSpikesPre./(ClusterBurstiness.burstSpikesPre+ClusterBurstiness.nonBurstSpikesPre);
    AP_belonging_to_bursts_post{ses} = ClusterBurstiness.burstSpikesPost./(ClusterBurstiness.burstSpikesPost+ClusterBurstiness.nonBurstSpikesPost);

    % Check for each unit whether the proportion of burst spikes changes significantly
    spikeProportions = cellfun(@(c,d) cellfun(@(x,y) [sum(x)/numel(x),sum(y)/numel(y)], c,d,'UniformOutput',false),...
        ClusterBurstiness.spikeTypePerTrialPre,ClusterBurstiness.spikeTypePerTrialPost, 'UniformOutput',false);
    for unit = 1:height(spikeProportions)
        % Substitute NaN values with zeros for statistical testing
        spikeProportions{unit} = cell2mat(spikeProportions{unit});
        spikeProportions{unit}(isnan(spikeProportions{unit})) = 0;
    end
    burstAPs_signif{ses} = spikeProportions;

    % Check for each unit whether the proportion of burst spikes changes significantly
    spikeProportions = cellfun(@(c,d) cellfun(@(x,y) [sum(~x),sum(~y)], c,d,'UniformOutput',false),...
        ClusterBurstiness.spikeTypePerTrialPre,ClusterBurstiness.spikeTypePerTrialPost, 'UniformOutput',false);
    for unit = 1:height(spikeProportions)
        % Substitute NaN values with zeros for statistical testing
        spikeProportions{unit} = cell2mat(spikeProportions{unit});
        spikeProportions{unit}(isnan(spikeProportions{unit})) = 0;
    end
    tonicAPs_signif{ses} = spikeProportions;

    non_burst_spikesPerTrial_pre{ses} = ClusterBurstiness.nonBurstSpikesPre./trialNum;
    non_burst_spikesPerTrial_post{ses} = ClusterBurstiness.nonBurstSpikesPost./trialNum;

    burstsPerTrial_pre{ses} = ClusterBurstiness.burstsPre./trialNum;
    burstsPerTrial_post{ses} = ClusterBurstiness.burstsPost./trialNum;

    spikesPerBurst_pre{ses} = ClusterBurstiness.avgSpikesPerBurstPre;
    spikesPerBurst_post{ses} = ClusterBurstiness.avgSpikesPerBurstPost;

    burstFrequency_pre{ses} = ClusterBurstiness.avgBurstFrequencyPre;
    burstFrequency_post{ses} = ClusterBurstiness.avgBurstFrequencyPost;

    trialNum_allFiles{ses} = trialNum;
    areas_allFiles{ses} = ClusterBurstiness.area;
    units_allFiles{ses} = ClusterBurstiness.unit;
    sessions_allFiles{ses} = repmat(fileSelection(ses),height(ClusterBurstiness),1);
end

% Remove empty cells
burstAPs_signif = vertcat(burstAPs_signif{:});
emptyIdx = cellfun(@isempty, burstAPs_signif);
burstAPs_signif(emptyIdx) = {[NaN, NaN]};

tonicAPs_signif = vertcat(tonicAPs_signif{:});
emptyIdx = cellfun(@isempty, tonicAPs_signif);
tonicAPs_signif(emptyIdx) = {[NaN, NaN]};

burstAPs_signif_all = NaN(numel(burstAPs_signif),1);
burstAPs_signif_all(~emptyIdx) = cellfun(@(x) signrank(x(:,1),x(:,2)), burstAPs_signif(~emptyIdx))<=0.05;

tonicAPs_signif_all = NaN(numel(tonicAPs_signif),1);
tonicAPs_signif_all(~emptyIdx) = cellfun(@(x) signrank(x(:,1),x(:,2)), tonicAPs_signif(~emptyIdx))<=0.05;

areas_all = vertcat(areas_allFiles{:});
sessions_all = vertcat(sessions_allFiles{:});
units_all = vertcat(units_allFiles{:});

putEx_allFiles = horzcat(putEx_allFiles{:})';

burstAPs_Info = table(sessions_all, areas_all, units_all,...
    putEx_allFiles, burstAPs_signif, burstAPs_signif_all,...
    'VariableNames',{'SessionName','Area','UnitID','PutExcitatory','CellResponses','Significance'});
tonicAPs_Info = table(sessions_all, areas_all, units_all,...
    putEx_allFiles, tonicAPs_signif, tonicAPs_signif_all,...
    'VariableNames',{'SessionName','Area','UnitID','PutExcitatory','CellResponses','Significance'});

avgCellResponse = table(sessions_all, areas_all, units_all,putEx_allFiles,...
    table('Size',[numel(areas_all),2],'VariableTypes',{'double','double'},'VariableNames',{'pre','post'}),...
    table('Size',[numel(areas_all),2],'VariableTypes',{'double','double'},'VariableNames',{'pre','post'}),...
    table('Size',[numel(areas_all),2],'VariableTypes',{'double','double'},'VariableNames',{'pre','post'}),...
    table('Size',[numel(areas_all),2],'VariableTypes',{'double','double'},'VariableNames',{'pre','post'}),...
    table('Size',[numel(areas_all),2],'VariableTypes',{'double','double'},'VariableNames',{'pre','post'}),...
    table('Size',[numel(areas_all),2],'VariableTypes',{'double','double'},'VariableNames',{'pre','post'}),...
    'VariableNames',{'SessionName','Area','UnitID','PutExcitatory',...
    'TonicSpikes','SpikesBelongingToBurst','BurstsToTonicRatio','BurstNum','avgSpikesPerBurst','avgBurstFrequency'});

%% Plot results

% Plot non-burst spikes per trial spontaneous vs. evoked
pre = vertcat(non_burst_spikesPerTrial_pre{:});
post = vertcat(non_burst_spikesPerTrial_post{:});

avgCellResponse.TonicSpikes.pre = pre;
avgCellResponse.TonicSpikes.post = post;

for ar = 1:numel(area_names)
    area_idx = contains(areas_all,area_names{ar});
    figure('Name',sprintf('NonBurstSpikesPerTrial_%s',area_names{ar}))
    hold on
    boxplot([pre(area_idx),post(area_idx)], 'Symbol', 'k.','Notch','on')
    parallelcoords([pre(area_idx & tonicAPs_signif_all==false),post(area_idx & tonicAPs_signif_all==false)], 'Color', 0.7*[1 1 1], 'LineStyle', '-');
    parallelcoords([pre(area_idx & tonicAPs_signif_all==true),post(area_idx & tonicAPs_signif_all==true)], 'Color', 0.2*[1 1 1], 'LineStyle', '-');
    axis auto

    % Paired sample t-test
    % consider using signrank for non-parametric testing
    p = signrank(pre(area_idx),post(area_idx));

    yt = get(gca, 'YTick');
    xt = get(gca, 'XTick');
    if p < 0.001
        plot([xt(1)+0.2,xt(2)-0.2], [1.1 1.1]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*1.125, '***','HorizontalAlignment','center')
    elseif p < 0.01
        plot([xt(1)+0.2,xt(2)-0.2], [1.1 1.1]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*1.125, '**','HorizontalAlignment','center')
    elseif p < 0.05
        plot([xt(1)+0.2,xt(2)-0.2], [1.1 1.1]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*1.125, '*','HorizontalAlignment','center')
    else
        plot([xt(1)+0.2,xt(2)-0.2], [1.1 1.1]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*1.125, 'n.s.','HorizontalAlignment','center')
    end

    title(sprintf('Amount of non-burst spikes per trial for given clusters from %s',area_names{ar}))
    ylabel('spikes')
    xticklabels({'spontaneous','evoked'})

    % Add a small pie chart, indicating the proportion of significantly
    % responding units
    axes('Position',[.75 .8 .1 .1]), box off
    piePlot = pie([sum(area_idx & tonicAPs_signif_all==true) sum(area_idx & tonicAPs_signif_all==false)]);
    patchHand = findobj(piePlot, 'Type', 'Patch');
    patchHand(2).FaceColor = 0.7*[1 1 1];
    patchHand(1).FaceColor = 0.2*[1 1 1];
    delete(findobj(piePlot,'Type','text'))

    % Save information on units and sessions for later retrieval
    Info = struct('data', [pre(area_idx & tonicAPs_signif_all==false),post(area_idx & tonicAPs_signif_all==false);...
        pre(area_idx & tonicAPs_signif_all==true),post(area_idx & tonicAPs_signif_all==true)],...
        'sessions', {[sessions_all(area_idx & tonicAPs_signif_all==false); sessions_all(area_idx & tonicAPs_signif_all==true)]},...
        'unitID', {[units_all(area_idx & tonicAPs_signif_all==false); units_all(area_idx & tonicAPs_signif_all==true)]});
    fig = gcf;
    fig.UserData = Info;

end

% Plot amount of spikes belonging to a burst spontaneous vs. evoked
pre = vertcat(AP_belonging_to_bursts_pre{:});
post = vertcat(AP_belonging_to_bursts_post{:});

avgCellResponse.SpikesBelongingToBurst.pre = pre;
avgCellResponse.SpikesBelongingToBurst.post = post;

for ar = 1:numel(area_names)
    area_idx = contains(areas_all,area_names{ar});
    figure('Name',sprintf('SpikesBelongingToBurst_%s',area_names{ar}))
    hold on
    boxplot([pre(area_idx),post(area_idx)], 'Symbol', 'k.','Notch','on')
    parallelcoords([pre(area_idx & burstAPs_signif_all==false),post(area_idx & burstAPs_signif_all==false)], 'Color', 0.7*[1 1 1], 'LineStyle', '-');
    parallelcoords([pre(area_idx & burstAPs_signif_all==true),post(area_idx & burstAPs_signif_all==true)], 'Color', 0.2*[1 1 1], 'LineStyle', '-');
    axis auto

    % Paired sample t-test
    % consider using signrank for non-parametric testing
    p = signrank(pre(area_idx),post(area_idx));

    yt = get(gca, 'YTick');
    xt = get(gca, 'XTick');
    if p < 0.001
        plot([xt(1)+0.2,xt(2)-0.2], [1.1 1.1]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*1.125, '***','HorizontalAlignment','center')
    elseif p < 0.01
        plot([xt(1)+0.2,xt(2)-0.2], [1.1 1.1]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*1.125, '**','HorizontalAlignment','center')
    elseif p < 0.05
        plot([xt(1)+0.2,xt(2)-0.2], [1.1 1.1]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*1.125, '*','HorizontalAlignment','center')
    else
        plot([xt(1)+0.2,xt(2)-0.2], [1.1 1.1]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*1.125, 'n.s.','HorizontalAlignment','center')
    end

    title(sprintf('Amount of spikes belonging to a burst for given clusters from %s',area_names{ar}))
    ylabel('burst proportion')
    xticklabels({'spontaneous','evoked'})

    % Add a small pie chart, indicating the proportion of significantly
    % responding units
    axes('Position',[.75 .8 .1 .1]), box off
    piePlot = pie([sum(area_idx & burstAPs_signif_all==true) sum(area_idx & burstAPs_signif_all==false)]);
    patchHand = findobj(piePlot, 'Type', 'Patch');
    patchHand(2).FaceColor = 0.7*[1 1 1];
    patchHand(1).FaceColor = 0.2*[1 1 1];
    delete(findobj(piePlot,'Type','text'))

    % Save information on units and sessions for later retrieval
    Info = struct('data', [pre(area_idx & burstAPs_signif_all==false),post(area_idx & burstAPs_signif_all==false);...
        pre(area_idx & burstAPs_signif_all==true),post(area_idx & burstAPs_signif_all==true)],...
        'sessions', {[sessions_all(area_idx & burstAPs_signif_all==false); sessions_all(area_idx & burstAPs_signif_all==true)]},...
        'unitID', {[units_all(area_idx & burstAPs_signif_all==false); units_all(area_idx & burstAPs_signif_all==true)]});
    fig = gcf;
    fig.UserData = Info;

end

% Plot burst ratio per area spontaneous vs. evoked
pre = vertcat(burstRatio_pre{:});
post = vertcat(burstRatio_post{:});

avgCellResponse.BurstsToTonicRatio.pre = pre;
avgCellResponse.BurstsToTonicRatio.post = post;

for ar = 1:numel(area_names)
    area_idx = contains(areas_all,area_names{ar});
    figure('Name',sprintf('BurstsToTonicRatio_%s',area_names{ar}))
    hold on
    boxplot([pre(area_idx),post(area_idx)], 'Symbol', 'k.','Notch','on')
    parallelcoords([pre(area_idx),post(area_idx)], 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
        'Marker', '.', 'MarkerSize', 10);
    axis auto

    % Paired sample t-test
    % consider using signrank for non-parametric testing
    p = signrank(pre(area_idx),post(area_idx));

    yt = get(gca, 'YTick');
    xt = get(gca, 'XTick');
    if p < 0.001
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '***','HorizontalAlignment','center')
    elseif p < 0.01
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '**','HorizontalAlignment','center')
    elseif p < 0.05
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '*','HorizontalAlignment','center')
    else
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, 'n.s.','HorizontalAlignment','center')
    end

    title(sprintf('Burst ratio of clusters from %s',area_names{ar}))
    ylabel('burst ratio')
    xticklabels({'spontaneous','evoked'})
end

% Plot total number of bursts per trial per area spontaneous vs. evoked
pre = vertcat(burstsPerTrial_pre{:});
post = vertcat(burstsPerTrial_post{:});

avgCellResponse.BurstNum.pre = pre;
avgCellResponse.BurstNum.post = post;

for ar = 1:numel(area_names)
    area_idx = contains(areas_all,area_names{ar});
    figure('Name',sprintf('BurstNum_%s',area_names{ar}))
    hold on
    boxplot([pre(area_idx),post(area_idx)], 'Symbol', 'k.','Notch','on')
    parallelcoords([pre(area_idx),post(area_idx)], 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
        'Marker', '.', 'MarkerSize', 10);
    axis auto

    % Paired sample t-test
    % consider using signrank for non-parametric testing
    p = signrank(pre(area_idx),post(area_idx));

    yt = get(gca, 'YTick');
    xt = get(gca, 'XTick');
    if p < 0.001
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '***','HorizontalAlignment','center')
    elseif p < 0.01
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '**','HorizontalAlignment','center')
    elseif p < 0.05
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '*','HorizontalAlignment','center')
    else
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, 'n.s.','HorizontalAlignment','center')
    end

    title(sprintf('Bursts per trial of clusters from %s',area_names{ar}))
    ylabel('bursts')
    xticklabels({'spontaneous','evoked'})
end

% Plot average number of spikes per burst per area spontaneous vs. evoked
pre = vertcat(spikesPerBurst_pre{:});
post = vertcat(spikesPerBurst_post{:});

avgCellResponse.avgSpikesPerBurst.pre = pre;
avgCellResponse.avgSpikesPerBurst.post = post;

for ar = 1:numel(area_names)
    area_idx = contains(areas_all,area_names{ar});
    figure('Name',sprintf('avgSpikesPerBurst_%s',area_names{ar}))
    hold on
    boxplot([pre(area_idx),post(area_idx)], 'Symbol', 'k.','Notch','on')
    parallelcoords([pre(area_idx),post(area_idx)], 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
        'Marker', '.', 'MarkerSize', 10);
    axis auto

    % Paired sample t-test
    % consider using signrank for non-parametric testing
    p = signrank(pre(area_idx),post(area_idx));

    yt = get(gca, 'YTick');
    xt = get(gca, 'XTick');
    if p < 0.001
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '***','HorizontalAlignment','center')
    elseif p < 0.01
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '**','HorizontalAlignment','center')
    elseif p < 0.05
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '*','HorizontalAlignment','center')
    else
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, 'n.s.','HorizontalAlignment','center')
    end

    title(sprintf('Average number of spikes per burst of clusters from %s',area_names{ar}))
    ylabel('spikes per burst')
    xticklabels({'spontaneous','evoked'})
end

% Plot average burst frequency per area spontaneous vs. evoked
pre = vertcat(burstFrequency_pre{:});
post = vertcat(burstFrequency_post{:});

avgCellResponse.avgBurstFrequency.pre = pre;
avgCellResponse.avgBurstFrequency.post = post;

for ar = 1:numel(area_names)
    area_idx = contains(areas_all,area_names{ar});
    figure('Name',sprintf('avgBurstFrequency_%s',area_names{ar}))
    hold on
    boxplot([pre(area_idx),post(area_idx)], 'Symbol', 'k.','Notch','on')
    parallelcoords([pre(area_idx),post(area_idx)], 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
        'Marker', '.', 'MarkerSize', 10);
    axis auto

    % Paired sample t-test
    % consider using signrank for non-parametric testing
    p = signrank(pre(area_idx),post(area_idx));

    yt = get(gca, 'YTick');
    xt = get(gca, 'XTick');
    if p < 0.001
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '***','HorizontalAlignment','center')
    elseif p < 0.01
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '**','HorizontalAlignment','center')
    elseif p < 0.05
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, '*','HorizontalAlignment','center')
    else
        plot([xt(1)+0.2,xt(2)-0.2], [0.8 0.8]*max(yt), '-k'),  text(mean(xt([1 2])), max(yt)*0.825, 'n.s.','HorizontalAlignment','center')
    end

    title(sprintf('Average burst frequency of clusters from %s',area_names{ar}))
    ylabel('frequency [Hz]')
    xticklabels({'spontaneous','evoked'})
end

%% Save figures
% Destination folder for matlab .fig files
destfold = fullfile(cohortPath,'Analysis-Figures','Burstiness-Analysis',sessionDescription{:});
if exist(destfold,"dir") == 0
    mkdir(destfold)
end

% Save single unit information about tonic and burst behavior 
save(fullfile(destfold,responseFileName),'burstAPs_Info', 'tonicAPs_Info','avgCellResponse','fileVersion')

figHandles = handle(sort(double(findall(0, 'type', 'figure'))));

fprintf("Saving Figures.\n")
if targetEstimate
    figure_suffix = [condition, '_', cell2mat(join(trialType,'&')), '_withTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f')];
else
    figure_suffix = [condition, '_', cell2mat(join(trialType,'&')), '_withoutTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f')];
end

% Figure data to add some more information
Info = struct('statTest', 'signrank',...
    'spontWindow',spontWindow,...
    'respWindow',respWindow,...
    'trigOffset',trigOffset,...
    'targetEstimate',targetEstimate);

mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
for i = 1:numel(figHandles)
    destfile = sprintf('%s\\%s_%s.fig', destfold, figure_suffix, figHandles(i).Name);
    if ~isempty(figHandles(i).UserData)
        figHandles(i).UserData = mergestructs(Info,figHandles(i).UserData);
    else
        figHandles(i).UserData = Info;
    end
    savefig(figHandles(i), destfile);
end
fprintf("\nDone!\n\n")

%% Helper functions
function [put_in_idx, put_ex_idx] = getWaveformDistribution(clWaveforms,areaCell)
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

% Choose trough to peak cutoff of 300 µs for thalamic nuclei
idx = find(cellfun(@(x) ismember(x, {'VPM','POm'}),areaCell));
put_in_idx(idx) = trough2peak_usec(idx) < 300;
put_ex_idx(idx) = trough2peak_usec(idx) >= 300;

% Choose trough to peak cutoff of 350 µs for cortex and ZI
idx = find(cellfun(@(x) ismember(x, {'BC','ZIv'}),areaCell));
put_in_idx(idx) = trough2peak_usec(idx) < 350;
put_ex_idx(idx) = trough2peak_usec(idx) >= 350;
end
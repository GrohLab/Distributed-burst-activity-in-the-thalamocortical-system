%% Analyze stimulus reaction activity

% Define animal individuals to analyze

close all
clearvars
clc

% Access correct individual
scriptFullPath = matlab.desktop.editor.getActiveFilename();
load(regexprep(scriptFullPath, 'Scripts.*', 'Scripts\userDataPath.mat'), 'cohortPath');
try
    load(fullfile(cohortPath,'animalData.mat'))
catch
    fprintf('\nNo animalData.mat file loaded.\n')
end

cohortPath = uigetdir(cohortPath,'Select an Animal');
if ~cohortPath
    fprintf(2,'\nNo start path was selected.\n\n')
    return
elseif exist(cohortPath,'dir')==0
    fprintf(2,'\nYour start path is non existing.')
    fprintf(2,'\nChoose a different one.\n\n')
    return
else
    [~,animalName,~] = fileparts(cohortPath);
    if isempty(regexp(animalName,'#\d','once')) || regexp(animalName,'#\d') ~= 1
        fprintf(2,'\nYour start path is not an animal directory.')
        fprintf(2,'\nChoose a different one.\n\n')
        return
    end
end

% Get the relevant variables of each animal from animalData.mat
cohort_str = getCohort(cohortPath);
cohort_num = str2double(regexp(cohort_str,'\d*','match'));

cohort_animals = dir(fullfile(fileparts(cohortPath),'#*'));
animal_num = find(strcmp({cohort_animals.name},animalName));

area_names = {'BC','VPM','POm','ZIv'};

% Analyze included or excluded target estimates
answer = questdlg('Do you want to include presumably off-target tetrodes?', ...
    'Include off-target tetrodes', ...
    'Yes','No','No');
switch answer
    case 'Yes'
        targetEstimate = false;
    case 'No'
        targetEstimate = true;
        try
            load(fullfile(cohortPath,'targetHit.mat'),'targetHit')
        catch
            fprintf(2,'\nNo targetHit.mat file. Progress with unfiltered analysis...\n')
            targetEstimate = false;
        end
    otherwise
        % Return if no answer was picked
        return
end

fprintf('\nLooking for automatically curated files in ''%s'' ...\n',cohortPath)
FileInfo = dir(fullfile(cohortPath,'**\*automatedCuration'));
fprintf('Files collected!\n')

% Choose files to analyze
answer = listdlg('ListString',fileparts({FileInfo.folder}),...
    'PromptString','Choose sessions to analyze.',...
    'ListSize',[600 350]);
FileInfo = FileInfo(answer,:);

mean_fr = table('Size',[1,4],'VariableTypes',repmat({'cell'},1,4),'VariableNames',area_names);

condition_names = {'Reward','Punishment','Lick','WhiskerContact_left','WhiskerContact_right',...
    'Middlepoint','Random','WhiskerContact_onlyLeftFirst','WhiskerContact_onlyRightFirst','onlyFirstLick'};
% Condition to calculate responses over
answer = listdlg('PromptString','Which condition do you want to analyze?', ...
        'ListString',condition_names,...
        'ListSize',[300 200],'InitialValue',4,'SelectionMode','single');
chCond = condition_names{answer};

% All registered lick events are usually analyzed.
% To analyze only the first lick event of each lick series, set the
% respective flag to true.
if isequal(chCond,'Lick')
    lickFlag = questdlg('Do you only want to extract the first lick event in a series of lick events?', ...
        'Lick event extraction', ...
        'Yes','No, consider all lick events','Yes');
    if isequal(lickFlag,'Yes')
        lickFlag = true;
    else
        lickFlag = false;
    end
end

spontWindow = [-0.6,-0.4]; % spontaneous window in sec; default = [-0.6,-0.4]
respWindow = [0,0.2]; % response window in sec; default = [0,0.2]
trigOffset = false; % for response on 'offset', set to true; default = false

assert(round(diff(respWindow),3)==round(diff(spontWindow),3),...
    'Salience_function:WindowsNotSameLength',...
    'responseWindow and spontaneousWindow must be of the same length.')
windowLength = round(diff(spontWindow),3);

if isequal(chCond,'Lick') && lickFlag && targetEstimate
    outputFileName = ['StimulusResponse_onlyFirst', chCond, '_withTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f.mat')];
elseif targetEstimate
    outputFileName = ['StimulusResponse_', chCond, '_withTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f.mat')];
else
    outputFileName = ['StimulusResponse_', chCond, '_withoutTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f.mat')];
end

%% Calculate activity per unit in given windows

for ses = 1:height(FileInfo)
    try
        basePath = fullfile(FileInfo(ses).folder,FileInfo(ses).name);
        fprintf('\nAnalyzing session %d/%d : %s\n',ses,height(FileInfo),basePath)
        %         if targetEstimate && ...
        %                 exist(fullfile(basePath,outputFileName),'file')
        %             fprintf('%s was already analyzed. Moving on...\n',basePath)
        %             continue
        %         elseif ~targetEstimate && ...
        %                 exist(fullfile(basePath,outputFileName),'file')
        %             fprintf('%s was already analyzed. Moving on...\n',basePath)
        %             continue
        %         end

        tempDir = basePath;
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
        intanPath = fullfile(intanDir,rhdFile(1).name);
        sampNum = Intan_sampleNum(intanPath);

        % For session name
        stageCount = regexp(basePath,'P3.[0-9.]*','match','once');
        sesCount = regexp(basePath,'session\d*','match','once');
        sesName = ['stage',stageCount(4:end),'_',sesCount];

        qlMetrics = readtable(fullfile(basePath,'metrics_test.csv'));
        qlMetrics = qlMetrics(:,2:end);
        qlMetrics = removevars(qlMetrics,'epoch_name');

        [~,~,contamPct] = tsvread(fullfile(basePath,'cluster_ContamPct.tsv'));

        % Assign clusters to area and calculate mean z-scores
        try
            dataInfo = dir(fullfile(basePath,'*all_channels.mat'));
            load(fullfile(dataInfo(1).folder,dataInfo(1).name),'sortedData','fs')
        catch
            try
                [sortedData, fs] = importPhyFiles(basePath);
            catch
                fprintf(2,'Error importing the phy files into Matlab format\n')
                return
            end
        end

        % Sometimes empty clusters are registered, which are excluded here
        sortedData = sortedData(all(~cellfun(@isempty, sortedData),2),:);

        clAutoInfo = cell(height(sortedData),1);
        clAutoInfo(cell2mat(sortedData(:,3))==1) = deal({'good'});
        clAutoInfo(cell2mat(sortedData(:,3))==2) = deal({'mua'});
        clAutoInfo(cell2mat(sortedData(:,3))==3) = deal({'noise'});

        qlMetrics = addvars(qlMetrics,clAutoInfo,str2double(contamPct(2:end,2)),'NewVariableNames',{'group','contamPct'},'After','cluster_id');

        clInfo = getClusterInfo(fullfile(basePath,'cluster_info.tsv'));

        if targetEstimate
            include = ismember(qlMetrics.group,{'good','mua'}) & qlMetrics.isolation_distance > 15 & qlMetrics.isi_viol < 3 & ismember(clInfo.ch,targetHit);
        else
            include = ismember(qlMetrics.group,{'good','mua'}) & qlMetrics.isolation_distance > 15 & qlMetrics.isi_viol < 3;
        end
        qlMetrics = addvars(qlMetrics,include,'NewVariableNames','include','After','group');

        % Assign recording area to channels
        area = cell(height(clInfo),1);
        for cl = 1:height(clInfo)
            depth = clInfo.depth(cl);
            switch depth
                case 1
                    area{cl} = 'BC';
                case 1400
                    area{cl} = 'POm';
                case 1700
                    area{cl} = 'VPM';
                case 2400
                    area{cl} = 'ZIv';
                otherwise
                    area{cl} = 'NaN';
                    fprintf('\nUnregistered depth. Assigning NaN...\n')
            end
        end

        qlMetrics = addvars(qlMetrics,clInfo.ch,area,'NewVariableNames',{'ch','area'},'After','cluster_id');

        % Get the Triggers
        % Set warning to error in order to catch it
        s = warning('error', 'MATLAB:load:variableNotFound');
        condInfo = dir(fullfile(intanDir,'*analysis.mat'));
        if ~isempty(condInfo)
            try load(fullfile(condInfo.folder,condInfo.name),'Conditions','Triggers','update')
                if update < 5
                    [Conditions, Triggers] = getConditions(fileparts(intanDir));
                end
            catch
                [Conditions, Triggers] = getConditions(fileparts(intanDir));
            end
        else
            [Conditions, Triggers] = getConditions(fileparts(intanDir));
        end
        % Reset warning
        warning(s);

        % Check if all Conditions are the same in every session
        if ~isequal(sort(condition_names),sort(unique({Conditions(~contains({Conditions.name},'AllTriggers')).name})))
            fprintf(2, '\nThe Conditions in %s are different from the default Conditions:\n', basePath)
            fprintf(2, '\n- %s', Conditions.name)
            fprintf(2, '\n\nAborting...\n')
            return
        end

        % Calculate moving response index
        StimulusResponse = table('Size',[3,4],'VariableTypes',repmat({'cell'},1,4),'VariableNames',area_names,'RowNames',{'unitId','stimulusFr','baselineFr'});
        chCond_trig = sort(Conditions(find(cellfun(@(x) isequal(x,chCond), {Conditions.name}),1)).Triggers)./fs; % triggers in sec

        load(fullfile(fileparts(intanDir),'videos\HispeedTrials.mat'),'HispeedTrials')
        % For logical indexing NaN values must be set to 0
        HispeedTrials.Lick(isnan(HispeedTrials.Lick))=0;

        if isequal(chCond,'Lick') && lickFlag
            % This returns the first lick event after every middle point
            idx = unique(sort(cell2mat(arrayfun(@(x) find(chCond_trig(:,1)*1000 > x, 1), HispeedTrials.PreviousMP, 'UniformOutput', false))));
            chCond_trig = chCond_trig(idx,:);
        end

        Go_NoGo_Neutral = nan(size(chCond_trig,1),1);
        Wide_Narrow_Intermediate = nan(size(chCond_trig,1),1);
        Lick = nan(size(chCond_trig,1),1);

        % Transcribes trial types into aperture settings
        % Only applicable for cohort 12
        fprintf('\nBe aware that the assignment of aperture setting, i.e., narrow, wide, or intermediate,\n')
        fprintf('is based on the stage settings of cohort 12 go = wide and no-go = narrow in all stages,\n')
        fprintf('except for the rule-switch and extinction stage.\n')
        for ar = 1:numel(area_names)
            StimulusResponse{'unitId',area_names{ar}}{:} = qlMetrics.cluster_id(ismember(qlMetrics.area,area_names{ar}) & qlMetrics.include);
            StimulusResponse{'stimulusFr',area_names{ar}}{:} = NaN(size(StimulusResponse{'unitId',area_names{ar}}{:},1), size(chCond_trig,1));
            StimulusResponse{'baselineFr',area_names{ar}}{:} = NaN(size(StimulusResponse{'unitId',area_names{ar}}{:},1), size(chCond_trig,1));

            for trial = 1:size(chCond_trig,1)
                for unit = 1:size(StimulusResponse{'unitId',area_names{ar}}{:},1)
                    pickedId = StimulusResponse{'unitId',area_names{ar}}{:}(unit);
                    pickedSpikes = sortedData{ismember(sortedData(:,1),num2str(pickedId)),2};

                    spontFr = sum(pickedSpikes >= chCond_trig(trial,trigOffset+1)+spontWindow(1) & pickedSpikes <= chCond_trig(trial,trigOffset+1)+spontWindow(2))/windowLength;
                    if chCond_trig(trial,trigOffset+1)+respWindow(2)*fs > sampNum
                        respFr = NaN;
                    else
                        respFr = sum(pickedSpikes >= chCond_trig(trial,trigOffset+1)+respWindow(1) & pickedSpikes <= chCond_trig(trial,trigOffset+1)+respWindow(2))/windowLength;
                    end

                    StimulusResponse{'stimulusFr',area_names{ar}}{:}(unit, trial) = respFr;
                    StimulusResponse{'baselineFr',area_names{ar}}{:}(unit, trial) = spontFr;
                end

                diffTrigger = chCond_trig(trial,1)*1000 - HispeedTrials.PreviousMP;
                diffTrigger(diffTrigger<0) = Inf;
                [~,idx] = min(diffTrigger);

                if ~isempty(idx)
                    Go_NoGo_Neutral(trial) = HispeedTrials.Go_NoGo_Neutral(idx);
                    Lick(trial) = HispeedTrials.Lick(idx);
                    if contains(basePath,'extinction') || contains(basePath,'ruleswitch')
                        switch HispeedTrials.Go_NoGo_Neutral_settingBased(idx)
                            case 1
                                Wide_Narrow_Intermediate(trial) = 2; % "Go" = "Narrow"
                            case 2
                                Wide_Narrow_Intermediate(trial) = 1; % "No-Go" = "Wide"
                            case 3
                                Wide_Narrow_Intermediate(trial) = 3;
                            otherwise
                                Wide_Narrow_Intermediate(trial) = NaN;
                        end
                    else
                        Wide_Narrow_Intermediate(trial) = HispeedTrials.Go_NoGo_Neutral_settingBased(idx); % "Go" = "Wide" and "No-Go" = "Narrow"
                    end
                end
            end
        end

        % Save data in table with structure (n = number of units | m = number of trials):
        %                       Area1   Area2   Area3
        % unitId                n x 1   ...
        % stimulusFr            n x m   ...
        % baselineFr            n x m   ...
        % syntax: "StimulusResponse_stimulus_withTargetEstimate_responseWind_baselineWind.mat"
        % Go_NoGo_Neutral = m x 1
        % Lick            = m x 1
        save(fullfile(basePath,outputFileName),'StimulusResponse','Go_NoGo_Neutral','Wide_Narrow_Intermediate','Lick')

    catch except
        fprintf(2, '\nError in: %s\n',fullfile(FileInfo(ses).folder,FileInfo(ses).name))
        fprintf(2, 'Error in Line: %d\n',except.stack.line)
        fprintf(2, 'Error: %s\n',except.message);
    end
end

fprintf('\nAll files analyzed!\n\n')


%% Single Cell Reactivity
% Split cells into groups of aperture responding units
% Reponses can either be single-targeted ("wide" or "narrow"), "both" or
% "none"

% Define animal individuals to analyze

close all
clearvars
clc

% This is a quick solution for providing only the automatedCuration
% directories. For more general usage, create an "analyzed_files.mat" file
% for the StimulusResponse_Analysis script.
scriptFullPath = matlab.desktop.editor.getActiveFilename();
load(regexprep(scriptFullPath, 'Scripts.*', 'Scripts\\userDataPath.mat'), 'cohortPath');
load(fullfile(cohortPath,'allFiles.mat'),'FileInfo')
updatedFolders = fullfile(cohortPath,{FileInfo.folder});
[FileInfo.folder] = deal(updatedFolders{:});
AllFiles = fullfile(fileparts({FileInfo.folder}), 'intan-signals\automatedCuration')';

% Access correct individual
try
    load(fullfile(cohortPath, 'animalData.mat'))
catch
    fprintf('\nNo animalData.mat file loaded.\n')
end

% Define brain areas
area_names = {'BC','VPM','POm','ZIv'};
area_colors = {'#377eb8','#4daf4a','#984ea3','#ff7f00'};

answer = questdlg('Would you like to implement a time shift for the spike data (for Matthias-Martin-Experiment)?', ...
	'Spike time shift', ...
	'Static time shift','Dynamic time shift','No time shift','No time shift');

if isequal(answer,'Static time shift')
    timeShift.type = 'static';
    
    answer = inputdlg(sprintf('Enter the desired static time shift (in msec) for the spike data.\n"None" if you don''t want to shift data points.'),...
        'Spike time shift',[1 50],{'None'});
    if str2double(answer{:})==0 || isnan(str2double(answer{:}))
        timeShift.type = 'none';
        timeShift.value = 0;
    else
        timeShift.value = str2double(answer{:});
    end
elseif isequal(answer,'Dynamic time shift')
    timeShift.type = 'dynamic';
    
    answer = inputdlg(sprintf('Enter the desired dynamic time shift (in msec per second of recording) for the spike data.\n"None" if you don''t want to shift data points.'),...
        'Spike time shift',[1 50],{'None'});
    if str2double(answer{:})==0 || isnan(str2double(answer{:}))
        timeShift.type = 'none';
        timeShift.value = 0;
    else
        timeShift.value = str2double(answer{:});
    end
else
    timeShift.type = 'none';
    timeShift.value = 0;
end

% Pick condition to calculate responses over
condition_names = {'Reward', 'Punishment', 'Lick', 'onlyFirstLick',...
    'WhiskerContact_left', 'WhiskerContact_right','Middlepoint',...
    'Random', 'WhiskerContact_onlyLeftFirst', 'WhiskerContact_onlyRightFirst'};
[condPick,tf] = listdlg('Name','Condition choice','PromptString',{'Which condition do you want to analyze?'},...
    'ListString',condition_names,'ListSize',[250 250],'InitialValue',5);
if tf
    chCond = condition_names{condPick};
else
    return
end

spontWindow = [-0.6,-0.4]; % spontaneous window in sec (default: [-0.6,-0.4])
respWindow = [0,0.2]; % response window in sec (default: [0,0.2])
trigOffset = false; % for response on 'offset', set to true (default: false)

assert(round(diff(respWindow),3)==round(diff(spontWindow),3),...
    'Salience_function:WindowsNotSameLength',...
    'responseWindow and spontaneousWindow must be of the same length.')
windowLength = round(diff(spontWindow),3);

outputFileName = ['SingleCellReactivity_', chCond, '_withTargetEstimate_', num2str(respWindow,'respond%.3f-%.3f_'), num2str(spontWindow,'baseline%.3f-%.3f.mat')];

if isequal(timeShift.type,'static')
    outputFileName = strrep(outputFileName,'.mat',sprintf('_timeShift%gmsec.mat',timeShift.value));
elseif isequal(timeShift.type,'dynamic')
    outputFileName = strrep(outputFileName,'.mat',sprintf('_dynamicTimeShift%gmsec.mat',timeShift.value));
end

% Pick sessions to analyze
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
    stageDescription = inputdlg('Enter a session description for saving the files:','Session description',[1 100],{definput});
else
    stageDescription = inputdlg('Enter a session description for saving the files:','Session description',[1 100]);
end
% Strip the string from any special characters in the beginning or end
stageDescription = regexprep(stageDescription,'[^a-zA-Z0-9]*$','');
stageDescription = regexprep(stageDescription,'^[^a-zA-Z0-9]*','');
stageDescription = strrep(stageDescription,' ','');

%% Calculate amount of responsive cells

if contains(chCond,{'Lick','Whisker'})
    condNames = {'wide','narrow','intermediate'};
else
    condNames = {'response'};
end
Response_Table = table('Size',[numel(condNames),numel(area_names)],'VariableTypes',repmat({'cell'},1,4),'VariableNames',area_names,'RowNames',condNames);

for ses = 1:height(fileSelection)

    basePath = fileSelection{ses};
    str_idx = regexp(basePath,'#\d*','end');
    animalPath = basePath(1:str_idx);

    fprintf('\nAnalyzing session %d/%d : %s\n',ses,height(fileSelection),basePath)
    if isfile(fullfile(basePath,outputFileName))
        load(fullfile(basePath,outputFileName),'SingleCellResponse')

        if contains(chCond,{'Lick','Whisker'})
            % Create table with all sessions
            for ar = 1:numel(area_names)
                idx = ismember(SingleCellResponse.area,area_names{ar}) & SingleCellResponse.include;
                Response_Table{'wide',area_names{ar}}{:} = [Response_Table{'wide',area_names{ar}}{:}; SingleCellResponse.respWide(idx)];
                Response_Table{'narrow',area_names{ar}}{:} = [Response_Table{'narrow',area_names{ar}}{:}; SingleCellResponse.respNarrow(idx)];
                Response_Table{'intermediate',area_names{ar}}{:} = [Response_Table{'intermediate',area_names{ar}}{:}; SingleCellResponse.respInterm(idx)];
            end
        else
            Response_Table{'response',area_names{ar}}{:} = [Response_Table{'response',area_names{ar}}{:}; SingleCellResponse.respStim(idx)];
        end
    else

        try

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

            % For session name
            stageCount = regexp(basePath,'P3.[0-9.]*','match','once');
            sesCount = regexp(basePath,'session\d*','match','once');
            sesName = ['stage',stageCount(4:end),'_',sesCount];

            SingleCellResponse = readtable(fullfile(basePath,'metrics_test.csv'));
            SingleCellResponse = SingleCellResponse(:,2:end);
            SingleCellResponse = removevars(SingleCellResponse,'epoch_name');

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

            SingleCellResponse = addvars(SingleCellResponse,clAutoInfo,str2double(contamPct(2:end,2)),'NewVariableNames',{'group','contamPct'},'After','cluster_id');

            clInfo = getClusterInfo(fullfile(basePath,'cluster_info.tsv'));

            try
                load(fullfile(animalPath,'targetHit.mat'),'targetHit')
            catch
                fprintf(2,'\nNo targetHit.mat file in path: "%s"\n',animalPath)
                return
            end
            include = ismember(SingleCellResponse.group,{'good','mua'}) & SingleCellResponse.isolation_distance > 15 & SingleCellResponse.isi_viol < 3 & ismember(clInfo.ch,targetHit);
            SingleCellResponse = addvars(SingleCellResponse,include,'NewVariableNames','include','After','group');

            SingleCellResponse = SingleCellResponse(:,1:5);

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

            SingleCellResponse = addvars(SingleCellResponse,clInfo.ch,area,'NewVariableNames',{'ch','area'},'After','cluster_id');

            % For whisker triggers calculate responses on every aperture state
            % For all other triggers only for the stimulus itself
            if contains(chCond,{'Lick','Whisker'})
                SingleCellResponse = addvars(SingleCellResponse,false(height(SingleCellResponse),1),...
                    false(height(SingleCellResponse),1),false(height(SingleCellResponse),1),'NewVariableNames',{'respWide','respNarrow','respInterm'},'After','include');
                SingleCellResponse = SingleCellResponse(:,1:8);

                stateList = {'wide','narrow','intermediate'};
                for state = 1:numel(stateList)
                    try
                        % Check for significant responses
                        DE_Salience_function(basePath,'responseWindow',respWindow,...
                            'spontaneousWindow',spontWindow,'condition',chCond,...
                            'adjustConditions',stateList(state),'interruptAfterPSTH',true,...
                            'timeShiftType',timeShift.type,'timeShiftValue',timeShift.value);
                        close all
                        signClusters = goodClusters(significantResponse);
                        SingleCellResponse{ismember(SingleCellResponse.cluster_id,cellfun(@str2num, signClusters)),5+state} = true;
                    catch
                    end
                end
            else
                SingleCellResponse = addvars(SingleCellResponse,false(height(SingleCellResponse),1),'NewVariableNames','respStim','After','include');
                SingleCellResponse = SingleCellResponse(:,1:6);
            end

            if contains(chCond,{'Lick','Whisker'})
                % Create table with all sessions
                for ar = 1:numel(area_names)
                    idx = ismember(SingleCellResponse.area,area_names{ar}) & SingleCellResponse.include;
                    Response_Table{'wide',area_names{ar}}{:} = [Response_Table{'wide',area_names{ar}}{:}; SingleCellResponse.respWide(idx)];
                    Response_Table{'narrow',area_names{ar}}{:} = [Response_Table{'narrow',area_names{ar}}{:}; SingleCellResponse.respNarrow(idx)];
                    Response_Table{'intermediate',area_names{ar}}{:} = [Response_Table{'intermediate',area_names{ar}}{:}; SingleCellResponse.respInterm(idx)];
                end
            else
                Response_Table{'response',area_names{ar}}{:} = [Response_Table{'response',area_names{ar}}{:}; SingleCellResponse.respStim(idx)];
            end

            % Save data as table
            save(fullfile(basePath,outputFileName),'SingleCellResponse')

        catch except
            fprintf(2, '\nError in: %s\n',fileSelection{ses})
            fprintf(2, 'Error in Line: %d\n',except.stack.line)
            fprintf(2, 'Error: %s\n',except.message);
        end

    end
end

%% Calculate amount of responsive cells for chosen sessions

% Pie charts for each area with mean responsiveness
figure('Name',sprintf('%s_ResponiveCells_PieChart',stageDescription{:}))

if isequal(timeShift.type,'static')
    userTitle = inputdlg('Enter a super title for the pie charts:','Figure Title',[1 100],{[stageDescription{:}, sprintf(' - %g msec time shift',timeShift.value)]});
elseif isequal(timeShift.type,'dynamic')
    userTitle = inputdlg('Enter a super title for the pie charts:','Figure Title',[1 100],{[stageDescription{:}, sprintf(' - %g msec dynamic time shift',timeShift.value)]});
else
    userTitle = inputdlg('Enter a super title for the pie charts:','Figure Title',[1 100],stageDescription);
end

sgtitle(userTitle{:})
for ar = 1:numel(area_names)
    subplot(2,2,ar)
    if contains(chCond,{'Lick','Whisker'})
        respLabels = ['no response',condNames,'wide & narrow','wide & intermediate','narrow & intermediate','all apertures'];
        % no response
        respVals(1) = sum(~Response_Table{'wide',area_names{ar}}{:} & ~Response_Table{'narrow',area_names{ar}}{:} & ~Response_Table{'intermediate',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});
        % wide response
        respVals(2) = sum(Response_Table{'wide',area_names{ar}}{:} & ~Response_Table{'narrow',area_names{ar}}{:} & ~Response_Table{'intermediate',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});
        % narrow response
        respVals(3) = sum(~Response_Table{'wide',area_names{ar}}{:} & Response_Table{'narrow',area_names{ar}}{:} & ~Response_Table{'intermediate',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});
        % intermediate response
        respVals(4) = sum(~Response_Table{'wide',area_names{ar}}{:} & ~Response_Table{'narrow',area_names{ar}}{:} & Response_Table{'intermediate',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});
        % wide & narrow response
        respVals(5) = sum(Response_Table{'wide',area_names{ar}}{:} & Response_Table{'narrow',area_names{ar}}{:} & ~Response_Table{'intermediate',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});
        % wide & intermediate response
        respVals(6) = sum(Response_Table{'wide',area_names{ar}}{:} & ~Response_Table{'narrow',area_names{ar}}{:} & Response_Table{'intermediate',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});
        % narrow & intermediate response
        respVals(7) = sum(~Response_Table{'wide',area_names{ar}}{:} & Response_Table{'narrow',area_names{ar}}{:} & Response_Table{'intermediate',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});
        % all apertures response
        respVals(8) = sum(Response_Table{'wide',area_names{ar}}{:} & Response_Table{'narrow',area_names{ar}}{:} & Response_Table{'intermediate',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});

        nonZero = false(1,8);
    else
        respLabels = {'no response','response'};
        % no response
        respVals(1) = sum(~Response_Table{'response',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});
        % response
        respVals(2) = sum(Response_Table{'response',area_names{ar}}{:})/...
            numel(Response_Table{'wide',area_names{ar}}{:});

        nonZero = false(1,2);
    end
    if all(~isnan(respVals))
        p = pie(respVals);
        fprintf('\nArea: %s| Perc of responsive cells: %.1f | n of cells: %i\n',area_names{ar},(1-respVals(1))*100,numel(Response_Table{1,ar}{:}))
        nonZero = logical(respVals) | nonZero;
        delete(findobj(p,'String','0%')); % Delete 0 displays
    else
        axis off
        text(0.5,0.5,'No values available.',...
            'HorizontalAlignment','center','VerticalAlignment','middle')
    end
    title(area_names{ar})

    patchHand = findobj(p, 'Type', 'Patch');
    if contains(chCond,{'Lick','Whisker'})
        patchHand(1).FaceColor = '#e6e6e6'; % no response
        patchHand(2).FaceColor = '#2fb62f'; % wide response
        patchHand(2).FaceAlpha = 0.5;
        patchHand(3).FaceColor = '#80b3ff'; % narrow response
        patchHand(3).FaceAlpha = 0.5;
        patchHand(4).FaceColor = '#cc00ff'; % intermediate response
        patchHand(4).FaceAlpha = 0.5;
        patchHand(5).FaceColor = '#1a75ff'; % wide & narrow response
        patchHand(5).FaceAlpha = 0.5;
        patchHand(6).FaceColor = '#ff9933'; % wide & intermediate response
        patchHand(6).FaceAlpha = 0.5;
        patchHand(7).FaceColor = '#006699'; % narrow & intermediate response
        patchHand(7).FaceAlpha = 0.5;
        patchHand(8).FaceColor = '#666633'; % all apertures response
        patchHand(8).FaceAlpha = 0.5;
    else
        patchHand(1).FaceColor = '#e6e6e6'; % no response
        patchHand(2).FaceColor = '#2fb62f'; % response
        patchHand(2).FaceAlpha = 0.5;
    end
end
lgd = legend(patchHand(nonZero),respLabels(nonZero),'Box','off');
lgd.Position = [0.1 0.465 0.1 0.1];

%% Save figure

% Destination folder for matlab .fig files
destfold = fullfile(cohortPath, 'Analysis-Figures','SingleCellReactivity-Analysis');
if exist(destfold,"dir") == 0
    mkdir(destfold)
end

figHandles = handle(sort(double(findall(0, 'type', 'figure'))));

for i = 1:numel(figHandles)
    destfile = fullfile(destfold,sprintf('%s_ResponiveCells_PieChart_%s',stageDescription{:},strrep(outputFileName,'.mat','.fig')));
    savefig(figHandles(i), destfile);
end

fprintf("Done!\n")

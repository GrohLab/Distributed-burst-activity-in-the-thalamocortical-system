%% Correlating decoding accuracy with number of input units
% Script that runs multiple decoding analyses with increasing number of
% units. By bootstrapping each step, one can define the average decoding
% accuracy of each step.
close all; clearvars; clc
scriptFullPath = matlab.desktop.editor.getActiveFilename();
load(regexprep(scriptFullPath, 'Scripts.*', 'Scripts\\userDataPath.mat'), 'cohortPath');

% Choose sessions to merge together
try
    load(fullfile(cohortPath,'animalData.mat'))
    load(fullfile(cohortPath,'allFiles.mat'),'FileInfo')
    updatedFolders = fullfile(cohortPath,{FileInfo.folder});
    [FileInfo.folder] = deal(updatedFolders{:});
catch
end

spontWindow = [-600,-400]; % spontaneous window in msec (default: [-600,-400])
respWindow = [0,200]; % response window in msec (default: [0,200])

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
    % Output folder
    sessionDescription = inputdlg('Enter a session description for saving the files:','Output Folder',[1 100],{definput});
    % Raster data folder
    rasterDataName = inputdlg('Name of raster data folder:','Raster Data Folder',[1 100],{definput});
    raster_data_dir = fullfile(cohortPath,'NeuralDecoding',rasterDataName{:});
else
    % Output folder
    sessionDescription = inputdlg('Enter a session description for saving the files:','Output Folder',[1 100]);
    % Raster data folder
    rasterDataName = inputdlg('Name of raster data folder:','Raster Data Folder',[1 100]);
    raster_data_dir = fullfile(cohortPath,'NeuralDecoding',rasterDataName{:});
end

% Choose condition to analyze
% Available conditions: 'Reward', 'Punishment', 'Lick','onlyFirstLick',
% 'WhiskerContact_left', 'WhiskerContact_right',
% 'WhiskerContact_onlyLeftFirst','WhiskerContact_onlyRightFirst'
chCond = 'WhiskerContact_left';

% Specify which trial types you want to analyze
trialType = {'allTrials','onlyGo','onlyNoGo','onlyNeutral','onlyNarrow','onlyWide','onlyIntermediate','onlyLick','onlyNoLick'};
[trial_idx,tf] = listdlg('PromptString',{'Specify trial type'},...
    'ListSize',[300 200],...
    'ListString',trialType);
if tf == 0
    % If nothing is picked, move on with unfiltered trial analysis
    trial_idx = 1;
end
trialType = trialType{trial_idx};

% Create raster_data.mat files
create_raster_data_files(fileSelection, raster_data_dir, chCond, trialType);

%% Set decoder parameters
% Define brain areas
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
    'Callback',{@doneExe,f,ui_field,area_names});

waitfor(doneButton)

answer = questdlg('What do you want to decode?','Decoding object',...
    'Go - No-go - Neutral','Wide - Narrow - Interm.','Lick - No-lick','Wide - Narrow - Interm.');
if isequal(answer, 'Go - No-go - Neutral')
    classifier_labels = 'Go_NoGo_Neutral';
    if isequal(trialType,'alltrials') 
        trialCompare = {'narrow','wide'};
    elseif isequal(trialType,'onlyLick')
        trialCompare = {'narrow&lick','wide&lick'};
    elseif isequal(trialType,'onlyNoLick')
        trialCompare = {'narrow&no-lick','wide&no-lick'};
    else
        trialCompare = {''};
    end
elseif isequal(answer, 'Wide - Narrow - Interm.')
    classifier_labels = 'Wide_Narrow_Intermediate';
    if isequal(trialType,'alltrials') 
        trialCompare = {'narrow','wide'};
    elseif isequal(trialType,'onlyLick')
        trialCompare = {'narrow&lick','wide&lick'};
    elseif isequal(trialType,'onlyNoLick')
        trialCompare = {'narrow&no-lick','wide&no-lick'};
    else
        trialCompare = {''};
    end
elseif isequal(answer, 'Lick - No-lick')
    classifier_labels = 'Lick';
    if isequal(trialType,'onlyNarrow')
        trialCompare = {'narrow&lick','narrow&no-lick'};
    elseif isequal(trialType,'onlyWide')
        trialCompare = {'wide&lick','wide&no-lick'};
    elseif isequal(trialType,'onlyIntermediate')
        trialCompare = {'intermediate&lick','intermediate&no-lick'};
    else
        trialCompare = {''};
    end
end

% Define the classifier to use
% 'max_correlation_coefficient_CL','poisson_naive_bayes_CL','libsvm_CL'
classifierName = 'libsvm_CL';

% Set a desired number of splits. A split of 10 means that 9 repetitions of
% each event are used for training and 1 example is used for testing.
% To get reasonable results you usually need at least 5 repetitions of each
% event (i.e., at least 5 splits)
splitNumber = inputdlg(sprintf('Enter the desired number of splits\n(NaN if you want to estimate it\nbased on the unit and event count):'),'Split number',[1 50],{'NaN'});
splitNumber = str2double(splitNumber);

% Define pool of cells, you want to choose from:
% 'All cells', 'Only responding with burst bias', 
% 'Only responding with tonic bias', 'Only touch-modulated cells'
cellSpecs = listdlg('PromptString','Which cells would you like to analyze?', ...
    'Name','Cell types','ListString',{'All cells', 'Only responding with burst bias', 'Only responding with tonic bias',...
    'Only touch-modulated cells'}, ...
    'ListSize',[250,150],'SelectionMode','single');

if ismember(cellSpecs,[2,3,4]) % For tonic and burst modulation
    trialCompareList = {'all', 'wide', 'wide&lick', 'wide&no-lick',...
        'narrow', 'narrow&lick', 'narrow&no-lick',...
        'intermediate', 'intermediate&lick', 'intermediate&no-lick'};
    if isempty(trialCompare{1})
        initVals = [2 5];
    else
        initVals = [find(isequal(trialCompareList,trialCompare{1})) ...
            find(isequal(trialCompareList,trialCompare{2}))];
    end
    [trialPick, tf] = listdlg('ListString',trialCompareList,...
        'PromptString','Which trial types do you want to compare?',...
        'SelectionMode','multiple','ListSize', [250 250],'InitialValue',initVals);
    if tf==0
        return
    end
    assert(numel(trialPick)==2,...
        'singleUniteBurstiness:selectInputNumber',...
        'You have to pick 2 trial types.')
    trialCompare = trialCompareList(trialPick);
    
    % Get the responseTypes.mat file (narrow&wide) to assess touch responses 
    responseDir = fullfile(cohortPath,'Analysis-Figures\Burstiness-Scatter');
    responseFile = fullfile(responseDir,definput,sprintf('responseTypes_%s_%s.mat',chCond,strjoin(trialCompare,'&')));

    if ~exist(responseFile,'file')
        % Flip the trial types, to see if that file does exist
        if exist(fullfile(responseDir,definput,sprintf('responseTypes_%s_%s.mat',chCond,strjoin(flip(trialCompare),'&'))),'file')
            responseFile = fullfile(responseDir,definput,sprintf('responseTypes_%s_%s.mat',chCond,strjoin(flip(trialCompare),'&')));
            load(responseFile,'responseTypes')
        else
            fprintf('\nNo responseTypes.mat file found for trial types %s. Run the ApertureResponseTypes.m script first.\n', strjoin(trialCompare,' & '))
            return
        end
    else
        load(responseFile,'responseTypes')
    end
    
    % Filter for picked waveforms
    idx = false(height(responseTypes),1);
    for ar = 1:numel(area_names)
        if contains(area_names{ar},'RS')
            idx = idx | (strncmp(responseTypes.Area,area_names{ar},2) & responseTypes.PutExcitatory);
        elseif contains(area_names{ar},'FS')
            idx = idx | (strncmp(responseTypes.Area,area_names{ar},2) & ~responseTypes.PutExcitatory);
        else % All waveforms
            idx = idx | strncmp(responseTypes.Area,area_names{ar},2);
        end
    end
    responseTypes = responseTypes(idx,:);

    % Display how many cells of that cell type there are
    fprintf('\n%i cells with burst bias in stage "%s"',sum(responseTypes.BurstSign),definput)
    fprintf('\n     - %i cells in BC',sum(responseTypes.BurstSign & contains(responseTypes.Area,'BC')))
    fprintf('\n     - %i cells in VPM',sum(responseTypes.BurstSign & contains(responseTypes.Area,'VPM')))
    fprintf('\n     - %i cells in POm',sum(responseTypes.BurstSign & contains(responseTypes.Area,'POm')))
    fprintf('\n     - %i cells in ZIv\n',sum(responseTypes.BurstSign & contains(responseTypes.Area,'ZIv')))
    
    % Display how many cells of that cell type there are
    fprintf('\n%i cells with tonic bias in stage "%s"',sum(responseTypes.TonicSign),definput)
    fprintf('\n     - %i cells in BC',sum(responseTypes.TonicSign & contains(responseTypes.Area,'BC')))
    fprintf('\n     - %i cells in VPM',sum(responseTypes.TonicSign & contains(responseTypes.Area,'VPM')))
    fprintf('\n     - %i cells in POm',sum(responseTypes.TonicSign & contains(responseTypes.Area,'POm')))
    fprintf('\n     - %i cells in ZIv\n',sum(responseTypes.TonicSign & contains(responseTypes.Area,'ZIv')))

    % Display how many cells of that cell type there are
    fprintf('\n%i cells with touch modulation (%.2f%% of all cells) in stage "%s"',sum(responseTypes.TouchModul),100*sum(responseTypes.TouchModul)/height(responseTypes),definput)
    fprintf('\n     - %i cells in BC (%.2f%% of BC cells)',sum(responseTypes.TouchModul & contains(responseTypes.Area,'BC')), ...
        100*sum(responseTypes.TouchModul & contains(responseTypes.Area,'BC'))/sum(contains(responseTypes.Area,'BC')))
    fprintf('\n     - %i cells in VPM (%.2f%% of VPM cells)',sum(responseTypes.TouchModul & contains(responseTypes.Area,'VPM')), ...
        100*sum(responseTypes.TouchModul & contains(responseTypes.Area,'VPM'))/sum(contains(responseTypes.Area,'VPM')))
    fprintf('\n     - %i cells in POm (%.2f%% of POm cells)',sum(responseTypes.TouchModul & contains(responseTypes.Area,'POm')), ...
        100*sum(responseTypes.TouchModul & contains(responseTypes.Area,'POm'))/sum(contains(responseTypes.Area,'POm')))
    fprintf('\n     - %i cells in ZIv (%.2f%% of ZIv cells)\n',sum(responseTypes.TouchModul & contains(responseTypes.Area,'ZIv')), ...
        100*sum(responseTypes.TouchModul & contains(responseTypes.Area,'ZIv'))/sum(contains(responseTypes.Area,'ZIv')))
    switch cellSpecs
        case 2
            unitType = 'burstUnits';
        case 3
            unitType = 'tonicUnits';
        case 4
            unitType = 'touchModul';
    end

else
    unitType = 'allUnits';
end

% Define the number of bootstraps and number of runs (i.e., unit number)
prompt = {'Enter desired unit number of analyze (2 inputs define an array of increasing values):','Enter desired bootstrap number per run:'};
answer = inputdlg(prompt,'Decoding parameters',[1 60; 1 60],{'1 80','20'});
numbers = regexp(answer, '\d+', 'match');
numbers = cellfun(@str2double, numbers, 'UniformOutput', false);
assert(all(cellfun(@isnumeric, numbers)),'Inputs must be numeric.')
if numel(numbers{1})==2
    % If two numbers were given, treat them as an array
    numRuns = numbers{1}(1):numbers{1}(2);
else
    numRuns = numbers{1};
end
numBootstrap = numbers{2};

%% Run the decoder

answer = questdlg('Do you only want to run all areas combined (to speed up the process)?','All areas',...
    'Yes','No, run everything','No, do not analyse combined areas','Yes');
onlyCombined = false;
dontRunCombined = false;
if isequal(answer, 'Yes')
    onlyCombined = true;
elseif isequal(answer, 'No, do not analyse combined areas')
    dontRunCombined = true;
end

for ar = 1:numel(area_names)+1
    if ar~=3 % Only VPM
        continue
    end
    if dontRunCombined && ar==1
        continue
    elseif onlyCombined && ar>1
        continue
    end
    if ar==1
        areaIdx = true(height(responseTypes),1);
        tempAreaName = 'allAreas';
    else
        tempAreaName = area_names{ar-1};
        areaIdx = strncmp(responseTypes.Area,tempAreaName,2);
    end
    switch cellSpecs
        case 1 % all cells
            maxCells = sum(areaIdx);
        case 2 % only burst
            maxCells = sum(responseTypes.BurstSign & areaIdx);
        case 3 % only tonic
            maxCells = sum(responseTypes.TonicSign & areaIdx);
        case 4 % touch-modulated
            maxCells = sum(responseTypes.TouchModul & ~responseTypes.BurstSign & areaIdx);
    end

    for run = numRuns
        % Run neural decoding toolbox on respective cell types
        if run > maxCells
            continue
        end
        
        fprintf('\nSubsampling %i %s cells from %i (run %i of %i)\n',...
            run,tempAreaName,maxCells,find(run==numRuns),min([numel(numRuns),maxCells]))
        
        % Save bootstrapped files in a subfolder for each cell count
        outputFolder = fullfile(cohortPath,'Decoding-with-increasing-cellNum',sessionDescription{:},tempAreaName,sprintf('%03icells',run));
        if exist(outputFolder,"dir") == 0
            mkdir(outputFolder)
        end

        % Check how many bootstraps are already analyzed
        % If bootstraps are already available, built upon that
        FileInfo = dir(fullfile(outputFolder,sprintf('Binned_data_results_%s_allSpikes_%s_%s_allTrials_*',chCond,unitType,tempAreaName)));
        if isempty(FileInfo)
            alreadyAnal = 0;
        elseif height(FileInfo)>=numBootstrap
            continue
        else
            alreadyAnal = height(FileInfo);
        end

        for bootstrap = alreadyAnal+1:numBootstrap
            % Run runNeuralDecodingToolbox 
            fprintf('\nStarting bootstrap %i of %i\n',bootstrap,numBootstrap)
            if ar==1 % Input area names as a cell array to analyze all
                [RastFig, DecodeFig, ~, ~] = runNeuralDecodingToolbox(raster_data_dir,area_names,chCond,'allSpikes', ...
                    unitType,trialType,classifier_labels,classifierName,'num_cv_splits',splitNumber,...
                    'numUnits',run,'numBootstrap',bootstrap,'outputFolder',outputFolder,'shuffle',false,'trialCompare',trialCompare);
            else
                [RastFig, DecodeFig, ~, ~] = runNeuralDecodingToolbox(raster_data_dir,tempAreaName,chCond,'allSpikes', ...
                    unitType,trialType,classifier_labels,classifierName,'num_cv_splits',splitNumber,...
                    'numUnits',run,'numBootstrap',bootstrap,'outputFolder',outputFolder,'shuffle',false,'trialCompare',trialCompare);
            end
            close(RastFig), close(DecodeFig)
        end
    end
end

%% Helper functions
% Function executed by the "doneButton" in subpopulation refinement
function doneExe(~,~,f,ui_field,area_names)

for ar = 1:numel(area_names)
    switch ui_field(ar).Value
        case 2
            area_names{ar} = [area_names{ar} '-RS'];
        case 3
            area_names{ar} = [area_names{ar} '-FS'];
    end
end

% Overwrite variables
assignin('caller','area_names',area_names)

close(f)
end

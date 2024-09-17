function DE_Salience_function(dataDir, varargin)
% varargin as name-value pair, with timeLapse, responseWindow, binSz,
% spontaneousWindow, condition, onOffStr, getWaveFlag, signFilterFlag,
% overwriteFlag, PSTHorder, adjustConditions

% user inputs
% (1) timeLapse: Viewing window [s] - default: [-0.8, 0.4]
% (2) responseWindow: Response window [s] - default: [0, 0.2]
% (3) binSz: Bin size [s] - default: 0.01
% (4) spontaneousWindow: Spontaneous window [s] - default: [-0.8, -0.6]
% (5) condition: 'Reward', 'Punishment', 'Lick', 'WhiskerContact_left', 'WhiskerContact_right',
%       'WhiskerContact_onlyLeftFirst','WhiskerContact_onlyRightFirst','onlyFirstLick' - default: 'WhiskerContact_left'
% (6) onOffStr: Trigger (onset or offset): 'on', 'off' - default: 'on'
% (7) getWaveFlag: Get waveforms: all, sign, none - default: all
% (8) signFilterFlag: Significance filter: true, false - default: true
% (9) overwriteFlag: Overwrite existing csv files: true, false - default: true
% (10) PSTHorder: PSTH order: 'id', 'ActiveUnit', 'Amplitude', 'ContamPct', 'KSLabel', 'amp', 'ch', 'depth', 'fr', 'group', 'n_spikes', 'sh' - default: 'depth'
% (11) adjustConditions: Exclude some triggers, if you only want to analyze e.g. "Go Successes" - default: 'all'
% (12) interruptAfterPSTH: If you only want to retrieve significant responses set to true - default: false
% (13) timeShiftType: Only relevant for the experiment of Matthias-Martin: Time shift can either be 'static', 'dynamic' or 'none' - default: 'none'
% (14) timeShiftValue: Only relevant for the experiment of Matthias-Martin: Desired time shift of spike data in msec - default: 0

% Default option values:
options = inputParser;

checkonOffStr = @(x) any(validatestring(x,{'on','off'}));
checkgetWaveFlag = @(x) any(validatestring(x,{'none', 'sign', 'all'}));

addRequired(options,'dataDir',@ischar);
addParameter(options,'timeLapse',[-0.8, 0.4],@isnumeric)
addParameter(options,'responseWindow',[0, 0.2],@isnumeric)
addParameter(options,'binSz',0.01,@isnumeric)
addParameter(options,'spontaneousWindow',[-0.8, -0.6],@isnumeric)
addParameter(options,'condition','WhiskerContact_left',@ischar)
addParameter(options,'onOffStr','on',checkonOffStr)
addParameter(options,'getWaveFlag','all',checkgetWaveFlag)
addParameter(options,'signFilterFlag',true,@islogical)
addParameter(options,'overwriteFlag',true,@islogical)
addParameter(options,'PSTHorder',{'depth', 'responsiveness'},@iscell)
addParameter(options,'adjustConditions',{'all'},@iscell)
addParameter(options,'interruptAfterPSTH',false,@islogical)
addParameter(options,'timeShiftType','none',@ischar)
addParameter(options,'timeShiftValue',0,@isnumeric)

parse(options,dataDir,varargin{:})
timeLapse = options.Results.timeLapse;
responseWindow = options.Results.responseWindow;
binSz = options.Results.binSz;
spontaneousWindow = options.Results.spontaneousWindow;
condition = options.Results.condition;
onOffStr = options.Results.onOffStr;
getWaveFlag = options.Results.getWaveFlag;
signFilterFlag = options.Results.signFilterFlag;
overwriteFlag = options.Results.overwriteFlag;
PSTHorder = options.Results.PSTHorder;
adjustConditions = options.Results.adjustConditions;
interruptAfterPSTH = options.Results.interruptAfterPSTH;
timeShiftType = options.Results.timeShiftType;
timeShiftValue = options.Results.timeShiftValue;



% Make sure the time shift option is only used in combination with the
% early interruption to prevent overwriting of existing files

assert(ismember(timeShiftType,{'static','dynamic','none'}),...
    'Salience_function:timeShiftType',...
    'timeShiftType must be one of the following values: static, dynamic, none')

if timeShiftValue
    assert(all(interruptAfterPSTH, logical(timeShiftValue)),...
        'Salience_function:timeShiftAndInterruption',...
        'When a time shift is used, make sure the interruptAfterPSTH parameter is set to true.')
end

assert(round(diff(responseWindow),3)==round(diff(spontaneousWindow),3),...
    'Salience_function:WindowsNotSameLength',...
    'responseWindow and spontaneousWindow must be of the same length.')

assert(exist(dataDir,'dir'),...
    'Salience_function:NodataDir',...
    'Your data path does not exist.')

validPSTHorder = {'id', 'ActiveUnit',...
    'Amplitude', 'ContamPct', 'KSLabel', 'amp', 'ch',...
    'depth', 'fr', 'group', 'n_spikes', 'sh', 'responsiveness'};

assert(numel(PSTHorder) <= 2,...
    'Salience_function:PSTHorderInputNumber',...
    'PSTHorder can only have 2 or less inputs')

cellfun(@(x) assert(ismember(x,validPSTHorder),...
    'Salience_function:PSTHorderInputValues',...
    sprintf('PSTHorder must be one of the following values: %s', strjoin(validPSTHorder,', '))),...
    PSTHorder)

validAdjustConditions = {'all', 'wide', 'narrow', 'intermediate', 'lick', 'no-lick'};

assert(numel(adjustConditions) <= 2,...
    'Salience_function:adjustConditionsInputNumber',...
    'adjustConditions can only have 2 or less inputs')

cellfun(@(x) assert(ismember(x,validAdjustConditions),...
    'Salience_function:adjustConditionsInputValues',...
    sprintf('adjustConditions must be one of the following values: %s', strjoin(validAdjustConditions,', '))),...
    adjustConditions)

FileInfo = dir(fullfile(dataDir,sprintf('StimulusResponse_%s*',condition)));
if ~any(cellfun(@isempty, adjustConditions)) && ~isempty(FileInfo) && ~any(contains(adjustConditions,'all'))
    load(fullfile(FileInfo(1).folder,FileInfo(1).name),'Wide_Narrow_Intermediate','Lick')
    if any(ismember(adjustConditions,'wide'))
        logApert = Wide_Narrow_Intermediate==1;
    elseif any(ismember(adjustConditions,'narrow'))
        logApert = Wide_Narrow_Intermediate==2;
    elseif any(ismember(adjustConditions,'intermediate'))
        logApert = Wide_Narrow_Intermediate==3;
    else
        logApert = true(length(Wide_Narrow_Intermediate),1);
    end

    if any(ismember(adjustConditions,'lick'))
        logLick = Lick==1;
    elseif any(ismember(adjustConditions,'no-lick'))
        logLick = Lick==0;
    else
        logLick = true(length(Lick),1);
    end
    logicalFlag = logApert & logLick;
end
%% Load the data

% Creating the figure directory
figureDir = fullfile(dataDir,'Figures\');
if ~exist(figureDir,'dir')
    if ~mkdir(figureDir)
        fprintf(1,'There was an issue with the figure folder...\n');
        return
    end
end

tempDir = dataDir;
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
[~,expName,~] = fileparts(rhdFile(1).name);
intanDir = rhdFile(1).folder;

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

% Adjust Triggers accordingly
if ~any(cellfun(@isempty, adjustConditions)) && ~any(contains(adjustConditions,'all'))
    idx = find(cellfun(@(x) isequal(x,condition), {Conditions.name}));
    for i = 1:numel(idx)
        Conditions(idx(i)).Triggers = sort(Conditions(idx(i)).Triggers);
        Conditions(idx(i)).Triggers = Conditions(idx(i)).Triggers(logicalFlag,:);
    end
end

% Assign clusters to area and calculate mean z-scores
try
    dataInfo = dir(fullfile(dataDir,'*all_channels.mat'));
    load(fullfile(dataInfo(1).folder,dataInfo(1).name),'sortedData','fs')
catch
    try
        [sortedData, fs] = importPhyFiles(dataDir);
    catch
        fprintf(2,'Error importing the phy files into Matlab format\n')
        return
    end
end

% Sometimes empty clusters are registered, which are excluded here
sortedData = sortedData(all(~cellfun(@isempty, sortedData),2),:);

% Cave: only for Matthias-Martin experiment
if timeShiftValue && isequal(timeShiftType,'static')
    % Convert time shift from msec to sec and add to spike data
    sortedData(:,2) = cellfun(@(x) x+timeShiftValue/1000,sortedData(:,2),'UniformOutput',false);
elseif timeShiftValue && isequal(timeShiftType,'dynamic')
    % Convert time shift from msec to sec and for each second add the dynamic time shift
    sortedData(:,2) = cellfun(@(x) arrayfun(@(y) y+floor(y)*timeShiftValue/1000, x),sortedData(:,2),'UniformOutput',false);
end

clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));

clInfo.group(cell2mat(sortedData(:,3))==1) = deal({'good'});
clInfo.group(cell2mat(sortedData(:,3))==2) = deal({'mua'});
clInfo.group(cell2mat(sortedData(:,3))==3) = deal({'noise'});

try
    animalPath = fileparts(fileparts(fileparts(intanDir)));
    load(fullfile(animalPath,'targetHit.mat'),'targetHit')
catch
    fprintf(2,'\nNo targetHit.mat file.\n')
    return
end

include = ismember(clInfo.group,{'good','mua'}) & clInfo.isolation_distance > 15 & clInfo.isi_viol < 3 & ismember(clInfo.ch,targetHit);

% Add labels for brain regions, depending on recording depths
sortedData(:,4) = {'N/A'};
for i = 1:height(sortedData)
    depth = clInfo.depth(strcmp(clInfo.id,sortedData{i,1}));
    if ~isempty(depth)
        switch depth
            case 1
                subscript = 'BC';
            case 1400
                subscript = 'POm';
            case 1700
                subscript = 'VPM';
            case 2400
                subscript = 'ZIv';
        end
    else
        subscript = 'NaN';
    end
    sortedData{i,4} = subscript;
end

% Load Conditions and Triggers

% Map the spikes with synchronized syntalos timestamps
% intanFile = fullfile(dataDir,strcat(expName,'.rhd'));
% sortedData = synchronizeSpikes(sortedData,intanFile,fs);

%% Constructing the helper 'global' variables
% Sampling frequency of continuous traces (i.e., ADC data)
fsCont = 10000;

% Number of total samples
Ns = structfun(@numel,Triggers);
Ns = min(Ns(Ns>1))*fs/fsCont;
% Total duration of the recording
Nt = Ns/fs;
% Useless clusters (labeled as noise or they have very low firing rate)

%%% THIS IS JUST FOR TESTING
%%% REPLACES ALL BAD UNITS WITH MUA
%%% CHECK RESPONSIVENESS OF NOISE UNITS
% [sortedData{cellfun(@(x) x==3,sortedData(:,3)),3}] = deal(2);

% badsIdx = cellfun(@(x) x==3,sortedData(:,3));
% bads = find(badsIdx);
% totSpkCount = cellfun(@numel,sortedData(:,2));
% clusterSpikeRate = totSpkCount/Nt;
% silentUnits = clusterSpikeRate < 0.1;
% bads = union(bads,find(silentUnits));
% badsIdx = badsIdx | silentUnits;
%
% goodsIdx = setdiff(1:size(sortedData,1),bads);
goodsIdx = find(include);
act_units = zeros(size(clInfo,1),1);
act_units(ismember(clInfo.id(:),sortedData(goodsIdx,1))) = 1;
if ~any(ismember(clInfo.Properties.VariableNames,'ActiveUnit'))
    clInfo = addvars(clInfo,act_units,'After','id',...
        'NewVariableNames','ActiveUnit');
end
gclID = sortedData(goodsIdx,1);
% Logical spike trace for the first good cluster
spkLog = StepWaveform.subs2idx(round(sortedData{goodsIdx(1),2}*fs),Ns);
% Subscript column vectors for the rest good clusters
spkSubs = cellfun(@(x) round(x.*fs),sortedData(goodsIdx(2:end),2),...
    'UniformOutput',false);
% Number of good clusters
Ncl = numel(goodsIdx);
% Redefining the stimulus signals from the low amplitude to logical values
whStim = {'piezo','whisker','mech','audio'};
cxStim = {'laser','light','reward','punishment','lick'};
lfpRec = {'lfp','s1','cortex','s1lfp'};
trigNames = fieldnames(Triggers);
numTrigNames = numel(trigNames);
ctn = 1;
continuousSignals = cell(numTrigNames,1);
continuousNameSub = zeros(size(trigNames));
while ctn <= numTrigNames
    if contains(trigNames{ctn},whStim,'IgnoreCase',true)
        continuousSignals{ctn} = Triggers.(trigNames{ctn});
        continuousNameSub(ctn) = ctn;
    end
    if contains(trigNames{ctn},cxStim,'IgnoreCase',true)
        continuousSignals{ctn} = Triggers.(trigNames{ctn});
        continuousNameSub(ctn) = ctn;
    end
    if contains(trigNames{ctn},lfpRec,'IgnoreCase',true)
        continuousSignals{ctn} = Triggers.(trigNames{ctn});
        continuousNameSub(ctn) = ctn;
    end
    ctn = ctn + 1;
end
continuousSignals(continuousNameSub == 0) = [];
% continuousNameSub(continuousNameSub == 0) = [];
% trigNames = trigNames(continuousNameSub);
%% Inter-spike intervals
isiFile = fullfile(dataDir,[expName,'_ISIvars.mat']);
if ~exist(isiFile,'file')
    spkSubs2 = cellfun(@(x) round(x.*fs), sortedData(goodsIdx,2),...
        'UniformOutput', false);
    ISIVals = cellfun(@(x) [x(1)/fs; diff(x)/fs], spkSubs2,...
        'UniformOutput', 0);
    NnzvPcl = cellfun(@numel,ISIVals);
    %     Nnzv = sum(NnzvPcl);
    rows = cell2mat(arrayfun(@(x,y) repmat(x,y,1), (1:Ncl)', NnzvPcl,...
        'UniformOutput', 0));
    cols = cell2mat(spkSubs2);
    vals = cell2mat(ISIVals);
    try
        ISIspar = sparse(rows, cols, vals);
    catch
        fprintf(1, 'Not possible to create such a big array')
    end
else
    load(isiFile,'ISIspar')
end
% ISIsignal = zeros(Ncl,Ns,'single');
% for ccl = 1:Ncl
%     ISIsignal(ccl,spkSubs2{ccl}) = ISIVals{ccl};
% end
%% User controlling variables

fprintf(1,'\nTrial types to plot: %s\n',strjoin(adjustConditions, ' & '))
fprintf(1,'PSTH ordered according to: %s\n',strjoin(PSTHorder, ' & '))

if signFilterFlag
    fprintf(1,'Only significant units analyzed.\n')
else
    fprintf(1,'All units analyzed.\n')
end

fprintf(1,'Time window: %.2f - %.2f ms\n',timeLapse(1)*1e3, timeLapse(2)*1e3)
fprintf(1,'Response window: %.2f - %.2f ms\n',responseWindow(1)*1e3, responseWindow(2)*1e3)
fprintf(1,'Bin size: %.3f ms\n', binSz*1e3)
fprintf(1,'Spontaneous window: %.2f to %.2f ms before the trigger\n\n',...
    spontaneousWindow(1)*1e3, spontaneousWindow(2)*1e3)
%% Condition triggered stacks

condNames = arrayfun(@(x) x.name,Conditions,'UniformOutput',false);
% match chosen condition with condition list
match = find(cellfun(@(x) isequal(condition,x), condNames));
chCond = match(1);

%% Constructing the stack out of the user's choice
% discStack - dicrete stack has a logical nature
% cst - continuous stack has a numerical nature
% Both of these stacks have the same number of time samples and trigger
% points. They differ only in the number of considered events.

[discStack, cst] = getStacks(spkLog,Conditions(chCond).Triggers,onOffStr,...
    timeLapse,fs,fsCont,spkSubs,continuousSignals);
% ISI stack
try
    [~, isiStack] = getStacks(spkLog,Conditions(chCond).Triggers, onOffStr,...
        timeLapse,fs,fsCont,[],ISIspar);
catch
    fprintf(1,'Not able to do the ISI stack\n')
end
% [dst, cst] = getStacks(spkLog, allWhiskersPlusLaserControl,...
%     'on',timeLapse,fs,fs,[spkSubs;{Conditions(allLaserStimulus).Triggers}],...
%     continuousSignals);
if ~exist(isiFile,'file') && exist('ISIspar','var') && ~timeShiftValue
    fprintf(1,'Saving the inter-spike intervals for each cluster... ');
    save(isiFile,'ISIspar','ISIVals','-v7.3')
    fprintf(1,'Done!\n')
end
% Number of clusters + the piezo as the first event + the laser as the last
% event, number of time samples in between the time window, and number of
% total triggers.
[Ne, Nt, NTa] = size(discStack);
% Computing the time axis for the stack
tx = (0:Nt - 1)/fs + timeLapse(1);
%% Considered conditions selection
% Choose the conditions to look at
auxSubs = setdiff(1:numel(condNames), chCond);
ccondNames = condNames(auxSubs);

% match chosen condition with remaining condition list
match = cellfun(@(x) isequal(condition,x), ccondNames);
cchCond = find(match,1);

% Select the onset or the offset of a trigger
fprintf(1,'Condition(s):\n')
fprintf('- ''%s''\n', Conditions(auxSubs(cchCond)).name)

consideredConditions = auxSubs(cchCond);
Nccond = length(consideredConditions);

%% Boolean flags
delayFlags = false(NTa,Nccond);
counter2 = 1;
for ccond = consideredConditions
    delayFlags(:,counter2) = ismember(Conditions(chCond).Triggers(:,1),...
        Conditions(ccond).Triggers(:,1));
    counter2 = counter2 + 1;
end

%% Computing which units/clusters/putative neurons respond to the stimulus
% Logical indices for fetching the stack values
sponActStackIdx = tx >= spontaneousWindow(1) & tx <= spontaneousWindow(2);
respActStackIdx = tx >= responseWindow(1) & tx <= responseWindow(2);
% The spontaneous activity of all the clusters, which are allocated from
% the second until one before the last row, during the defined spontaneous
% time window, and the whisker control condition.

timeFlags = [sponActStackIdx;respActStackIdx];
% Time window
delta_t = diff(responseWindow);
% Statistical tests
[Results, Counts] = statTests(discStack, delayFlags, timeFlags);

indCondSubs = cumsum(Nccond:-1:1);
consCondNames = condNames(consideredConditions);
% Plotting statistical tests
[Figs, Results] = scatterSignificance(Results, Counts,...
    strrep(condNames(consideredConditions),'_','\_'), delta_t, sortedData(goodsIdx,1));
configureFigureToPDF(Figs);
stFigBasename = fullfile(figureDir,[expName,' ']);
stFigSubfix = sprintf(' Stat RW%.1f-%.1fms SW%.1f-%.1fms.fig',...
    responseWindow(1)*1e3, responseWindow(2)*1e3, spontaneousWindow(1)*1e3,...
    spontaneousWindow(2)*1e3);
ccn = 1;
%for cc = indCondSubs
for cc = 1:numel(Figs)
    if ~ismember(cc, indCondSubs)
        altCondNames = strsplit(Figs(cc).Children(2).Title.String,': ');
        altCondNames = altCondNames{2};
    else
        altCondNames = consCondNames{ccn};
        ccn = ccn + 1;
    end
    stFigName = [stFigBasename, altCondNames, stFigSubfix];
    if ~timeShiftValue % Save only when there is no time shift
        savefig(Figs(cc), stFigName)
    end
end
H = cell2mat(cellfun(@(x) x.Pvalues,...
    arrayfun(@(x) x.Activity, Results(indCondSubs), 'UniformOutput', 0),...
    'UniformOutput', 0)) < 0.05;

% Htc = sum(H,2);
CtrlCond = contains(consCondNames,'control','IgnoreCase',true);
if ~nnz(CtrlCond)
    CtrlCond = true(size(H,2),1);
end
wruIdx = any(H(:,CtrlCond),2);
Nwru = nnz(wruIdx);

fprintf('%d whisker responding clusters:\n', Nwru);
fprintf('- %s\n',gclID{wruIdx})

if interruptAfterPSTH
    % Interrupt script here, to retrieve significantly responding clusters
    assignin('base','goodClusters',gclID)
    assignin('base','significantResponse',wruIdx)
    assignin('base','clInfo',clInfo)
    return
end

%% Get waveforms of respective clusters
if isequal(getWaveFlag,'all')
    waveFileName = [expName,'_waveforms_all.mat'];
    waveFile = dir(fullfile(dataDir,waveFileName));
    if isempty(waveFile)
        getClusterWaveform(gclID, dataDir, waveFileName);
    else
        fprintf('\nWaveforms already analyzed.\n')
    end
elseif isequal(getWaveFlag,'sign')
    waveFileName = [expName,sprintf('_waveforms_sign_%s_%s_%.2fbz_%srespWind_%sspontWind_%s.mat',...
        altCondNames,onOffStr,binSz,[num2str(responseWindow(1)),'_',num2str(responseWindow(2))],...
        [num2str(spontaneousWindow(1)),'_',num2str(spontaneousWindow(2))],strjoin(adjustConditions,'&'))];
    waveFile = dir(fullfile(dataDir,waveFileName));
    if isempty(waveFile)
        getClusterWaveform(gclID(wruIdx), dataDir, waveFileName);
    else
        fprintf('\nWaveforms already analyzed.\n')
    end
end

%% Configuration structure
configStructure = struct('Experiment', fullfile(dataDir,expName),...
    'Viewing_window_s', timeLapse, 'Response_window_s', responseWindow,...
    'BinSize_s', binSz, 'Trigger', struct('Name', condNames{chCond},...
    'Edge',onOffStr), 'ConsideredConditions',{consCondNames});

%% Filter question
filterIdx = true(Ne,1);
filtStr = 'unfiltered';
if signFilterFlag
    filterIdx = [true; wruIdx];
    filtStr = 'filtered';
end

%% Getting the relative spike times for the whisker responsive units (wru)
% For each condition, the first spike of each wru will be used to compute
% the standard deviation of it.
cellLogicalIndexing = @(x,idx) x(idx);
isWithinResponsiveWindow =...
    @(x) x > responseWindow(1) & x < responseWindow(2);

firstSpike = zeros(Nwru,Nccond);
M = 16;
binAx = responseWindow(1):binSz:responseWindow(2);
condHist = zeros(size(binAx,2)-1, Nccond);
firstOrdStats = zeros(2,Nccond);
condParams = zeros(M,3,Nccond);
txpdf = responseWindow(1):1/fs:responseWindow(2);
condPDF = zeros(numel(txpdf),Nccond);
csvBase = fullfile(dataDir, expName);
csvSubfx = sprintf(' VW%.1f-%.1f ms.csv', timeLapse(1)*1e3, timeLapse(2)*1e3);
existFlag = false;
condRelativeSpkTms = cell(Nccond,1);
relativeSpkTmsStruct = struct('name',{},'SpikeTimes',{});
spkDir = fullfile(dataDir, 'SpikeTimes');
for ccond = 1:size(delayFlags,2)
    csvFileName = [csvBase,' ',consCondNames{ccond}, csvSubfx];
    relativeSpikeTimes = getRasterFromStack(discStack,~delayFlags(:,ccond),...
        filterIdx(3:end), timeLapse, fs, true, false);
    relativeSpikeTimes(:,~delayFlags(:,ccond)) = [];
    relativeSpikeTimes(~filterIdx(2),:) = [];
    condRelativeSpkTms{ccond} = relativeSpikeTimes;
    %     respIdx = cellfun(isWithinResponsiveWindow, relativeSpikeTimes,...
    %         'UniformOutput',false);
    clSpkTms = cell(size(relativeSpikeTimes,1),1);
    if exist(csvFileName, 'file') && ccond == 1
        existFlag = true;
        if overwriteFlag
            existFlag = false;
            fprintf(1,'Overwriting... ');
        end
    end
    fID = 1;
    if ~existFlag
        fID = fopen(csvFileName,'w');
        fprintf(fID,'%s, %s\n','Cluster ID','Relative spike times [ms]');
    end
    for cr = 1:size(relativeSpikeTimes, 1)
        clSpkTms(cr) = {sort(cell2mat(relativeSpikeTimes(cr,:)))};
        if fID > 2
            fprintf(fID,'%s,',gclID{cr});
            fprintf(fID,'%f,',clSpkTms{cr});fprintf(fID,'\n');
        end
    end
    if fID > 2
        fclose(fID);
    end
    relativeSpkTmsStruct(ccond).name = consCondNames{ccond};
    relativeSpkTmsStruct(ccond).SpikeTimes = condRelativeSpkTms{ccond};
    %{
    % First spike in the trial
    spikeTimesINRespWin = cellfun(cellLogicalIndexing,...
        relativeSpikeTimes, respIdx, 'UniformOutput',false);
    allSpikeTimes = cell2mat(spikeTimesINRespWin(:)');
    condParams(:,:,ccond) = emforgmm(allSpikeTimes, M, 1e-6, 0);
    condPDF(:,ccond) = genP_x(condParams(:,:,ccond), txpdf);
    firstOrdStats(:,ccond) = [mean(allSpikeTimes), std(allSpikeTimes)];
    hfig = figure('Visible', 'off'); h = histogram(allSpikeTimes, binAx,...
        'Normalization', 'probability');
    condHist(:,ccond) = h.Values;
    close(hfig)
    for ccl = 1:Nwru
        frstSpikeFlag = ~cellfun(@isempty,spikeTimesINRespWin(ccl,:));
        firstSpike(ccl,ccond) = std(...
            cell2mat(spikeTimesINRespWin(ccl,frstSpikeFlag)));
    end
    %}
end

relSpkFileName =...
    sprintf('%s RW%.2f - %.2f ms SW%.2f - %.2f ms %s exportSpkTms.mat',...
    expName, responseWindow*1e3, spontaneousWindow*1e3,...
    Conditions(chCond).name);
if ~exist(relSpkFileName,'file')
    save(fullfile(dataDir, relSpkFileName), 'relativeSpkTmsStruct',...
        'configStructure')
end
%% Standard Deviations of First Spikes After Each Trigger per Unit
% firstSpikes(relativeSpkTmsStruct, gclID, dataDir);

%% Plot PSTH

% Generating orderedStr
plotID = gclID(filterIdx(2:Ncl+1));
signID = gclID(wruIdx);
orderedStr = [strjoin(PSTHorder,' & '), ' ordered'];

% goodsIdx = logical(clInfo.ActiveUnit);
csNames = fieldnames(Triggers);
Nbn = diff(timeLapse)/binSz;
Nbn = round(Nbn);
% if (round(Nbn,3) - round(Nbn)) ~= 0
%     Nbn = ceil(Nbn);
% else
%     Nbn = round(Nbn);
% end
PSTH = zeros(nnz(filterIdx) - 1, Nbn, Nccond);
PSTH_all = zeros(Ncl, Nbn, Nccond);
PSTHn_all = zeros(Ncl, Nbn, Nccond);
psthTx = (0:Nbn-1) * binSz + timeLapse(1);
spontaneousBins = psthTx>=spontaneousWindow(1) & psthTx<spontaneousWindow(2);
responseBins = psthTx>=responseWindow(1) & psthTx<responseWindow(2);

psthFigs = gobjects(Nccond,1);
Ntc = size(cst,2);
if ~isempty(plotID)
    for ccond = 1:Nccond
        figFileName =...
            sprintf('%s %s VW%.1f-%.1f ms B%.1f ms RW%.1f-%.1f ms SW%.1f-%.1f ms %sset %s (%s).fig',...
            expName, consCondNames{ccond}, timeLapse*1e3, binSz*1e3,...
            responseWindow*1e3, spontaneousWindow*1e3, onOffStr, orderedStr,...
            filtStr);
        [PSTH(:,:,ccond), trig, sweeps] = getPSTH(discStack(filterIdx,:,:),timeLapse,...
            ~delayFlags(:,ccond),binSz,fs);

        [PSTH_all(:,:,ccond), ~, ~] = getPSTH(discStack,timeLapse,...
            ~delayFlags(:,ccond),binSz,fs);

        % Normalizing the PSTH to the maximal value on each cluster.
        %         PSTHn = (PSTH-min(PSTH,[],2))./max(PSTH,[],2);
        PSTHn = PSTH(:,:,ccond)./max(PSTH(:,:,ccond),[],2);
        PSTHn_all(:,:,ccond) = PSTH_all(:,:,ccond)./max(PSTH_all(:,:,ccond),[],2);

        responsiveness = nan(height(PSTH(:,:,ccond)),1);
        for unit = 1:height(PSTH(:,:,ccond))
            responsiveness(unit) = mean(PSTHn(unit,responseBins,ccond))-mean(PSTHn(unit,spontaneousBins,ccond));
        end

        if ~exist('clInfo','var')
            clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
        end
        clInfo = [clInfo, table(nan(height(clInfo),1),'VariableNames',{'responsiveness'})];
        clInfo{plotID,'responsiveness'} = responsiveness;

        [ordSubs, names] = getPSTHorder(PSTHorder,clInfo,plotID,filterIdx,Ncl,sortedData);

        if exist('cst', 'var') && ~isempty(cst)
            stims = mean(cst(:,:,delayFlags(:,ccond)),3);
            stims = stims - median(stims,2);
            for cs = 1:size(stims,1)
                if abs(log10(var(stims(cs,:),[],2))) < 13
                    [m,b] = lineariz(stims(cs,:),1,0);
                    stims(cs,:) = m*stims(cs,:) + b;
                else
                    stims(cs,:) = zeros(1,Ntc);
                end
            end
        else
            stims = zeros(1, Ntc);
        end
        psthFigs(ccond) = plotClusterReactivity(PSTH(ordSubs,:,ccond), trig,...
            sweeps, timeLapse, binSz, [consCondNames(ccond); names],...
            strrep(expName,'_','\_'), stims, csNames);
        psthFigs(ccond).Children(end).YLabel.String =...
            [psthFigs(ccond).Children(end).YLabel.String,...
            sprintf('^{%s}',orderedStr)];
        figFilePath = fullfile(figureDir, figFileName);
        savefig(psthFigs(ccond), figFilePath);
    end

    % Save the table for the significantly responding units
    signResp = cellfun(@(x) ismember(x, signID), sortedData(:,1));
    signTable = cell2table([sortedData(:,[1,3,4]),num2cell(include),num2cell(signResp)],"VariableNames",{'unitID','group','area','include','signResponse'});

    save(fullfile(dataDir,sprintf('allResponses_%s_%s_%.2fbz_%stimeLapse_%srespWind_%sspontWind_%s.mat',...
        altCondNames,onOffStr,binSz,[num2str(timeLapse(1)),'_',num2str(timeLapse(2))],...
        [num2str(responseWindow(1)),'_',num2str(responseWindow(2))],...
        [num2str(spontaneousWindow(1)),'_',num2str(spontaneousWindow(2))],strjoin(adjustConditions,'&'))),'signTable','PSTH_all','timeLapse','responseWindow','spontaneousWindow')
else
    % If there are no significantly responding units, create zero vector
    signResp = cellfun(@(x) ismember(x, signID), sortedData(:,1));
    signTable = cell2table([sortedData(:,[1,3,4]),num2cell(include),num2cell(signResp)],"VariableNames",{'unitID','group','area','include','signResponse'});

    save(fullfile(dataDir,sprintf('allResponses_%s_%s_%.2fbz_%stimeLapse_%srespWind_%sspontWind_%s.mat',...
        altCondNames,onOffStr,binSz,[num2str(timeLapse(1)),'_',num2str(timeLapse(2))],...
        [num2str(responseWindow(1)),'_',num2str(responseWindow(2))],...
        [num2str(spontaneousWindow(1)),'_',num2str(spontaneousWindow(2))],strjoin(adjustConditions,'&'))),'signTable','PSTH_all','timeLapse','responseWindow','spontaneousWindow')
end

end

function [ordSubs, names] = getPSTHorder(PSTHorder,clInfo,pclID,filterIdx,Ncl,sortedData)

ordSubs = 1:nnz(filterIdx(2:Ncl+1));
names = pclID(ordSubs);
if ~any(strcmp(PSTHorder,'id'))
    [~,ordSubs] = sortrows(clInfo(pclID,:),PSTHorder);

    if ismember(PSTHorder(1),{'depth','responsiveness'})
        names = cell(numel(pclID),1);
        for i = 1:numel(pclID)
            subscript = sortedData{strcmp(sortedData(:,1),char(pclID{ordSubs(i)})),4};
            names{i} = [pclID{ordSubs(i)},sprintf('_{%s}',subscript)];
        end
    elseif strcmp(PSTHorder(1),'group')
        names = cell(numel(pclID),1);
        for i = 1:numel(pclID)
            group_name = clInfo.group{strcmp(clInfo.id,char(pclID{ordSubs(i)}))};
            names{i} = [pclID{ordSubs(i)},sprintf('_{%s}',group_name)];
        end
    else
        names = pclID(ordSubs);
    end

end

end
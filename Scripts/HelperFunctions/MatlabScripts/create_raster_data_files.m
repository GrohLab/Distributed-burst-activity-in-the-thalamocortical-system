%% Create raster_data.mat files for The Neural Decoding Toolbox
% For each unit the code depends on a particular file, with three
% variables: raster_data, raster_labels, and raster_site_info
% raster_data: matrix as [num_trials x num_time_points]
% time points could be bins of milliseconds with either a spike (= 1) or no
% spike (= 0)
% raster_labels: structure that contains the different experimental conditions (e.g., Go, No-Go, Go-Success)
% raster_site_info: additional information (can also be empty)

function create_raster_data_files(fileSelection, outputFolder, chCond, trialType)

if exist(outputFolder,'dir') == 0
    mkdir(outputFolder)
end

for file = 1:height(fileSelection)
    filePath = fileSelection{file};
    fprintf('\nConverting file %d/%d\n',file,height(fileSelection))
    run_create_raster_data(filePath,outputFolder,chCond,trialType)
end

fprintf('\nDone!\n\n')
end

%% Create raster data
function run_create_raster_data(filePath,outputFolder,chCond,trialType)
animalName = regexp(filePath,'#\d*','once','match');
cohort_num = getCohort(filePath);
cohort_num = str2double(regexp(cohort_num,'\d*','match'));

stage_num = regexp(filePath,'P3.[0-9.]*','match','once');
stage_num = str2double(stage_num(4:end));

stage_str = regexp(filePath,'P3.[0-9.]*_[0-9.a-zA-Z_]*','match','once');
stage_str = regexp(stage_str,'_.*','match','once');
stage_str = stage_str(2:end);

ses_num = regexp(filePath,'session\d*','match','once');
ses_num = str2double(regexp(ses_num,'\d*','match','once'));

tempDir = filePath;
rhdFile = dir(fullfile(tempDir,'*analysis.mat'));
iter = 1;
while isempty(rhdFile) && iter <=3
    tempDir = fileparts(tempDir);
    rhdFile = dir(fullfile(tempDir,'*analysis.mat'));
    iter = iter + 1;
    if isempty(rhdFile) && iter > 3
        error('No rhd file found.')
    end
end
intanDir = rhdFile(1).folder;

FileInfo = dir(fullfile(filePath,sprintf('StimulusResponse_%s*.mat', chCond)));
load(fullfile(FileInfo(1).folder,FileInfo(1).name),'Go_NoGo_Neutral','Lick','Wide_Narrow_Intermediate')

% Define bin size in msec
binSize = 1;

% Define window before and after trigger onset [in msec]
before = 800;
after = 400;

% Extract trigger points and create raster_data for each cell and trial
% Get the Triggers
% Set warning to error in order to catch it
s = warning('error', 'MATLAB:load:variableNotFound');
condInfo = dir(fullfile(intanDir,'*analysis.mat'));
if ~isempty(condInfo)
    try load(fullfile(condInfo.folder,condInfo.name),'Conditions','update')
        if update < 5
            [Conditions, ~] = getConditions(fileparts(intanDir));
        end
    catch
        [Conditions, ~] = getConditions(fileparts(intanDir));
    end
else
    [Conditions, ~] = getConditions(fileparts(intanDir));
end
% Reset warning
warning(s);

try
    dataInfo = dir(fullfile(filePath,'*all_channels.mat'));
    load(fullfile(dataInfo(1).folder,dataInfo(1).name),'sortedData','fs')
    clInfo = getClusterInfo(fullfile(filePath,'cluster_info.tsv'));
catch
    try
        [sortedData, fs] = importPhyFiles(filePath);
    catch
        fprintf(2,'Error importing the phy files into Matlab format\n')
        return
    end
end
% Sometimes empty clusters are registered, which are excluded here
sortedData = sortedData(all(~cellfun(@isempty, sortedData),2),:);

str_idx = regexp(filePath,'#\d*','end');
animalPath = filePath(1:str_idx);
% Include only those units, which tetrodes were within the recording area
try
    load(fullfile(animalPath,'targetHit.mat'),'targetHit')
catch
    error(2,'\nNo targetHit.mat file. Progress with unfiltered analysis...\n')
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

% Check if the clusters from the BurstinessData file match those from the sortedData
load(fullfile(filePath,'BurstinessData.mat'),'TypeOfSpike','unitIDs')
if ~isequal(sortedData(:,1),unitIDs)
    error('burstsUponTrigger:unitMismatch',...
        'Units IDs of analyzed BurstinessData.mat file and sortedData do not match.')
end

% Load triggers in seconds
chCond_trig = sort(Conditions(find(cellfun(@(x) isequal(x,chCond), {Conditions.name}),1)).Triggers)./fs; % triggers in sec
if isequal(chCond,'onlyFirstLick')
    load(fullfile(fileparts(intanDir),'videos\HispeedTrials.mat'),'HispeedTrials')
    % This returns the first lick event after every middle point
    idx = unique(sort(cell2mat(arrayfun(@(x) find(chCond_trig(:,1)*1000 > x, 1), HispeedTrials.PreviousMP, 'UniformOutput', false))));
    chCond_trig = chCond_trig(idx,:);
end

% Create raster_labels with trial types
switch trialType
    case 'allTrials'
        raster_labels.Go_NoGo_Neutral = Go_NoGo_Neutral;
        raster_labels.Lick = Lick;
        raster_labels.Wide_Narrow_Intermediate = Wide_Narrow_Intermediate;
    case 'onlyGo'
        chCond_trig = chCond_trig(Go_NoGo_Neutral==1,:);
        raster_labels.Go_NoGo_Neutral = Go_NoGo_Neutral(Go_NoGo_Neutral==1);
        raster_labels.Lick = Lick(Go_NoGo_Neutral==1);
        raster_labels.Wide_Narrow_Intermediate = Wide_Narrow_Intermediate(Go_NoGo_Neutral==1);
    case 'onlyNoGo'
        chCond_trig = chCond_trig(Go_NoGo_Neutral==2,:);
        raster_labels.Go_NoGo_Neutral = Go_NoGo_Neutral(Go_NoGo_Neutral==2);
        raster_labels.Lick = Lick(Go_NoGo_Neutral==2);
        raster_labels.Wide_Narrow_Intermediate = Wide_Narrow_Intermediate(Go_NoGo_Neutral==2);
    case 'onlyNeutral'
        chCond_trig = chCond_trig(Go_NoGo_Neutral==3,:);
        raster_labels.Go_NoGo_Neutral = Go_NoGo_Neutral(Go_NoGo_Neutral==3);
        raster_labels.Lick = Lick(Go_NoGo_Neutral==3);
        raster_labels.Wide_Narrow_Intermediate = Wide_Narrow_Intermediate(Go_NoGo_Neutral==3);
    case 'onlyWide'
        chCond_trig = chCond_trig(Wide_Narrow_Intermediate==1,:);
        raster_labels.Go_NoGo_Neutral = Go_NoGo_Neutral(Wide_Narrow_Intermediate==1);
        raster_labels.Lick = Lick(Wide_Narrow_Intermediate==1);
        raster_labels.Wide_Narrow_Intermediate = Wide_Narrow_Intermediate(Wide_Narrow_Intermediate==1);
    case 'onlyNarrow'
        chCond_trig = chCond_trig(Wide_Narrow_Intermediate==2,:);
        raster_labels.Go_NoGo_Neutral = Go_NoGo_Neutral(Wide_Narrow_Intermediate==2);
        raster_labels.Lick = Lick(Wide_Narrow_Intermediate==2);
        raster_labels.Wide_Narrow_Intermediate = Wide_Narrow_Intermediate(Wide_Narrow_Intermediate==2);
    case 'onlyIntermediate'
        chCond_trig = chCond_trig(Wide_Narrow_Intermediate==3,:);
        raster_labels.Go_NoGo_Neutral = Go_NoGo_Neutral(Wide_Narrow_Intermediate==3);
        raster_labels.Lick = Lick(Wide_Narrow_Intermediate==3);
        raster_labels.Wide_Narrow_Intermediate = Wide_Narrow_Intermediate(Wide_Narrow_Intermediate==3);
    case 'onlyLick'
        chCond_trig = chCond_trig(Lick==1,:);
        raster_labels.Go_NoGo_Neutral = Go_NoGo_Neutral(Lick==1);
        raster_labels.Lick = Lick(Lick==1);
        raster_labels.Wide_Narrow_Intermediate = Wide_Narrow_Intermediate(Lick==1);
    case 'onlyNoLick'
        chCond_trig = chCond_trig(Lick==0,:);
        raster_labels.Go_NoGo_Neutral = Go_NoGo_Neutral(Lick==0);
        raster_labels.Lick = Lick(Lick==0);
        raster_labels.Wide_Narrow_Intermediate = Wide_Narrow_Intermediate(Lick==0);
end

WaveformInfo = dir(fullfile(filePath, '*waveforms_all.mat'));
try load(fullfile(WaveformInfo.folder,WaveformInfo.name),'clWaveforms')
catch
end
if all(cellfun(@(x,y) isequal(x,y), sortedData(:,1),clWaveforms(:,1)))
    [~, put_ex_idx] = getWaveformDistribution(clWaveforms,sortedData(:,4));
elseif all(cellfun(@(x) any(ismember(sortedData(:,1),x)), clWaveforms(:,1)))
    % Order the units according to the ClusterBurstiness table
    unitOrder = nan(height(clWaveforms),1);
    for i = 1:height(sortedData)
        unitOrder(i) = find(ismember(clWaveforms(:,1),sortedData{i,1}));
    end
    clWaveforms = clWaveforms(unitOrder,:);

    [~, put_ex_idx] = getWaveformDistribution(clWaveforms,sortedData(:,4));
else
    error('burstUponTrigger:waveformUnitsError','Units in the clWaveforms variable do not match the ones from the sortedData.')
end

for cl = 1:height(sortedData)

    spikeTypeNames = {'allSpikes','burstSpikes','tonicSpikes'};
    for spikeType = 1:3
        unitId = sortedData{cl,1};

        if put_ex_idx(cl)
            unitArea = [sortedData{cl,4},'-RS'];
        else
            unitArea = [sortedData{cl,4},'-FS'];
        end

        fileName = fullfile(outputFolder,sprintf('%s_%s_%s_%s_%s_stage%d_ses%d_cl%s_raster_data.mat',chCond,spikeTypeNames{spikeType},unitArea,trialType,animalName,stage_num,ses_num,unitId));
        if ~exist(fileName,'file')
            raster_data = zeros(height(chCond_trig),(before+after)/binSize);

            switch spikeType
                case 1
                    cl_spikeTimes = sortedData{cl,2};
                case 2
                    cl_spikeTimes = sortedData{cl,2}(logical(TypeOfSpike{cl}));
                case 3
                    cl_spikeTimes = sortedData{cl,2}(~logical(TypeOfSpike{cl}));
            end

            for trial = 1:height(chCond_trig)
                adjustSpikeTimes = round((cl_spikeTimes-chCond_trig(trial,1))*1000)+before;
                spikePoints = adjustSpikeTimes(adjustSpikeTimes > 0 & adjustSpikeTimes <= (before+after));
                raster_data(trial,spikePoints) = 1;
            end

            raster_site_info.unitId = unitId;
            raster_site_info.unitArea = unitArea;

            raster_site_info.alignment_event_time = before+1;
            raster_site_info.animalName = animalName;
            raster_site_info.cohort_num = cohort_num;
            raster_site_info.stage_num = stage_num;
            raster_site_info.stage_str = stage_str;
            raster_site_info.ses_num = ses_num;
            raster_site_info.spikeType = spikeTypeNames{spikeType};

            % Save data in output folder
            % Include area names and condition in file name, in order to be able to
            % decode them seperately later on
            save(fileName,'raster_data','raster_labels','raster_site_info')
        end
    end
end
end

%% Helper functions
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

% Choose trough to peak cutoff of 300 µs for thalamic nuclei
idx = find(cellfun(@(x) ismember(x, {'VPM','POm'}),cellArea));
put_in_idx(idx) = trough2peak_usec(idx) < 300;
put_ex_idx(idx) = trough2peak_usec(idx) >= 300;

% Choose trough to peak cutoff of 350 µs for cortex and ZI
idx = find(cellfun(@(x) ismember(x, {'BC','ZIv'}),cellArea));
put_in_idx(idx) = trough2peak_usec(idx) < 350;
put_ex_idx(idx) = trough2peak_usec(idx) >= 350;
end
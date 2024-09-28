function clWaveforms = getClusterWaveform(clusterID, dataDir, baseName)
%GETCLUSTERWAVEFORM reads the raw binary file and compiles the voltage
%traces for the given cluster (the cluster number should be the one
%assigned by Kilosort/Phy). The output is a cell (or a structure)
%containing the mean waveform and its standard deviation.
%   waveform = getClusterWaveform(clusterID)
%       INPUTS
%           - clusterID - character array, double or cell array containing
%           the clusters from which the waveforms are required.
%       OUTPUT
%           - waveform - cell array or double vector containing the mean
%           wavefrom for the given cluster(s).
% Emilio Isaias-Camacho @GrohLab 2019 (modified by Filippo Heimburg)

%% Input validation
clWaveforms = cell(1,3);
checkNature = @(x) [iscell(x), ischar(x), isnumeric(x)];
getLastCell = @(x) x{numel(x)};
if ~any(checkNature(clusterID))
    fprintf(1,'Unsupported input!\n')
    fprintf(1,'Be sure you input either the cluster ID as a ')
    fprintf(1,'character vector, a number or\nas a cell array')
    return
end
if ~exist(dataDir, 'dir')
    fprintf(1,'Not possible to retrieve waveforms without the data!\n')
    fprintf(1,'Please provide the data directory\n')
    return
end
% Converting the ID(s) to cell arrays according to their nature
% switch bit2int(checkNature(clusterID),'left-msb')
switch bi2de(checkNature(clusterID),'left-msb')
    case 1
        fprintf(1,'Numeric ID detected\n')
        clusterID = {num2str(clusterID)};
    case 2
        fprintf(1,'Character ID detected\n')
        clusterID = {clusterID};
    case 4
        fprintf(1,'Cell ID detected\n')
        charFlag = all(cellfun(@ischar, clusterID));
        if ~charFlag
            fprintf(1,'The cluster ID(s) should be only character within')
            fprintf(1,' a cell array\n')
            fprintf(1,'Please provide the ID(s) as required and try again\n')
            return
        end
end
clusterID = unique(clusterID);
%% Getting ready for the file reading
% Reading the cluster summary
clTable = getClusterInfo(fullfile(dataDir, 'cluster_info.tsv'));
% % Reading the channel map
fP = fopen(fullfile(dataDir,'params.py'),'r');
fgetl(fP);
ln = fgetl(fP);
fclose(fP);
Nch = double(getLastCell(textscan(ln,'%s = %d'))-1);
necessaryFiles = {fullfile(dataDir, 'channel_map.npy');...
    fullfile(dataDir, 'spike_templates.npy');...
    fullfile(dataDir, 'spike_clusters.npy')};
allOk = cellfun(@(x) exist(x,'file'), necessaryFiles);
if ~all(allOk)
    fprintf(1,'The following files were not found:\n')
    for cf = find(~allOk)
        fprintf(1,'%s\n', necessaryFiles{cf})
    end
    fprintf(1,'Cannot continue without these files. Aborting...\n')
    return
end
%% Reading the necessary files
% Reading the channel order
chanMap = readNPY(necessaryFiles{1});
% Preparatory variables for organising the output
spkTmls = readNPY(necessaryFiles{2});
spkCls = readNPY(necessaryFiles{3});
%% Input arguments verification
% Logical variables for clusters (clIdx) and spikes (spkIdx)
clIdx = false(size(clTable, 1), numel(clusterID));
spkIdx = false(size(spkCls,1), numel(clusterID));
clTempSubs = cell(numel(clusterID),1);
% Verifying if the given clusters exist in this experiment
for ccl = 1:numel(clusterID)
    clIdx(:,ccl) = strcmp(clTable.id, clusterID(ccl));
    spkIdx(:,ccl) = spkCls == str2double(clusterID(ccl));
    % Determining the template for the given cluster
    clTempSubs{ccl} = mode(spkTmls(spkIdx(:,ccl)));
end
missClustFlag = ~any(clIdx,1);
if ~all(~missClustFlag)    
    fprintf(1,'Some of the given clusters do not exist in this experiment\n')
    fprintf(1,'Clusters not found:\n')
    fprintf(1,'%s\n', clusterID{missClustFlag})
    if sum(missClustFlag) < numel(clusterID)
        contAns = questdlg('Continue without these clusters?', 'Continue?',...
            'Yes','No','Yes');
        if strcmp(contAns, 'No')
            fprintf(1,'Aborting...\n')
            return
        end
    else
        fprintf(1,'No valid cluster ID provided!\n')
        return
    end
end    
clusterID(missClustFlag) = [];
clIdx(:,missClustFlag) = [];
spkIdx(:, missClustFlag) = [];
clTempSubs(missClustFlag) = [];
[clSub,~] = find(clIdx);
% Determining hosting channels
% ch2read = chanMap(clTable{clusterID, 'channel'} + 1);
chanInTable = cellfun(@(x) ismember(x,{'ch','channels'}), clTable.Properties.VariableNames);
ch2read = clTable{clusterID, clTable.Properties.VariableNames{chanInTable}};

%% Verifying if the waveform(s) for the given cluster(s) was/were computed

waveFile = dir(fullfile(dataDir,baseName));

if ~isempty(waveFile)
    waveFileName = fullfile(dataDir, waveFile.name);
    load(waveFileName,'clWaveforms')
    % N_exCl = size(clWaveforms, 1);
    [exIdx, exSub] = ismember(clWaveforms(:,1),clusterID);
    % 'Removing' unwanted cluster waveforms and retrieving missing clusters
    clWaveforms = clWaveforms(exIdx,:);
    % Number of stored clusters
    Nstcl = nnz(exIdx);
    % Number of requested clusters
    Nrqcl = length(clusterID);
    missingClSubs = setdiff(1:Nrqcl, exSub);
    Nbncl = numel(missingClSubs);
    if Nbncl
        % Allocation of the bigger space
        clOutput = cell(Nstcl+Nbncl,3);
        % Reading the clusters from the file
        cl_fromBin = fetchWaveforms_fromBin(dataDir, clusterID(missingClSubs),...
            chanMap, Nch, clSub(missingClSubs), clTempSubs(missingClSubs),...
            spkIdx(:,missingClSubs), ch2read(missingClSubs));
        % Ordering the clusters alfabetically due to the string nature of the
        % ID
        reqClID = [clWaveforms(:,1); cl_fromBin(:,1)];
        [~, reqOrdSub] = sort(reqClID);
        % Assigning the clusters from the cell-array
        clOutput(reqOrdSub <= Nstcl,:) = clWaveforms;
        % Assigning the clusters from the bin file
        clOutput(reqOrdSub > Nstcl,:) = cl_fromBin;
        % Re-assigning the variable names
        clWaveforms = clOutput;
        save(waveFileName, 'clWaveforms','-append')
    end
else
    clWaveforms = fetchWaveforms_fromBin(dataDir, clusterID,...
        chanMap, Nch, clSub, clTempSubs, spkIdx, ch2read);
    waveFileName = fullfile(dataDir, baseName);
    save(waveFileName, 'clWaveforms')
end
end

function clWaveforms =...
    fetchWaveforms_fromBin(dataDir, clusterID,...
    chanMap, Nch, clSub, clTempSubs, spkIdx, ch2read)
%% Reading the binary file
clWaveforms = cell(numel(clusterID),3);
% Determinig the spike times for the given clusters
spikeFile = dir(fullfile(dataDir,'*_all_channels.mat'));
if ~isempty(spikeFile)
    load(fullfile(dataDir, spikeFile.name), 'sortedData', 'fs')
    if ~exist('fs','var')
        fsFile = dir(fullfile(dataDir,'*_sampling_frequency.mat'));
        load(fullfile(dataDir, fsFile.name), 'fs')
    end
end

% Sometimes empty clusters are registered, which are excluded here
sortedData = sortedData(all(~cellfun(@isempty, sortedData),2),:);

% Taking ~1.25 ms around the spike.
spikeWaveTime = 2*round(1.25e-3 * fs) + 1;
spikeSamples = (spikeWaveTime - 1)/2;

tempDir = dataDir;
binFile = dir(fullfile(tempDir,'*.bin'));
iter = 1;
while isempty(binFile) && iter <=3
    tempDir = fileparts(tempDir);
    binFile = dir(fullfile(tempDir,'*.bin'));
    iter = iter + 1;
    if isempty(binFile) && iter > 3
        fprintf(1, 'No binary file found.\n')
        fprintf(1, 'Without a binary file it is impossible to get the waveforms.\n')
        return
    end
end

pcFeat = readNPY(fullfile(dataDir, 'pc_features.npy'));
pcInd = readNPY(fullfile(dataDir, 'pc_feature_ind.npy'));
spkSubs = cellfun(@(x) round(x.*fs),sortedData(clSub,2),...
    'UniformOutput',false);
% [ch2read, readOrder, repeatChs] = unique(ch2read);
fID = fopen(fullfile(binFile.folder, binFile.name), 'r');
cchan = 1;
% Main loop
while ~feof(fID) && cchan <= numel(clusterID)
    % Computing the location of the channel features
    pcIdx = ch2read(cchan) == chanMap(pcInd(clTempSubs{cchan}+1,:)+1);
    clFeat = pcFeat(spkIdx(:,cchan), :, pcIdx);
    fprintf(1,'Reading channel %d ',ch2read(cchan))
    % Jumping to the channel
    fseek(fID, 2*(ch2read(cchan)), 'bof');
    % Computing the distance from spike to spike
    spkDists = [spkSubs{cchan}(1);diff(spkSubs{cchan})];
    fprintf(1,'looking for cluster %s...', clusterID{cchan})
    % Allocating space for the spikes
    waveform = zeros(spikeWaveTime, numel(spkSubs{cchan}));
    %fig = figure('Color',[1,1,1],'Visible', 'off');
    %ax = axes('Parent', fig); ax.NextPlot = 'add';
    %subSet = 1:floor(numel(spkDists)*0.1);
    for cspk = 1:numel(spkSubs{cchan})
        % Jumping to 1 ms before the time when the spike occured
        fseek(fID, 2*((Nch+1)*(spkDists(cspk) - spikeSamples)), 'cof');
        % Reading the waveform
        cwf = fread(fID, [spikeWaveTime, 1], 'int16=>single', 2*Nch);
        % Assigning the waveform to the saving variable. If the waveform is
        % cut, then assign only the gathered piece.
        try
            waveform(:,cspk) = cwf;
        catch
            waveform(1:length(cwf),cspk) = cwf;
        end
        % Jumping back to the exact time of the spike
        fseek(fID, -2*((Nch+1)*(spikeSamples+1)), 'cof');
        %    if ismember(cspk,subSet)
        %        plot(ax,waveform(:,cspk),'DisplayName',num2str(cspk));
        %    end
    end
    fprintf(1,' done!\n')
    clWaveforms(cchan,:) = [clusterID(cchan), {waveform}, {clFeat}];
    cchan = cchan + 1;
    frewind(fID);
    %fig.Visible = 'on';
end
fclose(fID);
end

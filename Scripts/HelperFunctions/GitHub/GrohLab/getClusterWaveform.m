function clWaveforms = getClusterWaveform(clusterID, dataDir)
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
% Emilio Isaias-Camacho @GrohLab 2019

%% Input validation
clWaveforms = cell(1,2);
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
try
    fP = fopen(fullfile(dataDir,'params.py'),'r');
    fname = fgetl(fP); fname = strsplit(fname, '= '); fname = fname{2};
    ln = fgetl(fP);
    fclose(fP);
catch
    fprintf(2, 'There''s no params.py file in this folder.\n')
    fprintf(1, 'Perhaps you joined 2 or more simultaneously recorded areas?\n')
    fprintf(1, 'Try joining all_channels.mat files and place it in this folder\n')
    fprintf(1, 'It is also likely that there''s no .bin file to read the waveforms from.\n')
    fprintf(1, 'For now, this function cannot continue\n')
    return
end

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
    clIdx(:,ccl) = strcmp(clTable{:,1}, clusterID(ccl));
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
chanInTable = contains(clTable.Properties.VariableNames,{'ch','channels'});
ch2read = clTable{clusterID, clTable.Properties.VariableNames{chanInTable}};

%% Verifying if the waveform(s) for the given cluster(s) was/were computed

waveFile = dir(fullfile(dataDir,'*_waveforms.mat'));

if ~isempty(waveFile)
    waveFileName = fullfile(dataDir, waveFile.name);
    load(waveFileName,'clWaveforms')
    % N_exCl = size(clWaveforms, 1);
    [exIdx, exSub] = ismember(clWaveforms(:,1),clusterID);
    % 'Removing' unwanted cluster waveforms and retrieving missing clusters
    clWaveforms = clWaveforms(exIdx,1:2);
    % Number of stored clusters
    Nstcl = nnz(exIdx);
    % Number of requested clusters
    Nrqcl = length(clusterID);
    missingClSubs = setdiff(1:Nrqcl, exSub);
    Nbncl = numel(missingClSubs);
    if Nbncl
        % Allocation of the bigger space
        clOutput = cell(Nstcl+Nbncl,2);
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
        switch size(cl_fromBin,2)
            case 3
                clOutput(reqOrdSub > Nstcl,:) = cl_fromBin;
            case 2
                clOutput(reqOrdSub > Nstcl,1:2) = cl_fromBin;
            otherwise
                if ~isempty(cl_fromBin)
                    fprintf(2, 'Something went wrong! Please verify file!\n')
                end
        end
        % Re-assigning the variable names
        clWaveforms = clOutput;
        save(waveFileName, 'clWaveforms','-append')
    end
else
    clWaveforms = fetchWaveforms_fromBin(dataDir, clusterID,...
        chanMap, Nch, clSub, clTempSubs, spkIdx, ch2read);
    % Naming the waveform file
    binFile = dir(fullfile(dataDir,'*.bin'));
    smrxFile = dir(fullfile(dataDir, '*.smrx'));
    if numel(binFile) == 1
        % As the binary file
        [~,baseName] = fileparts(fullfile(binFile.folder,binFile.name));
    elseif numel(smrxFile) == 1
        % As the smrx file
        [~,baseName] = fileparts(fullfile(smrxFile.folder,smrxFile.name));
    elseif numel(binFile) > 1 || numel(smrxFile) > 1
        % User dependent name
        [~, bbaseNames] = arrayfun(@(x) fileparts(x.name), binFile,...
            'UniformOutput', 0);
        [~, sbaseNames] = arrayfun(@(x) fileparts(x.name), smrxFile,...
            'UniformOutput', 0);
        bsBaseNames = [bbaseNames;sbaseNames];
        [bnSub, iOk] = listdlg('ListString', bsBaseNames,...
            'PromptString', 'Output name:',...
            'CancelString', 'None of these');
        if ~iOk
            slfnans = questdlg('Would you like to name it yourself?',...
                'File name','Yes','No','Yes');
            if strcmp(slfnans,'Yes')
                slfName = inputdlg('File name:','Write a title',[16,200],...
                    bsBaseNames{1});
                [~,baseName] = fileparts(slfName);
            end
        else
            baseName = bsBaseNames{bnSub};
        end
    end
    waveFileName = fullfile(dataDir, [baseName, '_waveforms.mat']);
    save(waveFileName, 'clWaveforms')
end
end

function clWaveforms =...
    fetchWaveforms_fromBin(dataDir, clusterID,...
    chanMap, Nch, clSub, clTempSubs, spkIdx, ch2read)
afOpt = {'UniformOutput', 0};
%% Reading the binary file
clWaveforms = cell(numel(clusterID),2);
% Determinig the spike times for the given clusters
spikeFile = dir(fullfile(dataDir,'*_all_channels.mat'));
if ~isempty(spikeFile)
    load(fullfile(dataDir, spikeFile.name), 'sortedData', 'fs')
    if ~exist('fs','var')
        fsFile = dir(fullfile(dataDir,'*_sampling_frequency.mat'));
        load(fullfile(dataDir, fsFile.name), 'fs')
    end
end
% Taking ~1.25 ms around the spike.
spikeWaveTime = 2*round(1.25e-3 * fs) + 1;
spikeSamples = (spikeWaveTime - 1)/2;
binFile = dir(fullfile(dataDir, '*.bin'));
if isempty(binFile)
    fprintf(1, 'Without a binary file it is impossible to get the waveforms')
    fprintf(1,'\n');
    return
end
% pcFeat = readNPY(fullfile(dataDir, 'pc_features.npy'));
% pcInd = readNPY(fullfile(dataDir, 'pc_feature_ind.npy'));
spkSubs = cellfun(@(x) round(x.*fs),sortedData(clSub,2),...
    'UniformOutput',false);
% [chs2read, readOrder, repeatChs] = unique(ch2read);
chs2read = unique(ch2read);
answ = 1; 
if numel(binFile) > 1
    [answ, iOk] = listdlg(...
        'ListString', arrayfun(@(x) x.name, binFile,'UniformOutput', 0),...
        'SelectionMode', 'single');
    if ~iOk
        strAns = questdlg('Are you sure you want to cancel?','Quit?',...
            'Yes','No','No');
        if strcmpi(strAns,'Yes')
            fprintf(1, 'Quitting!\n')
            return
        end
        answ = 1;
    end
end
fID = fopen(fullfile(dataDir, binFile(answ).name), 'r');
cchan = 1;
% Main loop
while ~feof(fID) && cchan <= size(chs2read,1)
    % Computing the location of the channel features
    % pcIdx = ch2read(cchan) == chanMap(pcInd(clTempSubs{cchan}+1,:)+1);
    % clFeat = pcFeat(spkIdx(:,cchan), :, pcIdx);
    fprintf(1,'Reading channel %d ',chs2read(cchan))
    % Jumping to the channel
    fseek(fID, 2*(ch2read(cchan)), 'bof');
    % Getting all clusters from the considered channel
    clustChanIdx = ch2read == chs2read(cchan); Nccl = sum(clustChanIdx);
    Nspks = cellfun(@(x) size(x,1), spkSubs(clustChanIdx));
    spkLbls = arrayfun(@(x) x*ones(Nspks(x),1), 1:Nccl,...
        'UniformOutput', 0); spkLbls = cat(1, spkLbls{:});
    chSpks = [cat(1, spkSubs{clustChanIdx}), spkLbls];
    [ordSpks, spkOrd] = sort(chSpks(:,1), 'ascend');
    % Computing the distance from spike to spike
    spkDists = [ordSpks(1);diff(ordSpks)];
    fprintf(1,'looking for cluster%s...', sprintf(' %s', clusterID{clustChanIdx}))
    % Allocating space for the spikes
    waveform = zeros(spikeWaveTime, sum(Nspks), 'single');
    %fig = figure('Color',[1,1,1],'Visible', 'off');
    %ax = axes('Parent', fig); ax.NextPlot = 'add';
    %subSet = 1:floor(numel(spkDists)*0.1);
    for cspk = 1:sum(Nspks)
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
    clWaveforms(clustChanIdx,:) = [clusterID(clustChanIdx),...
        arrayfun(@(x) waveform(:,spkOrd(chSpks(:,2) == x)), (1:Nccl)', afOpt{:})];
    cchan = cchan + 1;
    frewind(fID);
    %fig.Visible = 'on';
end
fclose(fID);
end

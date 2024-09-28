function [ sortedData, fs ] = ...
    importPhyFiles(pathToData, outputName, removeNoise)
%IMPORTPHYFILES reads the files created by both Kilosort and Phy or Phy2
%and, depending on the flag configurations, it saves a .mat file. The
%contents of this file are constructed with the output files (.tsv, .npy).
%           importPhyFiles( pathToData, outputName, removeNoise,
%               ChAndAmpFlag)
%       INPUTS: no marking means required, [] mean optional but needed for
%                                             later input arguments.
%           pathToData - Path string to where the .tsv and .npy files are
%                           located
%           [outputName] - Name for the .mat file to be created.
%                           WARNING! if the file already exists, the
%                           function will overwrite it
%           [removeNoise] - Default false; flag to indicate the inclusion
%                           of the clusters labelled as noise.
%       OUTPUTS:
%           No variable imported to the workspace but a .mat files written
%           in the given folder; either in the pathToData or in the
%           outputName.
%  Alexandra Merkel & Emilio Isaias-Camacho @ GrohLab 2019
% Maximum number of groups. 3 so far: 1 - good, 2 - MUA, and 3 - noise

MX_CLSS = 3;
% Take Kilosort Data and create *all_channels.map for further analysis

% What is your outputfile called
if exist(fullfile(pathToData,'recording.dat'),'file') && ~exist('outputName','var')
    outputName = 'recording';
elseif ~exist('outputName','var')

    tempDir = pathToData;
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
    outputName = strrep(rhdFile(1).name,'analysis.mat','');
end
% Input sampling frequency (not written in the phy files from KiloSort)
try
    load(fullfile( pathToData, [char(outputName), '_sampling_frequency.mat']),...
        'fs')
catch
    
    try
        load(fullfile(pathToData,'rez.mat'),'rez')
        fs = rez.ops.fs;
    catch
        
        try
            settings = parseXML(fullfile(pathToData,'settings.xml'));
            idx = contains({settings.Attributes.Name},'SampleRateHertz');
            fs = str2double(settings.Attributes(idx).Value);
        catch
            
            try
                text = fileread(fullfile(pathToData,'params.py'));
                fs = str2double(regexp(regexp(text, 'sample_rate = \d*', 'match','once'),'\d*','match'));
            catch LE
                disp( LE.message)
                fsString = inputdlg('Please provide a sampling frequency:');
                
                try
                    fs = str2double(fsString{1});
                    if isnan(fs)
                        fprintf('I''m sorry... you should put in ONLY NUMBERS :)\nStart again\n')
                        fprintf('No output file written\n')
                        return
                    end
                catch
                    fprintf(1,'Cancel button pressed. ')
                    fprintf(1,'Transforming the files might help you get the sampling ')
                    fprintf(1,'frequency.\n');
                    return
                end
            end
        end
    end
    
    save(fullfile( pathToData, [char(outputName), '_sampling_frequency.mat']),...
        'fs')
end

% Do you want to remove Noise Data?
if ~exist('removeNoise','var')
    removeNoise = false;
end

% Reading the KiloSort datafiles (phy updates these files)
spkTms = readNPY(fullfile(pathToData, 'spike_times.npy'));
clID = readNPY(fullfile(pathToData, 'spike_clusters.npy'));

if exist(fullfile(pathToData,'cluster_group.tsv.v2'),'file')
    clGr = readTSV(fullfile(pathToData,'cluster_group.tsv.v2'));
else
    clGr = readTSV(fullfile(pathToData,'cluster_group.tsv'));
end

if ~exist(fullfile(pathToData, 'cluster_info.tsv'),'file')
    fprintf('\nNo cluster_info.tsv file.\n')
    fprintf('\nAssuming automated curation...\n')
end

% Creating a logical index for the user labels
nIdx = cellfun(@strcmp,clGr(:,2),repmat("noise",size(clGr,1),1));
gIdx = cellfun(@strcmp,clGr(:,2),repmat("good",size(clGr,1),1));
mIdx = cellfun(@strcmp,clGr(:,2),repmat("mua",size(clGr,1),1));
grIdx = [gIdx, mIdx, nIdx];

allClusters = cellfun(@num2str,clGr(:,1),'UniformOutput',false);

spkCell = cell(size(clGr,1),1);
for ccln = 1:numel(allClusters)
    spkIdx = clID == clGr{ccln,1};
    spkCell(ccln) = {double(spkTms(spkIdx))/fs};
end

row = find(grIdx');
[r,col2] = ind2sub(size(grIdx),row);
[~,rearrange] = sort(r);
clGroup = col2(rearrange);
[row2, ~] = find(grIdx);
cuts = sum(grIdx,1);
limits = [0,cumsum(cuts)];
for cgroup = 1:MX_CLSS
    clGroup(row2(limits(cgroup)+1:limits(cgroup+1))) = cgroup;
end
sortedData = cat(2,allClusters,spkCell,num2cell(clGroup));
% Removes the noise from the matrix and saves an alternative output file
if removeNoise
    sortedData(nIdx,:) = [];
end
filename = [char(outputName), '_all_channels.mat'];
fname = fullfile(pathToData, filename);
if exist(fname, 'file')
    saveAns = questdlg(sprintf('%s exists already! Overwrite?',filename),...
        'Overwrite?', 'Yes', 'No', 'Yes');
    if strcmp(saveAns, 'No')
        fprintf(1, 'No file written!\n')
        return
    end
end
save(fname, 'sortedData', 'fs');
end
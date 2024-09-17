function clusterInfo = getClusterInfo(filename)
% GETCLUSTERINFO reads the phy generated file 'cluster_info.tsv' and
% returns a table with the variables assigned into an output table
% 'clusterInfo'. It will have size NxM, with N clusters and M variables.
%       clusterInfo = getClusterInfo(filename)
%           INPUTS
%               - filename - string indicating the relative or absolute
%               path to the 'cluster_info.tsv' file.
%           OUTPUTS
%               - clusterInfo - an NxM table with N rows corresponding to
%               each cluster and M columns corresponding to each variable
%               in the file (e.g. id, amplitude, channel, etc.)
% Emilio Isaias-Camacho @ GrohLab 2020
%% Initial and auxiliary variables
getContFirstCell = @(x) x{1};
clusterInfo = table;

%% Validation of the input argument; the filename
% Checking for existence
if ~exist(filename,'file')
    fprintf(1,'The given file doesn''t exist\n')
    return
end
% Checking for validity
[~, ~, fext] = fileparts(filename);
fext = char(fext);
if ~strcmp(fext,'.tsv')
    fprintf(1,'Wrong file fed to the function!\n')
    fprintf(1,'It must be the phy-generated file ''cluster_info.tsv''\n')
    fprintf(1,'Try again with the correct file\n')
    return
end
%% Reading the file
% Open the file for reading
fID = fopen(filename,'r');
% Read the header for each column
heads = getContFirstCell(textscan(fgetl(fID), '%s\t'));
Nv = length(heads);
% Read all the values in the file
conts = getContFirstCell(textscan(fID, '%s', 'Delimiter', '\t',...
    'Whitespace', ' '));
emptyIdx = cellfun(@isempty, conts); conts(emptyIdx) = {'0'};
% Close the file and process the contents
fclose(fID);
%% Scanning the extracted contents 
Ne = length(conts);
conts_rs = reshape(conts,Nv,Ne/Nv)';
for cv = 1:Nv
    % Proceeding slightly different for each variable
    switch heads{cv}
        case {'cluster_id','id','KSLabel','group','NeuronType','Region'}
            % Read as string if the header indicates strings
            vals = cellfun(@(x) getContFirstCell(textscan(x,'%s')),...
                conts_rs(:,cv));
        case 'firing_rate'
            % Read specially as the values are accompanied by the units
            vals = cell2mat(cellfun(@(x) textscan(x,'%f spk/s'),...
                conts_rs(:,cv)));
        otherwise
            % Read simply as a number for the rest of the variables
            vals = cell2mat(cellfun(@(x) textscan(x,'%f'),...
                conts_rs(:,cv)));
    end
    % Assigning the values to the corresponding variable for all custers
    clusterInfo.(heads{cv}) = vals;
end
% Naming the rows and columns appropiately and respectively.
clusterInfo.Properties.DimensionNames = {'Clusters', 'Measures'};
clusterInfo.Properties.RowNames = clusterInfo.(heads{1});
end

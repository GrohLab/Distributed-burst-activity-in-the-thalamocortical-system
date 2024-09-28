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
% Filippo Heimburg @ GrohLab 2020 (adapted from Emilio Isaias-Camacho 2020)

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

[data, header, raw] = tsvread(filename);
data = data(:,~cellfun(@isempty, header));
raw = raw(:,~cellfun(@isempty, header));
header = header(~cellfun(@isempty, header));

strIdx = all(isnan(data)) | ismember(header,{'cluster_id','id'});

varTypes = repmat({'double'},1,numel(header));
varTypes(strIdx) = {'char'};

% Suppress character warning, since they need to be char arrays
warning('off','MATLAB:table:PreallocateCharWarning')
clusterInfo = table('Size',[size(data,1)-1 size(data,2)],'VariableTypes',varTypes,'VariableNames',header);
warning('on','MATLAB:table:PreallocateCharWarning')

clusterInfo(:,strIdx) = raw(2:end,strIdx);
clusterInfo{:,~strIdx} = data(2:end,~strIdx);

% Rename the cluster_id variable for historic reasons
clusterInfo.Properties.DimensionNames = {'Clusters', 'Measures'};
try
    clusterInfo.Properties.RowNames = clusterInfo.id;
catch
    clusterInfo = renamevars(clusterInfo,"cluster_id","id");
    clusterInfo.Properties.RowNames = clusterInfo.id;
end

end

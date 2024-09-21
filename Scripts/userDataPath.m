%% Add path to the data directory

% This assumes that the Data dir is saved in the default Documents dir
documentsPath = fullfile(getenv('USERPROFILE'), 'Documents');
cohortPath = fullfile(documentsPath, 'Data');

% Extract the folder path of this script to save the variable
scriptFullPath = matlab.desktop.editor.getActiveFilename();
scriptFolder = fileparts(scriptFullPath);
save(fullfile(scriptFolder,'userDataPath.mat'),'cohortPath');
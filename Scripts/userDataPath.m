%% Add path to the data directory

% This assumes that the Data dir is saved in the default Documents dir
documentsPath = fullfile(getenv('USERPROFILE'), 'Documents');
cohortPath = fullfile(documentsPath, 'Distributed-burst-activity-in-the-thalamocortical-system/Data');

% Extract the folder path of this script to save the variable
scriptFullPath = matlab.desktop.editor.getActiveFilename();
scriptFolder = fileparts(scriptFullPath);
save(fullfile(scriptFolder,'userDataPath.mat'),'cohortPath');

% Update the session paths in the animalData.mat file
load(fullfile(cohortPath, 'animalData.mat'))
for animal = 1:height(animalData.cohort(12).animal(:))
    name = animalData.cohort(12).animal(animal).animalName;
    nameIdx = regexp(animalData.cohort(12).animal(animal).session_names{1},name,'start');
    animalData.cohort(12).animal(animal).session_names = cellfun(@(x) fullfile(cohortPath, x(nameIdx:end)), animalData.cohort(12).animal(animal).session_names, 'UniformOutput', false);
end
save(fullfile(cohortPath, 'animalData.mat'),'animalData')
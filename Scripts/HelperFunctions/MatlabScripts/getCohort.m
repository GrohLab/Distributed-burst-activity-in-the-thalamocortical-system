function cohort = getCohort(animalPath)

% [~,cohortString] = fileparts(fileparts(animalPath));
cohortString = regexp(animalPath,'Cohort\d*','match','once');
cohort_num = regexp(cohortString,'\d*','match','once');

cohort = ['cohort_',cohort_num];

end
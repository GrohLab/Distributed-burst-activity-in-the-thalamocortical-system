function [isFieldResult, substruct] = fieldInStruct (inStruct, fieldName)
% inStruct is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
isFieldResult = 0;
substruct = NaN;
f = fieldnames(inStruct(1));
for i=1:length(f)
if(strcmp(f{i},strtrim(fieldName)))
isFieldResult = 1;
return;
elseif isstruct(inStruct(1).(f{i}))
isFieldResult = fieldInStruct(inStruct(1).(f{i}), fieldName);
if isFieldResult
    substruct = f{i};
return;
end
end
end
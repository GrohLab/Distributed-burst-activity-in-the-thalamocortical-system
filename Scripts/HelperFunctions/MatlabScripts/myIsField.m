function fieldExist = myIsField(inStruct, fieldName)
% inStruct is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
fieldExist = 0;
f = fieldnames(inStruct(1));
for i=1:length(f)
    if(strcmp(f{i},strtrim(fieldName)))
        fieldExist = 1;
        %         fieldOut = f{i};
        return
    elseif isstruct(inStruct(1).(f{i}))
        fieldExist = myIsField(inStruct(1).(f{i}), fieldName);
        %         fieldOut = inStruct(1).(f{i});
        if fieldExist
            return
        end
    end
end
end
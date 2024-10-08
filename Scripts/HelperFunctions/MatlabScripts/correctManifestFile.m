function correctManifestFile(manifest_file, attributes_file)
% Some manifest files, created by syntalos (especially the ones that were
% created between 2021-05-18 and 2021-06-15), have false timestamp indices (data_aux.parts).
% This function serves as to correct these errors.

manifest_struct = read(manifest_file);

% Check if two data field were created for data_aux (timestamps)
if myIsField(manifest_struct,'data_aux') && numel(manifest_struct.data_aux) > 1
    % Add time zone to the datetime
    manifest_struct.time_created.TimeZone = '+02:00';
    
    % Add the missing description
    Order = {'file_type','summary','parts'};
    manifest_struct.data_aux{1,1}.summary = 'Video timestamps';
    manifest_struct.data_aux{1,1} = orderfields(manifest_struct.data_aux{1,1}, Order);
    
    % Correct the shifted index count
    for sec = 1:numel(manifest_struct.data_aux{1, 1}.parts)
        manifest_struct.data_aux{1, 1}.parts{1,sec}.index = ...
            manifest_struct.data_aux{1, 1}.parts{1,sec}.index + 1;
    end
    
    % Add the first video section
    manifest_struct.data_aux{1,1}.parts = {struct('fname', ...
        manifest_struct.data_aux{1,2}.parts{1,1}.fname,...
        'index', 0),manifest_struct.data_aux{1,1}.parts{1,1:end}};
    
    % Delete the second data_aux field
    manifest_struct.data_aux = manifest_struct.data_aux{1,1};
    
    % Sometimes the timestamps files were doubled, so delete the doubles,
    % if present.
    if numel(manifest_struct.data_aux.parts) > numel(manifest_struct.data.parts)
        manifest_struct.data_aux.parts = manifest_struct.data_aux.parts(1:numel(manifest_struct.data.parts));
    end
    
    % Rename old manifest.toml
    [fold, ~, ext] = fileparts(manifest_file);
    movefile(manifest_file,fullfile(fold,strcat('original_mani',ext)))
    write(manifest_file, manifest_struct)
    
elseif myIsField(manifest_struct,'data_aux') && numel(manifest_struct.data_aux{:}.parts) > numel(manifest_struct.data.parts)
    manifest_struct.data_aux{1,1}.parts = manifest_struct.data_aux{1,1}.parts(1:numel(manifest_struct.data.parts));
    
    % Add time zone to the datetime
    manifest_struct.time_created.TimeZone = '+02:00';
    
    % Rename old manifest.toml
    [fold, ~, ext] = fileparts(manifest_file);
    movefile(manifest_file,fullfile(fold,strcat('original_mani',ext)))
    write(manifest_file, manifest_struct)
end

% Some attributes files generated strange values without a key-pair.
% This function removes those loners.
attributes_struct = read(attributes_file);
[bool, substruct] = fieldInStruct(attributes_struct,'artefact');
if bool
    add = attributes_struct.video.artefact - floor(attributes_struct.video.artefact);
    attributes_struct.video.framerate = attributes_struct.video.framerate + add;
    % 12 digit precision
    attributes_struct.video.framerate = sprintf('%.12f',attributes_struct.video.framerate);
       
    attributes_struct.(substruct) = ...
        rmfield(attributes_struct.(substruct),'artefact');
    
    % Rename old attributes.toml
    [fold, ~, ext] = fileparts(attributes_file);
    movefile(attributes_file,fullfile(fold,strcat('original_attrib',ext)))
    write(attributes_file, attributes_struct)
end

end

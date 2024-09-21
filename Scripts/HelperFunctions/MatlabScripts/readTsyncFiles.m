function timestamps = readTsyncFiles(filePath)

% Add python script to path
pe = pyenv;
if logical(pe.Status) % If python environment is not loaded, then load
    pyenv('ExecutionMode','InProcess',...
        'Version','C:\Users\Groh\anaconda3\python.exe');
elseif ~strcmp(pe.Executable,"C:\Users\Groh\anaconda3\pythonw.exe")
    %     terminate(pyenv)
    pyenv('ExecutionMode','InProcess',...
        'Version','C:\Users\Groh\anaconda3\pythonw.exe');
    %     error('Python environment is not correct (base). You have to restart Matlab to load this environment.')
end

pathToReadFile = fileparts("C:/Users/Groh/PycharmProjects/LabWork/read_tsync_file.py");
if count(py.sys.path,pathToReadFile) == 0
    insert(py.sys.path,int32(0),pathToReadFile)
end
py.importlib.import_module('read_tsync_file');

if contains(filePath,'tis-camera')
    pathSession = fileparts(fileparts(fileparts(filePath)));
    ts = py.read_tsync_file.get_timestamps_overview(pathSession);
    timestamps = double(ts{1});
elseif contains(filePath,'intan-signals')
    fprintf('\nLoading synchronized timestamps...\n')
    pathSession = fileparts(fileparts(filePath));
    sync = py.read_tsync_file.get_timestamps_intan(pathSession);
    if isequal(sync.units.char,'millisecond')
        timestamps = double(sync.magnitude);
    elseif isequal(sync.units.char,'second')
        timestamps = double(sync.magnitude).*1000;
    else
        error('The unit of the timestamps must be millisecond or second.')
    end
    
    % For plotting the registered syntalos tsync timestamps
%     figure, plot(timestamps), title('Tsync timestamps'), ylabel('msec'), xlabel('Sample num')
    %
    
    offset_msec = timestamps(1);
    tsync_samp_count = numel(timestamps);
    % ADD VARIABLE 'continuous' true/false in order to distinguish the
    % continuously mapped and the starting_offset corrected timestamps
    if sum(timestamps==0) > 1
        % Some tsync files seem to stop recording at some point, which
        % results in many zero-value timestamps. These are corrected with
        % fitted values, based on the slope of the trajectory.
        tsync_complete = false;
        
        zero_idx = find(timestamps==0,1);
        P = polyfit((1:zero_idx),timestamps(1:zero_idx),1);
        timestamps(zero_idx:end) = polyval(P,(zero_idx:tsync_samp_count));
    else
        tsync_complete = true;
    end
    
    FileInfo = dir(fullfile(fileparts(filePath),'**\*_sampling_frequency.mat'));
    
    if ~isempty(FileInfo)
        load(fullfile(FileInfo.folder,FileInfo.name),'fs')
    else
        try
            settings = parseXML(fullfile(fileparts(filePath),'settings.xml'));
            idx = contains({settings.Attributes.Name},'SampleRateHertz');
            fs = str2double(settings.Attributes(idx).Value);
        catch
            fs = 30000;
        end
    end
    
    rhd_samp_count = Intan_sampleNum(filePath); 
    
    continuous_ts_msec = ((1:rhd_samp_count).*1000)./fs; 
    sync_ts_msec = timestamps;
%     offset_ts_msec = continuous_ts_msec - (1/fs) + offset_msec;
	Timestamps = table(continuous_ts_msec,sync_ts_msec);
    
    % For plotting the time-delay between intan data and syntalos timestamp
%     diff = sync_ts_msec - continuous_ts_msec;
%     figure, plot(diff), title('Difference between cont and sync timestamps'), ylabel('Difference in msec'), xlabel('Sample num')
    
    % For tracking the various changes, made to the edlio software
    update = 1;

    if rhd_samp_count == tsync_samp_count
        same_samp_count = true;
    else
        same_samp_count = false;
        fprintf('\nrhd_samp_count and tsync_samp_count not equal\n')
    end
    
    save(fullfile(fileparts(filePath),'intanTimestamps.mat'),'Timestamps','offset_msec','tsync_complete','same_samp_count','update')
    
elseif contains(filePath,'hispeed1')
    pathSession = fileparts(fileparts(fileparts(filePath)));
    ts = py.read_tsync_file.get_timestamps_highspeed(pathSession,1);
    timestamps = cell(1,length(ts));
    for i = 1:length(ts)
        timestamps{i} = double(ts{i});
    end
elseif contains(filePath,'hispeed2')
    pathSession = fileparts(fileparts(fileparts(filePath)));
    ts = py.read_tsync_file.get_timestamps_highspeed(pathSession,2);
    timestamps = cell(1,length(ts));
    for i = 1:length(ts)
        timestamps{i} = double(ts{i});
    end
end

end

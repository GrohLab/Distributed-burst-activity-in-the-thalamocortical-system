function HispeedTrials = assignTrialsToVideos(videoDir)
addpath 'Z:\Filippo\Scripts\MatlabScripts\general functions'

% Check if the correct directory was chosen
% The user has to pick the overall 'videos' folder, not a hispeed folder

if ~exist(videoDir,'dir')
    error('The given folder doesn''t exist.')
elseif ~endsWith(videoDir,'videos') 
    error('Incorrect directory. \n%s','Choose the overall ''videos'' folder, not a hispeed folder')
end

% Check if HispeedTrials table was already generated in the past
if exist(fullfile(videoDir,'HispeedTrials.mat'),'file')
    load(fullfile(videoDir,'HispeedTrials.mat'),'HispeedTrials')
    return
end

% Count the number of hispeed sections for preallocation
count = 0;
for highspeedNum = 1:2
    hispeedDir = fullfile(videoDir,['hispeed',num2str(highspeedNum)]);

    % Some manifest files were corrupted, while a new version of Syntalos was
    % tested out. This function corrects the scripts.
    manifest_file = fullfile(hispeedDir,'manifest.toml');
    attributes_file = fullfile(hispeedDir,'attributes.toml');
    correctManifestFile(manifest_file, attributes_file)
    
    list = strtrim(string(ls(hispeedDir)));
    new = sum((endsWith(list,[".mkv",".avi",".mp4"])) & ~contains(list,'DLC')); % Get amount of the hispeed recordings (exclude videos, created with DLC)
    count = count+new;
end

% Make a table with the video paths, numeric section, corresponding hispeed
% camera, logical Go and logical Success
% I had to set the Go or Success marks as doubles to include NaNs

% The variable 'Go_NoGo_Neutral_settingBased' is in most cases equivalent
% to the variable 'Go_NoGo_Neutral'. Only in extinction phases these two
% differ in the actual width of the aperture ('settingBased'),
% and the outcome.
HispeedTrials = table('Size',[count,9],'VariableTypes',...
    {'string','double','double','double','double','double','cell','double','double'},...
    'VariableNames',{'VideoPath','Hispeed','Section','Go_NoGo_Neutral','Go_NoGo_Neutral_settingBased','Lick','Timestamps','Event_Index','PreviousMP'});

% Check the next Go/No-Go event

% Get the event list and the timestamps of that session
[ts,Events,~]=xlsread(fullfile(fileparts(videoDir),'events','table.csv'));
Events=Events(2:end,2:3);
Events(:,1) = num2cell(ts);

% Get the indices of all beam and trial events (the pattern 'Go' will
% index both Go and No-Go Trials).
beam_idx = find(contains(Events(:,2), 'inner side'));
beam_ts = cell2mat(Events(beam_idx,1));
trial_idx = find(contains(Events(:,2), {'Go','Neutral','Diagonal','Medium'}));
drumset_idx{1} = find(contains(Events(:,2), 'Drum 1 set to'));
drumset_idx{2} = find(contains(Events(:,2), 'Drum 2 set to'));
mp_idx = find(contains(Events(:,2), 'Middle point'));

for highspeedNum = 1:2
    empty = find(cellfun('isempty', HispeedTrials{:,'VideoPath'} ),1);
    hispeedDir = fullfile(videoDir,['hispeed',num2str(highspeedNum)]);
    list = strtrim(string(ls(hispeedDir)));
    hispeedPath = list((endsWith(list,[".mkv",".avi",".mp4"])) & ~contains(list,'DLC')); % Get the name of the hispeed recordings (exlude videos, made with DLC)
    hispeedPath = string(natsort(cellstr(hispeedPath))); % This naturally sorts the sections
    
    % One could probably also chronologically count the sections, but that
    % would give errors, when deleting certain files.
    sec_num = regexp(hispeedPath,'\d*\.','Match');
    sec_num = cellfun(@(x) str2double(x),sec_num);
    hispeedPath = fullfile(hispeedDir,hispeedPath);
    
    HispeedTrials{empty:empty+numel(hispeedPath)-1,'VideoPath'} = hispeedPath;
    HispeedTrials{empty:empty+numel(hispeedPath)-1,'Hispeed'} = highspeedNum;
    HispeedTrials{empty:empty+numel(hispeedPath)-1,'Section'} = sec_num;
    
    % Get the timestamps of each section in a cell array
    % Older versions saved the timestamps as .csv files, more recent ones as .tsync files
    if any(contains(list,'.tsync'))
        try timestamps = readTsyncFiles(hispeedPath{1});
        catch except
            % This is highly risky, as you are about to add computed time
            % stamps to a corrupted file. If this is feasible (e.g., if no
            % spike data is synched to the frames, you can go on and
            % comment the error section out.
            if isequal(except.message,'Python Error: Exception: Block terminator not found: Some data is likely corrupted.')
%                 fprintf(2,'\nIf you want to create your own time stamps, you can comment the error out\n')
%                 fprintf(2,'\nError in: %s | Line: %i\n',except.stack(end-1).name,except.stack(end-1).line)
%                 error(except.message)
                timestamps = cell(1,numel(sec_num));
            end
        end
    else
        timestamps = cell(1,numel(sec_num));
        for sec = 1:numel(sec_num)
            matrix_paths = regexprep(hispeedPath,{'.mkv','.avi'},'_timestamps.csv');
            timestamps{sec} = readmatrix(matrix_paths(sec));
        end
    end 
    
    % If there are empty timestamp files, infer the timestamps from the
    % frames and frame rate
    if any(cellfun(@isempty, timestamps))
        for empty_sec = find(cellfun(@isempty, timestamps))
            obj = VideoReader(hispeedPath{empty_sec});  
            fr = round(obj.FrameRate);
            numFrames = obj.NumFrames;
            if fr ~= 240
                fprintf(2,'\nThe manual timestamp creation, currently only works with a frame rate of 240FPS.\n')
                return
            end
            idx = find(contains(Events(:,2),['Started recording on HiSpeed Cam',num2str(highspeedNum)]),empty_sec);
            start_ts = Events{idx(end),1};
            % With a framerate of 240FPS the average delay between two
            % frames would be 4.1667msec. In other words, every sixth frame
            % diff needs to be 5msec.
            % This leads to the first frame being added an extra 4msec
            % which is fine, because the median delay to the event prompt, 
            % usually is 4msec in un-corrupted tsync files.
            add_time = 4*ones(numFrames,1);
            add_time(6:6:end) = 5;
            timestamps{empty_sec} = [(1:numFrames)',start_ts + cumsum(add_time)];
        end
    end

    for sec = 1:numel(sec_num)
        % Get timestamps of that section
        HispeedTrials{empty+sec-1,'Timestamps'} = {timestamps{sec}(:,2)};
        
        % Compare the closest 'Drums x inner side'
        [~,idx_min] = min(abs(beam_ts-timestamps{sec}(1,2)));
        idx_lick = find(trial_idx > beam_idx(idx_min),1);
        idx_previous = find(mp_idx < beam_idx(idx_min),1,'last');

        previous_time = Events{mp_idx(idx_previous),1};
        HispeedTrials{empty+sec-1,'PreviousMP'} = previous_time;
        
        % Access the last 'Drum x set to' event
        % This is especially relevant for assessing the correct state in
        % the extinction stage.
        idx_set = find(drumset_idx{highspeedNum} < beam_idx(idx_min),1,'last');
        last_event = Events{drumset_idx{highspeedNum}(idx_set),2};

        if isempty(idx_lick) % Could occur at the end of a session
            if contains(last_event,'reward','IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral'} = 1;
            elseif contains(last_event,'punishment','IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral'} = 2;
            elseif contains(last_event,{'neutral','diagonal','medium'},'IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral'} = 3;
            else
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral'} = NaN;
            end
            
            if contains(last_event,'reward','IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = 1;
            elseif contains(last_event,'punishment','IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = 2;
            elseif contains(last_event,{'neutral','diagonal','medium'},'IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = 3;
            else
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = NaN;
            end
            
            HispeedTrials{empty+sec-1,'Lick'} = NaN;
            HispeedTrials{empty+sec-1,'Event_Index'} = NaN;
        % Shouldn't be more than 100ms apart, otherwise probably artefact
        % Find the next trial outcome
        elseif Events{beam_idx(idx_min),1} - timestamps{sec}(1,2) <= 100
            next_event = Events{trial_idx(idx_lick),2};
            event_time = Events{trial_idx(idx_lick),1};
            
            % Find high-speed timestamp that is closest to the actual lick
            % event, since no-go trials are only interrupted 1sec after the
            % lick has occured
            [~,idx_lick_event] = min(abs(event_time-timestamps{sec}(:,2)));
            HispeedTrials{empty+sec-1,'Event_Index'} = idx_lick_event;
            
            % Assign Go/No-Go and Success/Failure
            if contains(next_event,'No-Go')
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral'} = 2;
            elseif contains(next_event,{'Neutral','Diagonal','Medium'})
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral'} = 3;
            else
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral'} = 1;
            end
            
            if contains(last_event,'reward','IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = 1;
            elseif contains(last_event,'punishment','IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = 2;
            elseif contains(last_event,{'neutral','diagonal','medium'},'IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = 3;
            else
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = NaN;
            end
            
            if contains(next_event,'Success') && contains(next_event,'No-Go')
                HispeedTrials{empty+sec-1,'Lick'} = 0;
            elseif contains(next_event,'Failure') && contains(next_event,'No-Go')
                HispeedTrials{empty+sec-1,'Lick'} = 1;            
            elseif contains(next_event,'Success')
                HispeedTrials{empty+sec-1,'Lick'} = 1;
            elseif contains(next_event,'Failure')
                HispeedTrials{empty+sec-1,'Lick'} = 0;
            elseif contains(next_event,{'No-Lick','No Lick'}) && contains(next_event,{'Neutral','Diagonal','Medium'})
                HispeedTrials{empty+sec-1,'Lick'} = 0;
            elseif contains(next_event,{'Neutral','Diagonal','Medium'})
                HispeedTrials{empty+sec-1,'Lick'} = 1;
            end
        else % Beam break and video start are more than 100ms apart
            HispeedTrials{empty+sec-1,'Go_NoGo_Neutral'} = NaN;
            
            % Set the aperture setting anyways
            if contains(last_event,'reward','IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = 1;
            elseif contains(last_event,'punishment','IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = 2;
            elseif contains(last_event,{'neutral','diagonal','medium'},'IgnoreCase',true)
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = 3;
            else
                HispeedTrials{empty+sec-1,'Go_NoGo_Neutral_settingBased'} = NaN;
            end

            HispeedTrials{empty+sec-1,'Lick'} = NaN;
        end
    end

end

end
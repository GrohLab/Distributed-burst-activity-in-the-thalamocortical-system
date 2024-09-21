function HispeedTrials = getWhiskerContacts(sessionDir)

%% Grab the DLC (filtered) csv file with the bodyparts

videoDir = fullfile(sessionDir,'videos');
HispeedTrials = assignTrialsToVideos(videoDir);

pathToScript = fileparts("C:\Users\Groh\PycharmProjects\pythonProject\DLCforMatlab.py");
if count(py.sys.path,pathToScript) == 0
    insert(py.sys.path,int32(0),pathToScript)
end

% Check for each video section whether it has already been analyzed and
% filtered. If not, do it now.
for i = 1:2
    HispeedTrials_temp = HispeedTrials(HispeedTrials.Hispeed==i,:);
    list = strtrim(string(ls(fullfile(videoDir,['hispeed',num2str(i)])))); % List of all files in the respective hispeed directory
    for k = 1:numel(HispeedTrials_temp.VideoPath)
        [~,videoName,~] = fileparts(HispeedTrials_temp.VideoPath{k});
        if any(startsWith(list,[videoName,'DLC']))
            if ~any(startsWith(list,[videoName,'DLC']) & endsWith(list,'filtered.csv'))
                error('The predictions of %s were not filtered yet. Filtering now...\n',videoName)
            end
        else
            error('You need to analyze the videos first, in order to get the whisker contact points.')
        end
    end
end

%% Assessing contact points with drums

% !!! Add: If screwplates are not visible, then average through drum points

% If the Contact points have already been extracted, load mat file

answer = [];
if exist(fullfile(videoDir,'HispeedTrials.mat'),'file')
    answer = questdlg('Whisker Contacts have already been analyzed. Would you like to load the .mat file?', ...
        'Import of previous Analysis', ...
        'Yes','No, re-analyze contact points','Yes');
end

if strcmp(answer,'Yes')
    load(fullfile(videoDir,'HispeedTrials.mat'),'HispeedTrials')
else
    HispeedTrials.ContactLeft = cell(size(HispeedTrials,1),1);
    HispeedTrials.ContactRight = cell(size(HispeedTrials,1),1);
    
    % Get points, where most distal whisker label proximizes the drum
    % |Whisker Label - Middle Point| - Radius < Tolerance (e.g. 3mm)
    for i = 1:2
        fprintf('Analyzing whisker contacts of %d. high-speed camera...\n',i)
        HispeedTrials_temp = HispeedTrials(HispeedTrials.Hispeed==i,:);
        list = strtrim(string(ls(fullfile(videoDir,['hispeed',num2str(i)])))); % List of all files in the respective hispeed directory
        dirc = dir(fullfile(videoDir,['hispeed',num2str(i)]));
        old_video = false;
        try
            date = dirc(ismember({dirc.name}, 'attributes.toml')).date;
            attributes = read(fullfile(videoDir,['hispeed',num2str(i)],'attributes.toml'));
            frame_dimensions(1) = attributes.frame_height; % Y-Coord
            frame_dimensions(2) = attributes.frame_width;  % X-Coord
        catch
            old_video = true;
            idx = find(endsWith(list,'meta.pickle'),1);
            metadata_file = fullfile(videoDir,['hispeed',num2str(i)],list(idx));
            fid = py.open(metadata_file,'rb');
            data = struct(py.pickle.load(fid));
            p = data.data{'frame_dimensions'};
            % Convert python tuple into array -> Dimensions in (y,x) because
            % high-speed cameras are turned 90 degrees
            frame_dimensions = cellfun(@double,cell(p));
        end
        
        for k = 1:size(HispeedTrials_temp,1)
            [~,videoName,~] = fileparts(HispeedTrials_temp.VideoPath{k});
            tokens = regexp(videoName, '_sec(\d+(\.\d+)?)', 'tokens');
            secNum = str2double(tokens{1});
            csvFile = list(startsWith(list,[videoName,'DLC']) & endsWith(list,'filtered.csv'));
            % If there is more than one filtered video, chose the latest one
            if numel(csvFile) > 1
               idx = find(ismember({dirc.name}, csvFile));
               [~,latest] = max([dirc(idx).datenum]);
               csvFile = dirc(idx(latest)).name;
               csvFile = fullfile(videoDir,['hispeed',num2str(i)],csvFile);
            else
                csvFile = fullfile(videoDir,['hispeed',num2str(i)],csvFile);
            end
            
            raw_table = readtable(csvFile,'Range','A2','VariableNamingRule','preserve');
            raw_table = raw_table(:,2:end); % Excludes the first column of indices
            headers = raw_table.Properties.VariableNames;
            headers = strrep(headers,'-',''); % DLC was trained with '-' assignments, which are not supported in Matlab
            raw_table.Properties.VariableNames = headers;
            
            % Mirror videos from hispeed2, when they were recorded before 2021,
            % as they were recorded with inverted x- and y-axis
            if (contains(date,'2020') || old_video) && i == 2
                for ii = 1:numel(headers)
                    if ~endsWith(headers{ii},["_1","_2"])
                        raw_table{:,ii} = frame_dimensions(1) - raw_table{:,ii}; % Y-coordinates
                        raw_table{:,ii+1} = frame_dimensions(2) - raw_table{:,ii}; % X-coordinates
                    end
                end
            end

            % Align the axis to the horizontal and set the middle point as zero
            % This is beeing done, by taking the axis through the lickports and the center
            % of the two backlight plates and the middle point markers
            
            % Get medians of stationary points with points >90% certainty
            screwplate_front = [median(raw_table.screwplate_front(raw_table.screwplate_front_2>0.9)),...
                median(raw_table.screwplate_front_1(raw_table.screwplate_front_2>0.9))];
            screwplate_back = [median(raw_table.screwplate_back(raw_table.screwplate_back_2>0.9)),...
                median(raw_table.screwplate_back_1(raw_table.screwplate_back_2>0.9))];
            
            % Fit an axis through the screw plate points
            p = polyfit([screwplate_front(2),screwplate_back(2)],...
                [screwplate_front(1),screwplate_back(1)],1);
            theta = atan(p(1));
            
            if theta > 0
                a = -pi;
            else
                a = pi;
            end
            
            % Rotation along a given rotation point(x1|y1)
            x_transformation = @(x,y,x1,y1,theta) cos(theta)*(x - x1)-sin(theta)*(y - y1);
            y_transformation = @(x,y,x1,y1,theta) sin(theta)*(x - x1)+cos(theta)*(y - y1);
            
            tilt_table = raw_table;
            for m = 1:numel(headers)
                if ~endsWith(headers{m},["_1","_2"])
                    tilt_table{:,m+1} = x_transformation(raw_table{:,m+1},raw_table{:,m},...
                        screwplate_front(2),screwplate_front(1),a-theta); % X-coordinates
                    tilt_table{:,m} = y_transformation(raw_table{:,m+1},raw_table{:,m},...
                        screwplate_front(2),screwplate_front(1),a-theta); % Y-coordinates
                end
            end
            
            % Grab medians again with the tilted table
            screwplate_front = [median(tilt_table.screwplate_front(tilt_table.screwplate_front_2>0.9)),...
                median(tilt_table.screwplate_front_1(tilt_table.screwplate_front_2>0.9))];
            screwplate_back = [median(tilt_table.screwplate_back(tilt_table.screwplate_back_2>0.9)),...
                median(tilt_table.screwplate_back_1(tilt_table.screwplate_back_2>0.9))];
                        
            % Conversion of pixel values into metric millimeters
            % The distance between the left and the right plate is 76mm
            conversion_fac = 76/(sqrt((screwplate_back(1)-screwplate_front(1))^2+(screwplate_back(2)-screwplate_front(2))^2));
            
            screwplate_front = screwplate_front*conversion_fac;
            screwplate_back = screwplate_back*conversion_fac;
            % Validate the X/Y transformation
            if screwplate_front(2) < screwplate_back(2) || screwplate_front(1)-screwplate_back(1) > 10
%                 error("The X/Y transformation seems to be faulty. Check the tilt angles again.")
            end
            
            metric_table = tilt_table;
            for m = 1:numel(headers)
                if ~endsWith(headers{m},["_1","_2"]) % !!!! WATCH OUT FOR HEADERS WITH UNDERSCORES
                    metric_table{:,m} = conversion_fac*tilt_table{:,m}; % X-coordinates
                    metric_table{:,m+1} = conversion_fac*tilt_table{:,m+1}; % Y-coordinates
                end
            end
            
            drumleft1 = [median(metric_table.drumleft1(metric_table.drumleft1_2>0.9)),...
                median(metric_table.drumleft1_1(metric_table.drumleft1_2>0.9))];
            drumleft2 = [median(metric_table.drumleft2(metric_table.drumleft2_2>0.9)),...
                median(metric_table.drumleft2_1(metric_table.drumleft2_2>0.9))];
            drumleft3 = [median(metric_table.drumleft3(metric_table.drumleft3_2>0.9)),...
                median(metric_table.drumleft3_1(metric_table.drumleft3_2>0.9))];
            drumleft4 = [median(metric_table.drumleft4(metric_table.drumleft4_2>0.9)),...
                median(metric_table.drumleft4_1(metric_table.drumleft4_2>0.9))];
            drumleft5 = [median(metric_table.drumleft5(metric_table.drumleft5_2>0.9)),...
                median(metric_table.drumleft5_1(metric_table.drumleft5_2>0.9))];
            
            drumright1 = [median(metric_table.drumright1(metric_table.drumright1_2>0.9)),...
                median(metric_table.drumright1_1(metric_table.drumright1_2>0.9))];
            drumright2 = [median(metric_table.drumright2(metric_table.drumright2_2>0.9)),...
                median(metric_table.drumright2_1(metric_table.drumright2_2>0.9))];
            drumright3 = [median(metric_table.drumright3(metric_table.drumright3_2>0.9)),...
                median(metric_table.drumright3_1(metric_table.drumright3_2>0.9))];
            drumright4 = [median(metric_table.drumright4(metric_table.drumright4_2>0.9)),...
                median(metric_table.drumright4_1(metric_table.drumright4_2>0.9))];
            drumright5 = [median(metric_table.drumright5(metric_table.drumright5_2>0.9)),...
                median(metric_table.drumright5_1(metric_table.drumright5_2>0.9))];
            
            
            dl = [drumleft1; drumleft2; drumleft3; drumleft4; drumleft5];
            dl = dl(all(~isnan(dl),2),:);
            dr = [drumright1; drumright2; drumright3; drumright4; drumright5];
            dr = dr(all(~isnan(dr),2),:);
            
            [x_leftdrum, y_leftdrum, rad_leftdrum] = circle_fit(dl(:,2),dl(:,1));
            [x_rightdrum, y_rightdrum, rad_rightdrum] = circle_fit(dr(:,2),dr(:,1));
            
            %         % Plot the fitted circles
                    figure
                    hold on
                    plot(dr(:,2),dr(:,1),'r*')
                    th = 0:pi/50:2*pi;
                    xunit = rad_rightdrum * cos(th) + x_rightdrum;
                    yunit = rad_rightdrum * sin(th) + y_rightdrum;
                    plot(xunit, yunit,'r');
            
                    plot(dl(:,2),dl(:,1),'m*')
                    th = 0:pi/50:2*pi;
                    xunit = rad_leftdrum * cos(th) + x_leftdrum;
                    yunit = rad_leftdrum * sin(th) + y_leftdrum;
                    plot(xunit, yunit,'m');
                    plot([screwplate_front(2) screwplate_back(2)],[screwplate_front(1) screwplate_back(1)],'go')
            
                    hold off
                    axis equal
            
            % Extract all frames, in which either the first or the second
            % whisker of each side is recognized with high probability.
            whisker_left_safe = find((metric_table.whisker2_left6_2 > 0.8 & metric_table.whisker2_left5_2 > 0.8 & metric_table.whisker2_left4_2 > 0.8) | ...
                (metric_table.whisker1_left6_2 > 0.8 & metric_table.whisker1_left5_2 > 0.8 & metric_table.whisker1_left4_2 > 0.8));
            whisker_right_safe = find((metric_table.whisker2_right6_2 > 0.8 & metric_table.whisker2_right5_2 > 0.8 & metric_table.whisker2_right4_2 > 0.8) | ...
                (metric_table.whisker1_right6_2 > 0.8 & metric_table.whisker1_right5_2 > 0.8 & metric_table.whisker1_right4_2 > 0.8));
            
            leftdrum_touch_msec = NaN(1,numel(whisker_left_safe));
            tolerance = 4; % Tolerance for whisker contact in mm
            for m = 1:numel(whisker_left_safe)
                % Check if the distance between first left whisker tip and left drum is < 2mm
                if sqrt((metric_table.whisker1_left6_1(whisker_left_safe(m))-x_leftdrum)^2+(metric_table.whisker1_left6(whisker_left_safe(m))-y_leftdrum)^2) ...
                        < (rad_leftdrum + tolerance)
                    touch_event = whisker_left_safe(m); % This is in frames
                    touch_event = HispeedTrials_temp.Timestamps{k}(touch_event); % This converts it to a timestamp in msec
                    leftdrum_touch_msec(m) = touch_event;
                    % Check if the distance between second left whisker tip and left drum is < 2mm
                elseif sqrt((metric_table.whisker2_left6_1(whisker_left_safe(m))-x_leftdrum)^2+(metric_table.whisker2_left6(whisker_left_safe(m))-y_leftdrum)^2) ...
                        < (rad_leftdrum + tolerance)
                    touch_event = whisker_left_safe(m); % This is in frames
                    touch_event = HispeedTrials_temp.Timestamps{k}(touch_event); % This converts it to a timestamp in msec
                    leftdrum_touch_msec(m) = touch_event;
                end
            end
            leftdrum_touch_msec = leftdrum_touch_msec(~isnan(leftdrum_touch_msec));
            if ~isempty(leftdrum_touch_msec)
                idx = find(HispeedTrials.Hispeed == i & HispeedTrials.Section == secNum);
                HispeedTrials.ContactLeft(idx) = {leftdrum_touch_msec};
            end
            
            rightdrum_touch_msec = NaN(1,numel(whisker_right_safe));
            for m = 1:numel(whisker_right_safe)
                % Check if the distance between first right whisker tip and right drum is < 2mm
                if sqrt((metric_table.whisker1_right6_1(whisker_right_safe(m))-x_rightdrum)^2+(metric_table.whisker1_right6(whisker_right_safe(m))-y_rightdrum)^2) ...
                        < (rad_rightdrum + tolerance)
                    touch_event = whisker_right_safe(m); % This is in frames
                    touch_event = HispeedTrials_temp.Timestamps{k}(touch_event); % This converts it to a timestamp in msec
                    rightdrum_touch_msec(m) = touch_event;
                    % Check if the distance between second right whisker tip and right drum is < 2mm
                elseif sqrt((metric_table.whisker2_right6_1(whisker_right_safe(m))-x_rightdrum)^2+(metric_table.whisker2_right6(whisker_right_safe(m))-y_rightdrum)^2) ...
                        < (rad_rightdrum + tolerance)
                    touch_event = whisker_right_safe(m); % This is in frames
                    touch_event = HispeedTrials_temp.Timestamps{k}(touch_event); % This converts it to a timestamp in msec
                    rightdrum_touch_msec(m) = touch_event;
                end
            end
            rightdrum_touch_msec = rightdrum_touch_msec(~isnan(rightdrum_touch_msec));
            if ~isempty(rightdrum_touch_msec)
                idx = find(HispeedTrials.Hispeed == i & HispeedTrials.Section == secNum);
                HispeedTrials.ContactRight(idx) = {rightdrum_touch_msec};
            end
        end
    end
    
    save(fullfile(videoDir,'HispeedTrials.mat'),'HispeedTrials')
    
end
end
function HispeedTrials = getWhiskerContactsApertures(sessionDir)
configPath = 'C:\Users\Groh\DLC-Modules\Whiskertrack_Apertures-Heimburg-2021-11-11\config.yaml';

%% Grab the DLC (filtered) csv file with the bodyparts

videoDir = fullfile(sessionDir,'videos');
HispeedTrials = assignTrialsToVideos(videoDir);

try
    pyenv('ExecutionMode','InProcess',...
        'Version','C:\Users\Groh\anaconda3\envs\DEEPLABCUT\pythonw.exe');
    pathToScript = 'C:\Users\Groh\PycharmProjects\LabWork\DLCforMatlab.py';
    if count(py.sys.path,fileparts(pathToScript)) == 0
        insert(py.sys.path,int32(0),fileparts(pathToScript))
    end
    py.importlib.import_module('DLCforMatlab');
catch
    fprintf('\nPython environment (DEEPLABCUT) could not be loaded. You have to restart Matlab to load this environment.\n')
    fprintf('Ignoring for now, assuming that all videos have been analyzed aleady.\n\n')
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
                % If not filtered, filter the data and display in terminal
                fprintf('The predictions of %s were not filtered yet. Filtering now...\n',videoName)

                % Filter predictions
                py.DLCforMatlab.filterpredictions(configPath,...
                    HispeedTrials_temp.VideoPath{k});
            end
        else
            % If not all videos are analzyed quit the script
            fprintf('\nThe video %s was not analyzed yet. Please analyze all videos, before you continue with this script.\n',videoName)
            return
        end
    end
end

%% Assessing contact points with drums

HispeedTrials.ContactLeft = cell(size(HispeedTrials,1),1);
HispeedTrials.ContactRight = cell(size(HispeedTrials,1),1);

% Get points, where most distal whisker label proximizes the aperture wings
% |Whisker Label - Wing| < Tolerance (e.g. 3mm)
for i = 1:2
    fprintf('Analyzing whisker contacts of %d. high-speed camera...\n',i)
    HispeedTrials_temp = HispeedTrials(HispeedTrials.Hispeed==i,:);
    list = strtrim(string(ls(fullfile(videoDir,['hispeed',num2str(i)])))); % List of all files in the respective hispeed directory
    dirc = dir(fullfile(videoDir,['hispeed',num2str(i)]));

    attributes = read(fullfile(videoDir,['hispeed',num2str(i)],'attributes.toml'));
    frame_dimensions(1) = attributes.video.frame_width;  % X-Coord
    frame_dimensions(2) = attributes.video.frame_height; % Y-Coord

    wing_right_reserve = cell(3,1);
    wing_left_reserve = cell(3,1);

    for k = 1:size(HispeedTrials_temp,1)
        [~,videoName,~] = fileparts(HispeedTrials_temp.VideoPath{k});
        tokens = regexp(videoName, '_sec(\d+(\.\d+)?)', 'tokens');
        secNum = str2double(tokens{1});
        csvFile_whisker = list(startsWith(list,[videoName,'DLC']) & endsWith(list,'filtered.csv') & contains(list,'Whiskertrack_Apertures'));
        csvFile_stationary = list(startsWith(list,[videoName,'DLC']) & endsWith(list,'filtered.csv') & contains(list,'Nosetrack_Apertures'));

        % If there is more than one filtered video, chose the latest one
        if numel(csvFile_whisker) == 0
            error('MyScript:NoCsvFile','There is no csv file for the whisker tracking of: %s',HispeedTrials_temp.VideoPath{k})
        elseif numel(csvFile_whisker) > 1
            idx = find(ismember({dirc.name}, csvFile_whisker));
            [~,latest] = max([dirc(idx).datenum]);
            csvFile_whisker = dirc(idx(latest)).name;
            csvFile_whisker = fullfile(videoDir,['hispeed',num2str(i)],csvFile_whisker);
        else
            csvFile_whisker = fullfile(videoDir,['hispeed',num2str(i)],csvFile_whisker);
        end

        if numel(csvFile_stationary) == 0
            error('MyScript:NoCsvFile','There is no csv file for the nose tracking of: %s',HispeedTrials_temp.VideoPath{k})
        elseif numel(csvFile_stationary) > 1
            idx = find(ismember({dirc.name}, csvFile_stationary));
            [~,latest] = max([dirc(idx).datenum]);
            csvFile_stationary = dirc(idx(latest)).name;
            csvFile_stationary = fullfile(videoDir,['hispeed',num2str(i)],csvFile_stationary);
        else
            csvFile_stationary = fullfile(videoDir,['hispeed',num2str(i)],csvFile_stationary);
        end

        raw_table = readtable(csvFile_whisker,'Range','A2','VariableNamingRule','preserve');
        raw_table = raw_table(:,2:end); % Excludes the first column of indices
        headers = raw_table.Properties.VariableNames;
        headers = strrep(headers,'-',''); % DLC was trained with '-' assignments, which are not supported in Matlab
        raw_table.Properties.VariableNames = headers;

        stationary_table = readtable(csvFile_stationary,'Range','A2','VariableNamingRule','preserve');
        stationary_table = stationary_table(:,2:end); % Excludes the first column of indices
        headers = stationary_table.Properties.VariableNames;
        headers = strrep(headers,'-',''); % DLC was trained with '-' assignments, which are not supported in Matlab
        stationary_table.Properties.VariableNames = headers;

        % Get medians of stationary points with points greater than the defined certainty
        median_table = table('Size',[2, 13],'VariableTypes',...
            repmat({'double'},1,13),'VariableNames',...
            {'wing_left_corner','wing_left_base_1','wing_left_base_2',...
            'wing_left_edge_1','wing_left_edge_2','wing_right_corner',...
            'wing_right_base_1','wing_right_base_2','wing_right_edge_1',...
            'wing_right_edge_2','lickport','lp_edge_left','lp_edge_right'},'RowNames',{'X','Y'});

        certainty = 0.8;

        median_table.wing_left_corner(:) = [median(stationary_table.wing_left_corner(stationary_table.wing_left_corner_2>certainty)),...
            median(stationary_table.wing_left_corner_1(stationary_table.wing_left_corner_2>certainty))];
        median_table.wing_left_base_1(:) = [median(stationary_table.wing_left_base_1(stationary_table.wing_left_base_1_2>certainty)),...
            median(stationary_table.wing_left_base_1_1(stationary_table.wing_left_base_1_2>certainty))];
        median_table.wing_left_base_2(:) = [median(stationary_table.wing_left_base_2(stationary_table.wing_left_base_2_2>certainty)),...
            median(stationary_table.wing_left_base_2_1(stationary_table.wing_left_base_2_2>certainty))];
        median_table.wing_left_edge_1(:) = [median(stationary_table.wing_left_edge_1(stationary_table.wing_left_edge_1_2>certainty)),...
            median(stationary_table.wing_left_edge_1_1(stationary_table.wing_left_edge_1_2>certainty))];
        median_table.wing_left_edge_2(:) = [median(stationary_table.wing_left_edge_2(stationary_table.wing_left_edge_2_2>certainty)),...
            median(stationary_table.wing_left_edge_2_1(stationary_table.wing_left_edge_2_2>certainty))];

        median_table.wing_right_corner(:) = [median(stationary_table.wing_right_corner(stationary_table.wing_right_corner_2>certainty)),...
            median(stationary_table.wing_right_corner_1(stationary_table.wing_right_corner_2>certainty))];
        median_table.wing_right_base_1(:) = [median(stationary_table.wing_right_base_1(stationary_table.wing_right_base_1_2>certainty)),...
            median(stationary_table.wing_right_base_1_1(stationary_table.wing_right_base_1_2>certainty))];
        median_table.wing_right_base_2(:) = [median(stationary_table.wing_right_base_2(stationary_table.wing_right_base_2_2>certainty)),...
            median(stationary_table.wing_right_base_2_1(stationary_table.wing_right_base_2_2>certainty))];
        median_table.wing_right_edge_1(:) = [median(stationary_table.wing_right_edge_1(stationary_table.wing_right_edge_1_2>certainty)),...
            median(stationary_table.wing_right_edge_1_1(stationary_table.wing_right_edge_1_2>certainty))];
        median_table.wing_right_edge_2(:) = [median(stationary_table.wing_right_edge_2(stationary_table.wing_right_edge_2_2>certainty)),...
            median(stationary_table.wing_right_edge_2_1(stationary_table.wing_right_edge_2_2>certainty))];

        median_table.lickport(:) = [median(stationary_table.lickport(stationary_table.lickport_2>certainty)),...
            median(stationary_table.lickport_1(stationary_table.lickport_2>certainty))];
        median_table.lp_edge_left(:) = [median(stationary_table.lp_edge_left(stationary_table.lp_edge_left_2>certainty)),...
            median(stationary_table.lp_edge_left_1(stationary_table.lp_edge_left_2>certainty))];
        median_table.lp_edge_right(:) = [median(stationary_table.lp_edge_right(stationary_table.lp_edge_right_2>certainty)),...
            median(stationary_table.lp_edge_right_1(stationary_table.lp_edge_right_2>certainty))];


        % Conversion of pixel values into metric millimeters
        % The distance between the left and the right lickport edge is 30mm

        conversion_fac = 30/(sqrt((median_table.lp_edge_left(1)-median_table.lp_edge_right(1))^2+(median_table.lp_edge_left(2)-median_table.lp_edge_right(2))^2));
        metric_table = raw_table;

        metric_table{:,1:3:end} = conversion_fac*raw_table{:,1:3:end}; % X-coordinates
        metric_table{:,2:3:end} = conversion_fac*raw_table{:,2:3:end}; % Y-coordinates

        median_table{1,:} = conversion_fac*median_table{1,:}; % X-coordinates
        median_table{2,:} = conversion_fac*median_table{2,:}; % Y-coordinates


        % Define boundaries of wings
        % NaN values have to be excluded, otherwise the polyfit will only return NaNs
        xvals = [median_table.wing_left_base_2(1),median_table.wing_left_base_1(1),median_table.wing_left_corner(1)];
        xvals = xvals(~isnan(xvals));
        yvals = [median_table.wing_left_base_2(2),median_table.wing_left_base_1(2),median_table.wing_left_corner(2)];
        yvals = yvals(~isnan(yvals));
        p_left.base = polyfit(xvals,yvals,1);

        xvals = [median_table.wing_left_edge_2(1),median_table.wing_left_edge_1(1),median_table.wing_left_corner(1)];
        xvals = xvals(~isnan(xvals));
        yvals = [median_table.wing_left_edge_2(2),median_table.wing_left_edge_1(2),median_table.wing_left_corner(2)];
        yvals = yvals(~isnan(yvals));
        p_left.edge = polyfit(xvals,yvals,1);

        if p_left.edge(2) > median_table.wing_left_corner(2)
            upper = -(p_left.edge(2)/p_left.edge(1));
            wing_left = [median_table.wing_left_corner(1) 0 0 upper;...
                median_table.wing_left_corner(2) p_left.base(2) 0 0];
        else
            wing_left = [median_table.wing_left_corner(1) 0 0;...
                median_table.wing_left_corner(2) p_left.base(2) p_left.edge(2)];
        end

        % If the wing is not recognized correctly, then use previous
        % stationary points
        if ~any(any(isnan(wing_left))) && ~isnan(HispeedTrials_temp.Go_NoGo_Neutral(k))
            wing_left_reserve{HispeedTrials_temp.Go_NoGo_Neutral(k)} = wing_left;
        elseif ~isnan(HispeedTrials_temp.Go_NoGo_Neutral(k))
            wing_left = wing_left_reserve{HispeedTrials_temp.Go_NoGo_Neutral(k)};
        end

        xvals = [median_table.wing_right_base_2(1),median_table.wing_right_base_1(1),median_table.wing_right_corner(1)];
        xvals = xvals(~isnan(xvals));
        yvals = [median_table.wing_right_base_2(2),median_table.wing_right_base_1(2),median_table.wing_right_corner(2)];
        yvals = yvals(~isnan(yvals));
        p_right.base = polyfit(xvals,yvals,1);

        xvals = [median_table.wing_right_edge_2(1),median_table.wing_right_edge_1(1),median_table.wing_right_corner(1)];
        xvals = xvals(~isnan(xvals));
        yvals = [median_table.wing_right_edge_2(2),median_table.wing_right_edge_1(2),median_table.wing_right_corner(2)];
        yvals = yvals(~isnan(yvals));
        p_right.edge = polyfit(xvals,yvals,1);

        fact = frame_dimensions(1)*conversion_fac;
        if p_right.edge(2) < median_table.wing_right_corner(2)
            lower = p_right.base(1)*fact + p_right.base(2);
            upper = -(p_right.edge(2)/p_right.edge(1));
            wing_right = [median_table.wing_right_corner(1) fact fact upper;...
                median_table.wing_right_corner(2) lower 0 0];
        else
            lower = p_right.base(1)*fact + p_right.base(2);
            upper = p_right.edge(1)*fact + p_right.edge(2);
            wing_right = [median_table.wing_right_corner(1) fact fact;...
                median_table.wing_right_corner(2) lower upper];
        end

        % If the wing is not recognized correctly, then use previous
        % stationary points
        if ~any(any(isnan(wing_right))) && ~isnan(HispeedTrials_temp.Go_NoGo_Neutral(k))
            wing_right_reserve{HispeedTrials_temp.Go_NoGo_Neutral(k)} = wing_right;
        elseif ~isnan(HispeedTrials_temp.Go_NoGo_Neutral(k))
            wing_right = wing_right_reserve{HispeedTrials_temp.Go_NoGo_Neutral(k)};
        end

        % Extract all frames, in which either the first or the second
        % whisker of each side is recognized with high probability.
        whisker_left_safe = find((metric_table.whisker2_left8_2 > certainty & metric_table.whisker2_left7_2 > certainty & metric_table.whisker2_left6_2 > certainty) | ...
            (metric_table.whisker1_left8_2 > certainty & metric_table.whisker1_left7_2 > certainty & metric_table.whisker1_left6_2 > certainty));
        whisker_right_safe = find((metric_table.whisker2_right8_2 > certainty & metric_table.whisker2_right7_2 > certainty & metric_table.whisker2_right6_2 > certainty) | ...
            (metric_table.whisker1_right8_2 > certainty & metric_table.whisker1_right7_2 > certainty & metric_table.whisker1_right6_2 > certainty));

        tolerance = 3; % Tolerance for whisker contact in mm

        left_touch_msec = NaN(1,numel(whisker_left_safe));
        for m = 1:numel(whisker_left_safe)
            % Get the distance between whisker tip and aperture with the
            % corresponding timestamp

            % Check if the distance between one of the last three whisker labels and
            % either of the aperture wings is < tolerance
            try
                if p_poly_dist(metric_table.whisker1_left8(whisker_left_safe(m)), metric_table.whisker1_left8_1(whisker_left_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_left8(whisker_left_safe(m)), metric_table.whisker2_left8_1(whisker_left_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_left7(whisker_left_safe(m)), metric_table.whisker1_left7_1(whisker_left_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_left7(whisker_left_safe(m)), metric_table.whisker2_left7_1(whisker_left_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_left6(whisker_left_safe(m)), metric_table.whisker1_left6_1(whisker_left_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_left6(whisker_left_safe(m)), metric_table.whisker2_left6_1(whisker_left_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_left8(whisker_left_safe(m)), metric_table.whisker1_left8_1(whisker_left_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_left8(whisker_left_safe(m)), metric_table.whisker2_left8_1(whisker_left_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_left7(whisker_left_safe(m)), metric_table.whisker1_left7_1(whisker_left_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_left7(whisker_left_safe(m)), metric_table.whisker2_left7_1(whisker_left_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_left6(whisker_left_safe(m)), metric_table.whisker1_left6_1(whisker_left_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_left6(whisker_left_safe(m)), metric_table.whisker2_left6_1(whisker_left_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance

                    %                     dist_left = p_poly_dist(metric_table.nosetip(whisker_left_safe(m)), metric_table.nosetip_1(whisker_left_safe(m)),wing_left(1,:),wing_left(2,:));
                    touch_event = whisker_left_safe(m); % This is in frames
                    touch_event = HispeedTrials_temp.Timestamps{k}(touch_event); % This converts it to a timestamp in msec
                    left_touch_msec(m) = touch_event;
                end
            catch
            end
        end
        left_touch_msec = left_touch_msec(~isnan(left_touch_msec));
        if ~isempty(left_touch_msec)
            idx = find(HispeedTrials.Hispeed == i & HispeedTrials.Section == secNum);
            HispeedTrials.ContactLeft(idx) = {left_touch_msec};
        end

        right_touch_msec = NaN(1,numel(whisker_right_safe));
        for m = 1:numel(whisker_right_safe)
            % Get the distance between whisker tip and aperture with the
            % corresponding timestamp

            % Check if the distance between one of the whisker tips and
            % either of the aperture wings is < tolerance
            try
                if p_poly_dist(metric_table.whisker1_right8(whisker_right_safe(m)), metric_table.whisker1_right8_1(whisker_right_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_right8(whisker_right_safe(m)), metric_table.whisker2_right8_1(whisker_right_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_right7(whisker_right_safe(m)), metric_table.whisker1_right7_1(whisker_right_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_right7(whisker_right_safe(m)), metric_table.whisker2_right7_1(whisker_right_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_right6(whisker_right_safe(m)), metric_table.whisker1_right6_1(whisker_right_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_right6(whisker_right_safe(m)), metric_table.whisker2_right6_1(whisker_right_safe(m)),wing_right(1,:),wing_right(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_right8(whisker_right_safe(m)), metric_table.whisker1_right8_1(whisker_right_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_right8(whisker_right_safe(m)), metric_table.whisker2_right8_1(whisker_right_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_right7(whisker_right_safe(m)), metric_table.whisker1_right7_1(whisker_right_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_right7(whisker_right_safe(m)), metric_table.whisker2_right7_1(whisker_right_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker1_right6(whisker_right_safe(m)), metric_table.whisker1_right6_1(whisker_right_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance || ...
                        p_poly_dist(metric_table.whisker2_right6(whisker_right_safe(m)), metric_table.whisker2_right6_1(whisker_right_safe(m)),wing_left(1,:),wing_left(2,:)) < tolerance
                    
                    touch_event = whisker_right_safe(m); % This is in frames
                    touch_event = HispeedTrials_temp.Timestamps{k}(touch_event); % This converts it to a timestamp in msec
                    right_touch_msec(m) = touch_event;
                end
            catch
            end
        end
        right_touch_msec = right_touch_msec(~isnan(right_touch_msec));
        if ~isempty(right_touch_msec)
            idx = find(HispeedTrials.Hispeed == i & HispeedTrials.Section == secNum);
            HispeedTrials.ContactRight(idx) = {right_touch_msec};
        end
    end
end

save(fullfile(videoDir,'HispeedTrials.mat'),'HispeedTrials')

end
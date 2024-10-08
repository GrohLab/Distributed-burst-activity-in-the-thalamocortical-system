function [Conditions, Triggers] = getConditions(sessionDir, varargin)
% Supplement the Highspeed Table with the individual on- and offset times
% of each trigger class, i.e. reward, punishment, drumcontact

% Default option values:
options = inputParser;

addRequired(options,'sessionDir',@ischar);
addParameter(options,'useAnalyzedData',true,@islogical)
parse(options,sessionDir,varargin{:})
useAnalyzedData = options.Results.useAnalyzedData;

intanDir = fullfile(sessionDir,'intan-signals');

rhdStruct = dir(fullfile(intanDir,'*.rhd'));
[~, rhdName, ~] = fileparts(rhdStruct.name);
rhdFile = fullfile(intanDir,rhdStruct.name);

if exist(fullfile(intanDir,'ADC_Data.mat'),'file') && useAnalyzedData
    load(fullfile(intanDir,'ADC_Data.mat'),'digitalReward', 'digitalPunish', 'digitalLick')
else
    [digitalReward, digitalPunish, digitalLick] = intanADC(rhdFile);
end


% Combine the left and right LP signals for the continuous Trigger signal
digitalReward_cont = digitalReward(2,:);
digitalReward_cont(logical(digitalReward(1,:))) = 1;

digitalPunish_cont = digitalPunish(2,:);
digitalPunish_cont(logical(digitalPunish(1,:))) = 1;

digitalLick_cont = digitalLick(2,:);
digitalLick_cont(logical(digitalLick(1,:))) = 1;

Triggers.Reward = digitalReward_cont;
Triggers.Punishment = digitalPunish_cont;
Triggers.Lick = digitalLick_cont;

% Convert continuous triggers to subscripts and then convert the sampling frequency,
% so that it matches the amplifier sampling frequency, which is usually 30kHz.
% Digital triggers are by default processed with 10kHz.
try
    load(fullfile(intanDir,[rhdName, '_sampling_frequency.mat']),...
        'fs')
catch
    
    try
        load(fullfile(intanDir,'rez.mat'),'rez')
        fs = rez.ops.fs;
    catch
        
        try
            settings = parseXML(fullfile(intanDir,'settings.xml'));
            idx = contains({settings.Attributes.Name},'SampleRateHertz');
            fs = str2double(settings.Attributes(idx).Value);
        catch
            
            try
                text = fileread(fullfile(intanDir,'params.py'));
                fs = str2double(regexp(regexp(text, 'sample_rate = \d*', 'match','once'),'\d*','match'));
            catch LE
                disp( LE.message)
                fsString = inputdlg('Please provide a sampling frequency:');
                
                try
                    fs = str2double(fsString{1});
                    if isnan(fs)
                        fprintf('I''m sorry... you should put in ONLY NUMBERS :)\nStart again\n')
                        fprintf('No output file written\n')
                        return
                    end
                catch
                    fprintf(1,'Cancel button pressed. ')
                    fprintf(1,'Transforming the files might help you get the sampling ')
                    fprintf(1,'frequency.\n');
                    return
                end
            end
        end
    end
    
    save(fullfile(intanDir,[rhdName, '_sampling_frequency.mat']),'fs')
end

fac = fs/10000; % Factor by which the digital signals have to multiplied, to extrapolate to the amplifier fs

rewardOnset = sort([(strfind(digitalReward(1,:),[0 1]) + 1) (strfind(digitalReward(2,:),[0 1]) + 1)])'*fac; % This finds all trigger onsets
rewardOffset = sort([strfind(digitalReward(1,:),[1 0]) strfind(digitalReward(2,:),[1 0])])'*fac; % This finds all trigger offsets
if numel(rewardOnset) ~= numel(rewardOffset)
    rewardOnset = rewardOnset(1:end-1);
end

punishOnset = sort([(strfind(digitalPunish(1,:),[0 1]) + 1) (strfind(digitalPunish(2,:),[0 1]) + 1)])'*fac;
punishOffset = sort([strfind(digitalPunish(1,:),[1 0]) strfind(digitalPunish(2,:),[1 0])])'*fac;
if numel(punishOnset) ~= numel(punishOffset)
    punishOnset = punishOnset(1:end-1);
end

lickOnset = sort([(strfind(digitalLick(1,:),[0 1]) + 1) (strfind(digitalLick(2,:),[0 1]) + 1)])'*fac;
lickOffset = sort([strfind(digitalLick(1,:),[1 0]) strfind(digitalLick(2,:),[1 0])])'*fac;
if numel(lickOnset) ~= numel(lickOffset)
    lickOnset = lickOnset(1:end-1);
end

% The whisker contacts are listed in timestamps [msec] and converted into sample numbers.
HispeedTrials = getWhiskerContactsApertures(sessionDir);

% Set warning to error in order to catch it
s = warning('error', 'MATLAB:load:variableNotFound');
try load(fullfile(intanDir,'intanTimestamps.mat'),'Timestamps','update')
    if update < 1
        readTsyncFiles(rhdFile);
        load(fullfile(intanDir,'intanTimestamps.mat'),'Timestamps')
    end
catch
    readTsyncFiles(rhdFile);
    load(fullfile(intanDir,'intanTimestamps.mat'),'Timestamps')
end
% Reset warning
warning(s);

whiskerLeftOnset = NaN(size(HispeedTrials,1),1);
whiskerLeftOffset = NaN(size(HispeedTrials,1),1);
whiskerRightOnset = NaN(size(HispeedTrials,1),1);
whiskerRightOffset = NaN(size(HispeedTrials,1),1);

leftTouchFirst = false(height(HispeedTrials),1);
rightTouchFirst = false(height(HispeedTrials),1);

for i = 1:size(HispeedTrials,1)
    if ~isempty(HispeedTrials.ContactLeft{i}) && isempty(HispeedTrials.ContactRight{i})
        leftTouchFirst(i) = true;
    elseif ~isempty(HispeedTrials.ContactLeft{i}) && HispeedTrials.ContactLeft{i}(1) < HispeedTrials.ContactRight{i}(1)
        leftTouchFirst(i) = true;
    end

    if ~isempty(HispeedTrials.ContactRight{i}) && isempty(HispeedTrials.ContactLeft{i})
        rightTouchFirst(i) = true;
    elseif ~isempty(HispeedTrials.ContactRight{i}) && HispeedTrials.ContactRight{i}(1) < HispeedTrials.ContactLeft{i}(1)
        rightTouchFirst(i) = true;
    end

    if numel(HispeedTrials.ContactLeft{i}) > 1
        [~,whiskerLeftOnset(i)] = min(abs(Timestamps.sync_ts_msec-HispeedTrials.ContactLeft{i}(1)));
        for m = 1:numel(HispeedTrials.ContactLeft{i})
            count = find(HispeedTrials.ContactLeft{i} > HispeedTrials.ContactLeft{i}(m) & HispeedTrials.ContactLeft{i} < HispeedTrials.ContactLeft{i}(m)+200, 1);
            if isempty(count) && HispeedTrials.ContactLeft{i}(m) ~= HispeedTrials.ContactLeft{i}(1)
                [~,whiskerLeftOffset(i)] = min(abs(Timestamps.sync_ts_msec-HispeedTrials.ContactLeft{i}(m)));
                break
            end
        end
    end
    if numel(HispeedTrials.ContactRight{i}) > 1
        [~,whiskerRightOnset(i)] = min(abs(Timestamps.sync_ts_msec-HispeedTrials.ContactRight{i}(1)));
        for m = 1:numel(HispeedTrials.ContactRight{i})
            count = find(HispeedTrials.ContactRight{i} > HispeedTrials.ContactRight{i}(m) & HispeedTrials.ContactRight{i} < HispeedTrials.ContactRight{i}(m)+200, 1);
            if isempty(count) && HispeedTrials.ContactRight{i}(m) ~= HispeedTrials.ContactRight{i}(1)
                [~,whiskerRightOffset(i)] = min(abs(Timestamps.sync_ts_msec-HispeedTrials.ContactRight{i}(m)));
                break
            end
        end
    end
end

whiskerLeftOnset_onlyFirstLeft = whiskerLeftOnset(~isnan(whiskerLeftOnset) & leftTouchFirst);
whiskerLeftOffset_onlyFirstLeft = whiskerLeftOffset(~isnan(whiskerLeftOffset) & leftTouchFirst);
whiskerRightOnset_onlyFirstRight = whiskerRightOnset(~isnan(whiskerRightOnset) & rightTouchFirst);
whiskerRightOffset_onlyFirstRight = whiskerRightOffset(~isnan(whiskerRightOffset) & rightTouchFirst);

whiskerLeftOnset = whiskerLeftOnset(~isnan(whiskerLeftOnset));
whiskerLeftOffset = whiskerLeftOffset(~isnan(whiskerLeftOffset));
whiskerRightOnset = whiskerRightOnset(~isnan(whiskerRightOnset));
whiskerRightOffset = whiskerRightOffset(~isnan(whiskerRightOffset));

idx = unique(sort(cell2mat(arrayfun(@(x) find(lickOnset/fs*1000 > x, 1), HispeedTrials.PreviousMP, 'UniformOutput', false))));

% Add 200 random triggers for shuffled analyses
sampleNum = Intan_sampleNum(rhdFile);
% Ignore time before first middle point and last 10sec
randOnset = randi([min(HispeedTrials.PreviousMP)*fs/1000, sampleNum-fs*10],200,1);
% Offset will be 100msec after the random onset
randOffset = randOnset+0.1*fs;

conditions = {'Reward','Punishment','Lick','WhiskerContact_left','WhiskerContact_right',...
    'Middlepoint','Random','WhiskerContact_onlyLeftFirst','WhiskerContact_onlyRightFirst','onlyFirstLick'};

triggers = {[rewardOnset rewardOffset], [punishOnset punishOffset], [lickOnset lickOffset], ...
    [whiskerLeftOnset whiskerLeftOffset], [whiskerRightOnset whiskerRightOffset], ...
    [HispeedTrials.PreviousMP*fs/1000 HispeedTrials.PreviousMP*fs/1000], ...
    [randOnset randOffset], ...
    [whiskerLeftOnset_onlyFirstLeft whiskerLeftOffset_onlyFirstLeft], ...
    [whiskerRightOnset_onlyFirstRight whiskerRightOffset_onlyFirstRight], ...
    [lickOnset(idx) lickOffset(idx)]};

% In the beginning of the session the lickports sometimes send out a digital pulses while initializing.
% To exclude such artefacts all triggers within in the first 3sec are removed.
triggers = cellfun(@(x) x(all(x > fs*3,2),:), triggers, 'UniformOutput', false);


Conditions = struct('name', cell(1, numel(conditions)), 'Triggers', cell(1, numel(conditions)));
for i = 1:numel(conditions)*2 % We have to copy the Conditions for the jittering script, as it is sensitive for paired stimuli
    if i > numel(conditions)
        j = i-numel(conditions);
        Conditions(i).name = conditions{j};
        Conditions(i).Triggers = triggers{j};
    else
        Conditions(i).name = conditions{i};
        Conditions(i).Triggers = triggers{i};
    end
end

Conditions(end+1).name = 'AllTriggers';
Conditions(end).Triggers = unique(sort(vertcat(Conditions(1:numel(conditions)).Triggers)),'rows');

% In order to track changes to the trigger calculation
update = 5;
save(fullfile(intanDir,[rhdName, 'analysis.mat']),'Conditions','Triggers','update')

fprintf('\nDone!\n')

end
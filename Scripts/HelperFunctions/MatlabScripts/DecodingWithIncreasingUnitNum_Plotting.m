%% Correlating decoding accuracy with number of input units
% Script that runs multiple decoding analyses with increasing number of
% units. By bootstrapping each step, one can define the average decoding
% accuracy of each step.
close all; clearvars; clc

% Choose sessions to plot
cohortDir = 'Z:\Filippo\Animals\Cohort12_33-38';
FileInfo = dir(fullfile(cohortDir,'Decoding-with-increasing-cellNum'));
FileInfo = FileInfo(3:end,:); FileInfo = FileInfo([FileInfo.isdir],:);

trialTypeList = {'allTrials','onlyGo','onlyNoGo','onlyNeutral','onlyNarrow','onlyWide','onlyIntermediate','onlyLick','onlyNoLick'};

% Pool of cells, you want to choose from:
% 'All cells', 'Only responding with burst bias',
% 'Only responding with tonic bias', 'Only touch-modulated cells'
unitTypeList = {'allUnits','burstUnits','tonicUnits','touchModul'};

% Decide how many conditions you want to plot in one figure
condNum = listdlg('Name','Conditions','PromptString','How many conditions do you want to plot in one figure?',...
    'ListString',arrayfun(@num2str, (1:4), 'UniformOutput', false),'ListSize',[300 150], 'SelectionMode','single');

% Create user interface to make the user pick the conditions
f = figure('Name','Compare Conditions');
ui_field = gobjects(condNum*3,1);
ui_text = gobjects(4+condNum*2,1);

ui_text(1) = uicontrol(f,'Style','text','Units','normalized',...
    'HorizontalAlignment','left','Position',[0.05 0.85 0.8 0.1],...
    'FontSize',10);
ui_text(1).String = 'For each condition, pick the desired parameters';

ui_text(2) = uicontrol(f,'Style','text','Units','normalized',...
    'HorizontalAlignment','left','Position',[0.15 0.75 0.25 0.1],...
    'FontSize',10,'String','Stage');
ui_text(3) = uicontrol(f,'Style','text','Units','normalized',...
    'HorizontalAlignment','left','Position',[0.42 0.75 0.25 0.1],...
    'FontSize',10,'String','Trial type');
ui_text(4) = uicontrol(f,'Style','text','Units','normalized',...
    'HorizontalAlignment','left','Position',[0.69 0.75 0.2 0.1],...
    'FontSize',10,'String','Response type');

count = 1;
for i = 1:condNum
    ui_text(i+4) = uicontrol(f,'Style','text','Units','normalized',...
        'HorizontalAlignment','left','Position',[0.05 0.85-0.1*i 0.2 0.05],...
        'FontSize',10,'String',sprintf('Cond %i',i));

    % Specify which stages you want to analyze
    ui_field(count) = uicontrol(f,'Style','popupmenu','Units','normalized',...
        'HorizontalAlignment','left','Position',[0.15 0.85-0.1*i 0.25 0.05],...
        'FontSize',10,'String',{FileInfo.name});
    count = count+1;

    % Specify which trial types you want to analyze
    ui_field(count) = uicontrol(f,'Style','popupmenu','Units','normalized',...
        'HorizontalAlignment','left','Position',[0.42 0.85-0.1*i 0.25 0.05],...
        'FontSize',10,'String',trialTypeList);
    count = count+1;

    % Specify which cell response types you want to analyze
    ui_field(count) = uicontrol(f,'Style','popupmenu','Units','normalized',...
        'HorizontalAlignment','left','Position',[0.69 0.85-0.1*i 0.2 0.05],...
        'FontSize',10,'String',unitTypeList);
    count = count+1;

end

doneButton = uicontrol(f,'Style','pushbutton','units','normalized',...
    'Position',[0.8 0.1 0.1 0.05],'String','Done',...
    'Callback',{@doneExeCond,f,ui_field,condNum,FileInfo,trialTypeList,unitTypeList});

waitfor(doneButton)

% Label the conditions
condLabels = cell(1,condNum);
for condIdx = 1:condNum
    inputStr = sprintf('%s_%s_%s',stageDescription{condIdx},trialType{condIdx},unitType{condIdx});
    answer = inputdlg(sprintf('How would you like to call your %i. condition?',condIdx),'Labeling',[1 80],{inputStr});
    condLabels{condIdx} = answer{:};
end

inputStr = fullfile(cohortDir,'Decoding-with-increasing-cellNum',sprintf('%s_%s',stageDescription{1},strjoin(condLabels,' & ')));
answer = inputdlg('Enter the name of the figure files, to be saved as.','Labeling',[1 150],{inputStr});
figureFile = answer{:};

% Choose condition to analyze
% Available conditions: 'Reward', 'Punishment', 'Lick','onlyFirstLick',
% 'WhiskerContact_left', 'WhiskerContact_right',
% 'WhiskerContact_onlyLeftFirst','WhiskerContact_onlyRightFirst'
chCond = 'WhiskerContact_left';

%% Decoder parameters
% Define brain areas
area_names = {'BC','VPM','POm','ZIv'};
area_colors = {'#377eb8','#4daf4a','#984ea3','#ff7f00'};

f = figure('Name','Waveform type');
set(gcf,'Position',[1000 600 420 300])
ui_field = gobjects(numel(area_names),1);
ui_text = gobjects(numel(area_names)+1,1);

ui_text(1) = uicontrol(f,'Style','text','Units','normalized',...
    'HorizontalAlignment','left','Position',[0.1 0.85 0.8 0.1],...
    'FontSize',10,'String','Choose the desired subpopulation for each area.');

for i = 1:numel(area_names)
    ui_text(i+1) = uicontrol(f,'Style','text','Units','normalized',...
        'HorizontalAlignment','left','Position',[0.1 0.85-0.1*i 0.1 0.05],...
        'FontSize',10,'String',area_names{i});

    if isequal(area_names{i},'ZIv')
        initVal = 3;
    else
        initVal = 2;
    end
    ui_field(i) = uicontrol(f,'Style','popupmenu','Units','normalized',...
        'HorizontalAlignment','left','Position',[0.2 0.85-0.1*i 0.2 0.05],...
        'FontSize',10,'String',{'all','RS units','FS units'},'Value', initVal);
end

doneButton = uicontrol(f,'Style','pushbutton','units','normalized',...
    'Position',[0.8 0.1 0.1 0.05],'String','Done',...
    'Callback',{@doneExe,f,ui_field,area_names});

waitfor(doneButton)

% answer = questdlg('What do you want to decode?','Decoding object',...
%     'Go - No-go - Neutral','Wide - Narrow - Interm.','Lick - No-lick','Wide - Narrow - Interm.');
% if isequal(answer, 'Go - No-go - Neutral')
%     classifier_labels = 'Go_NoGo_Neutral';
%     if isequal(trialType,'alltrials')
%         trialCompare = {'narrow','wide'};
%     elseif isequal(trialType,'onlyLick')
%         trialCompare = {'narrow&lick','wide&lick'};
%     elseif isequal(trialType,'onlyNoLick')
%         trialCompare = {'narrow&no-lick','wide&no-lick'};
%     else
%         trialCompare = {''};
%     end
% elseif isequal(answer, 'Wide - Narrow - Interm.')
%     classifier_labels = 'Wide_Narrow_Intermediate';
%     if isequal(trialType,'alltrials')
%         trialCompare = {'narrow','wide'};
%     elseif isequal(trialType,'onlyLick')
%         trialCompare = {'narrow&lick','wide&lick'};
%     elseif isequal(trialType,'onlyNoLick')
%         trialCompare = {'narrow&no-lick','wide&no-lick'};
%     else
%         trialCompare = {''};
%     end
% elseif isequal(answer, 'Lick - No-lick')
%     classifier_labels = 'Lick';
%     if isequal(trialType,'onlyNarrow')
%         trialCompare = {'narrow&lick','narrow&no-lick'};
%     elseif isequal(trialType,'onlyWide')
%         trialCompare = {'wide&lick','wide&no-lick'};
%     elseif isequal(trialType,'onlyIntermediate')
%         trialCompare = {'intermediate&lick','intermediate&no-lick'};
%     else
%         trialCompare = {''};
%     end
% end

% Define the classifier to use
% 'max_correlation_coefficient_CL','poisson_naive_bayes_CL','libsvm_CL'
classifierName = 'libsvm_CL';

% Set a desired number of splits. A split of 10 means that 9 repetitions of
% each event are used for training and 1 example is used for testing.
% To get reasonable results you usually need at least 5 repetitions of each
% event (i.e., at least 5 splits)
splitNumber = inputdlg(sprintf('Enter the desired number of splits\n(NaN if you want to set it to default):'),'Split number',[1 50],{'NaN'});
splitNumber = str2double(splitNumber);

% Define the number of bootstraps you want to average
answer = inputdlg('Enter desired unit number of bootstraps to average.','Bootstrap Selection',[1 60],{'20'});
numBootstrap = str2double(answer);

%% Extract decoding information and plot it
respWindow = [0,200]; % response window in msec (default: [0,200])
for ar = 1:numel(area_names)+1 % +1 for all areas
    if ar~=2 % Plot only VPM
        continue
    end
    if ar==numel(area_names)+1
        plotTitle = 'allAreas';
        plotColor = '#000000';
    else
        plotTitle = area_names{ar};
        plotColor = area_colors{ar};
    end
    p = gobjects(1,condNum);
    fig = figure('Name',sprintf('DecodingAccuracyXCellNum_%s',plotTitle));
    hold on

    if condNum==2 % With twoconditions perform statistical testing
        compareVals = cell(2,height(FileInfo));
    end
    for condIdx = 1:condNum
        fprintf('\nExtracting data of %s cells for stage %s and condition %s / %s\n',plotTitle,stageDescription{condIdx},unitType{condIdx},trialType{condIdx})
        
        if isnan(splitNumber)
            resultName = sprintf('Binned_data_results_%s_allSpikes_%s_%s_%s_%s_*splits_bootstrap*',chCond,unitType{condIdx},plotTitle,trialType{condIdx},classifierName);
        else
            resultName = sprintf('Binned_data_results_%s_allSpikes_%s_%s_%s_%s_%isplits_bootstrap*',chCond,unitType{condIdx},plotTitle,trialType{condIdx},classifierName,splitNumber);
        end
        decodingDir = fullfile(cohortDir,'Decoding-with-increasing-cellNum',stageDescription{condIdx},plotTitle,'*cells*');
        FileInfo = dir(decodingDir); % Get all the cell runs
        cellVals = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), {FileInfo.name});

        meanVals = NaN(1,height(FileInfo));
        stdVals = NaN(1,height(FileInfo));
        semVals = NaN(1,height(FileInfo));
        count = 1;

        for cellRun = 1:height(FileInfo)
            BootstrapInfo = dir(fullfile(FileInfo(cellRun).folder, FileInfo(cellRun).name, resultName));
            if isempty(BootstrapInfo)
                % Progress bar
                if cellRun==height(FileInfo)
                    fprintf('100%% done\n')
                elseif cellRun/height(FileInfo)>=0.2*count
                    fprintf('%i%% done\n',20*count)
                    count = count+1;
                end
                continue
            end
            bootstrapVals = cellfun(@(x) regexp(x, 'bootstrap(\d+)', 'tokens'), {BootstrapInfo.name});
            bootstrapVals = cellfun(@str2double, bootstrapVals);
            if ~all(diff(bootstrapVals)==1) % Check if there are missing bootstraps
                idx = find(diff(bootstrapVals)~=1);
                error('DecodingWithIncreasingNum:BootstrapContinuity','Missing bootstrap values between bootstrap %i and %i (and possibly others). Correct that first.',...
                    bootstrapVals(idx(1)), bootstrapVals(idx(1)+1))
            end
            % Check if chosen bootstrap values are represented
            % If not, correct the bootstrap value or abort
            % Only do this with the first run, otherwise you have different values whilst analyzing
            if ~all(ismember((1:numBootstrap), bootstrapVals)) && cellRun==1
                % Correct the number of bootstraps
                answer = questdlg(sprintf('There are less analyzed bootstraps (= %i) than your picked value(= %i). Want to lower the value?',numel(bootstrapVals),numBootstrap), ...
                    'Correct Bootstraps', ...
                    'Yes','No, abort','Yes');
                if isequal(answer,'Yes')
                    numBootstrap = numel(bootstrapVals);
                else
                    return
                end
            elseif ~all(ismember((1:numBootstrap), bootstrapVals))
                error('DecodingWithIncreasingNum:BootstrapNumberNotMatching','The existing bootstrap values can not be corrected twice. Otherwise you''d have different values whilst analyzing.')
            end

            % Load all recording results into one array
            for i = 1:numBootstrap
                idx = contains({BootstrapInfo.name},sprintf('bootstrap%03i',i));
                load(fullfile(BootstrapInfo(idx).folder, BootstrapInfo(idx).name), 'DECODING_RESULTS')
                if i==1
                    params = DECODING_RESULTS.DS_PARAMETERS.binned_site_info.binning_parameters;
                    decodingResults = NaN(numBootstrap, numel(params.the_bin_start_times));
                    xvals = params.the_bin_start_times + params.bin_width/2 - params.alignment_event_time;
                    decodingResults(1,:) = diag(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results);
                else
                    decodingResults(i,:) = diag(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results);
                end
            end

            % Extract mean decoding accuracy in the given response window
            % Register mean and sem of the cell run
            idx = xvals>=respWindow(1) & xvals<respWindow(2);
            meanVals(cellRun) = mean(decodingResults(:,idx),'all','omitmissing');
            stdVals(cellRun) = std(decodingResults(:,idx),[],'all','omitmissing');
            semVals(cellRun) = std(decodingResults(:,idx),[],'all','omitmissing')/sqrt(numel(decodingResults(:,idx)));

            if condNum==2 % With two conditions perform statistical testing
                compareVals{condIdx,cellRun} = decodingResults(:,idx);
            end
        end

        % Plot results
        colorGrad = min(hex2rgb(plotColor) * (1 + 0.2*(condIdx-1)),[1 1 1]); % Create darker shades (max value of 1)

        curve1 = meanVals + semVals;
        curve2 = meanVals - semVals;
        idx = find(~isnan(curve1));
        fill([idx fliplr(idx)], [curve1(idx) fliplr(curve2(idx))],[0 0 .85],...
            'FaceColor',colorGrad,'EdgeColor','none','FaceAlpha',0.2);
        p(condIdx) = plot(meanVals,'Color',colorGrad,'LineWidth',1.5);
        ylim([0.5 1])
        title('Decoding accuracy as a function of cell number')
        xlabel('Cell number')
        ylabel('Decoding accuracy')
        fig.UserData = struct('respondWind', respWindow, ...
            'meanVals', meanVals, ...clo
            'stdVals', stdVals, ...
            'semVals', semVals);
    end

    if condNum==2 % With twoconditions perform statistical testing
        % If some conditions have unequal number of units, only compare the present units
        idx = ~cellfun(@isempty, compareVals);
        compareVals = compareVals(:,all(idx,1));
        h = cellfun(@(x,y) ttest(reshape(x,1,[]),reshape(y,1,[])), compareVals(1,:),compareVals(2,:),'UniformOutput',true);
        
        barBegins = strfind(num2str(h), [0 1]);
        if h(1)
            barBegins = [1, barBegins+1];
        end
        barEnds = strfind(num2str(h), [1 0]);
        if h(end)
            barEnds = [barEnds, numel(h)]; %#ok<AGROW>
        end

        if ~isempty(barBegins)
            for i = 1:numel(barBegins)
                plot([barBegins(i),barEnds(i)],...
                    [0.55, 0.55], 'Color',plotColor,'LineWidth',4)
            end
        end
    end
    ylim([0.5 1])
    legend(p,condLabels,'Location','east')
    hold off

    savefig(fig, [figureFile,'_',fig.Name,'.fig'])
end


%% Helper functions
% Function executed by the "doneButton" in condition refinement
function doneExeCond(~,~,f,ui_field,condNum,FileInfo,trialTypeList,unitTypeList)

stageDescription = cell(1,condNum);
unitType = cell(1,condNum);
trialType = cell(1,condNum);
for condIdx = 1:condNum
    count = (condIdx-1)*3;
    stageDescription{condIdx} = FileInfo(ui_field(count+1).Value).name;
    trialType{condIdx} = trialTypeList{ui_field(count+2).Value};
    unitType{condIdx} = unitTypeList{ui_field(count+3).Value};
end

% Overwrite variables
assignin('caller','stageDescription',stageDescription)
assignin('caller','unitType',unitType)
assignin('caller','trialType',trialType)

close(f)
end

% Function executed by the "doneButton" in subpopulation refinement
function doneExe(~,~,f,ui_field,area_names)

for ar = 1:numel(area_names)
    switch ui_field(ar).Value
        case 2
            area_names{ar} = [area_names{ar} '-RS'];
        case 3
            area_names{ar} = [area_names{ar} '-FS'];
    end
end

% Overwrite variables
assignin('caller','area_names',area_names)

close(f)
end

%% Compares multiple sessions based on the contour maps of burst and tonic scatter maps
% Necessary to run the ApertureResponseTypes.m script for the respective
% stages first.

close all; clearvars; clc;
% Define areas to plot
area_names = {'BC','VPM','POm','ZIv'};
area_colors = {'#377eb8','#4daf4a','#984ea3','#ff7f00'};

% Choose condition to analyze
prompt = {'Reward', 'Punishment', 'Lick', 'onlyFirstLick',...
    'WhiskerContact_left', 'WhiskerContact_right','Middlepoint',...
    'Random', 'WhiskerContact_onlyLeftFirst', 'WhiskerContact_onlyRightFirst', 'onlyFirstLick'};
condPick = listdlg('PromptString','Pick condition to analyze.', ...
    'ListString',prompt,...
    'ListSize',[300 200],'InitialValue',5,'SelectionMode','single');
condition = prompt{condPick};

% Pick trial types to compare
trialType = {'All', 'Wide', 'Wide&lick', 'Wide&no-lick',...
    'Narrow', 'Narrow&lick', 'Narrow&no-lick',...
    'Intermediate', 'Intermediate&lick', 'Intermediate&no-lick'};
[trialPick, tf] = listdlg('ListString',trialType,...
    'PromptString','Which trial type(s) do you want to analyze?',...
    'SelectionMode','multiple','ListSize', [250 250]);
if tf==0
    return
end
assert(numel(trialPick)==2,...
    'singleUniteBurstiness:selectInputNumber',...
    'You have to pick 2 trial types.')
trialType = trialType(trialPick);

% Pick stages to compare
scriptFullPath = matlab.desktop.editor.getActiveFilename();
try load(regexprep(scriptFullPath, 'Scripts.*', 'Scripts\\userDataPath.mat'), 'cohortPath');
catch
    error('No userDataPath.mat file found. Run the userDataPath.m file first.')
end
startPath = fullfile(cohortPath,'Analysis-Figures\Burstiness-Scatter');
analyzedStages = dir(startPath);
analyzedStages = analyzedStages([analyzedStages.isdir]);
analyzedStages = {analyzedStages(3:end).name};

[stagePick,tf] = listdlg('PromptString',{'Which stages do you want to analyze?'},...
    'SelectionMode','multiple','ListString',analyzedStages,'ListSize',[350 250]);

stageLabels = cell(1,numel(stagePick));
for i = 1:numel(stagePick)
    stageLabels(i) = inputdlg(sprintf('Rename each stage for figures. %i. selected stage:',i),'Stage labels',[1 100],analyzedStages(stagePick(i)));
end

% Change order of items if necessary 
f = figure('Name','Stage ordering');
set(gcf,'color','w')
axis off
text(0.5,0.9, {'Choose the order in which you want';'to visualize the stages.'},'HorizontalAlignment','center','FontSize',12,'Units','normalized')
hListbox = reorderableListbox('String', stageLabels,'FontSize',10,'Units','normalized','Position',[0.2 0.2 0.6 0.5]);

doneButton = uicontrol(f,'Style','pushbutton','units','normalized',...
    'Position',[0.8 0.1 0.12 0.07],'FontSize',10,'String','Done',...
    'Callback',{@doneExeStage,f,hListbox,stageLabels,stagePick});

waitfor(doneButton)

% Choose which waveforms to analyze
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
        'FontSize',10,'String',{'all','RS units','FS units'},'Value',initVal);
end

doneButton = uicontrol(f,'Style','pushbutton','units','normalized',...
    'Position',[0.8 0.1 0.1 0.05],'String','Done',...
    'Callback',{@doneExeWave,f,ui_field,area_names});

waitfor(doneButton)

% Define which cells to analyze
cellSpecs = listdlg('PromptString','Which cells would you like to analyze?', ...
    'Name','Cell types','ListString',{'Only biased cells (BCNs and TCNs)','All touch-modulated cells'}, ...
    'ListSize',[250,150],'SelectionMode','single');

%% Check if respective stages and trial types were already analyzed
% Create zero array with dimensions n_areas x n_stages x 2 (tonic, burst, touch-modulated)
filePresent = zeros(numel(area_waveforms),numel(stagePick),3);
% Return 1, if trial types were analyzed (only) in the right order,
% and 2, if trial types were analyzed (only) in reversed order (important for
% signed values).
% Return 3, if both sets of trial types are present
for i = 1:numel(stagePick)
    stageDir = fullfile(startPath,analyzedStages{stagePick(i)});
    stageDir = dir(stageDir);
    for ar = 1:numel(area_waveforms)
        for cellResp = 1:3
            % Check for waveform pre-filtered and un-filtered files
            if cellResp==1
                figFile = sprintf('ContourmapBurst-%s.fig',area_waveforms{ar});
            elseif cellResp==2
                figFile = sprintf('ContourmapTonic-%s.fig',area_waveforms{ar});
            elseif cellResp==3
                figFile = sprintf('ContourmapTouchModul-%s.fig',area_waveforms{ar});
            end
            preAnalysis = {stageDir((endsWith({stageDir.name},figFile) | ...
                endsWith({stageDir.name},figFile)) & ...
                contains({stageDir.name},condition)).name};
            if any(contains(preAnalysis,strjoin(trialType,'&'),'IgnoreCase',true)) && any(contains(preAnalysis,strjoin(flip(trialType),'&'),'IgnoreCase',true))
                filePresent(ar,i,cellResp) = 3;
            elseif any(contains(preAnalysis,strjoin(trialType,'&'),'IgnoreCase',true))
                filePresent(ar,i,cellResp) = 1;
            elseif any(contains(preAnalysis,strjoin(flip(trialType),'&'),'IgnoreCase',true))
                filePresent(ar,i,cellResp) = 2;
            end
        end
    end
end

% Adjust trial types according to analysis results
spikeTypes = {'Burst spikes (ContourmapBurst)','Tonic spikes (ContourmapTonic)','Touch-modulated cells (ContourmapBurst)'};
if any(~filePresent,'all') % Not all files analyzed
    [a,b,c] = ind2sub(size(filePresent),find(~filePresent));
    fprintf(2,'\nMissing files:\n')
    arrayfun(@(c,b,a) fprintf(2,'- %s for area "%s" in stage "%s"\n',spikeTypes{c},area_waveforms{a},analyzedStages{stagePick(b)}),c,b,a)
    return
elseif all(ismember(filePresent,[2,3])) % All files in reverse order
    fprintf('\nAll files were analyzed in reversed order. Continue with reversed order...\n\n')
    trialType = flip(trialType);
elseif ~all(ismember(filePresent,[1,3])) % Some files in reversed order, some not
    [a,b,c] = ind2sub(size(filePresent),find(filePresent==2));
    fprintf(2,'\nSome files were analyzed in reversed order (i.e., %s instead of %s), while others were not.\n',strjoin(flip(trialType),'&'),strjoin(trialType,'&'))
    fprintf(2,'Consider re-analyzing them, before you compare them.\n')
    fprintf(2,'Files in reversed order:\n')
    arrayfun(@(c,a,b) fprintf(2,'- %s for area "%s" in stage "%s"\n',spikeTypes{c},area_waveforms{a},analyzedStages{stagePick(b)}),c,a,b)
    return
end

%% Load contourmaps and save data points

% Bias values for each area and stage (cell array with n_area x n_stage)
burstBias = cell(numel(area_waveforms),numel(stagePick));
tonicBias = cell(numel(area_waveforms),numel(stagePick));
for i = 1:numel(stagePick)
    stageDir = fullfile(startPath,analyzedStages{stagePick(i)});
    stageDir = dir(stageDir);
    for ar = 1:numel(area_waveforms)
        % Retrieve burst bias data
        if cellSpecs==1
            figFile = sprintf('%s_%s_ContourmapBurst-%s.fig',condition,strjoin(trialType,'&'),area_waveforms{ar});
        elseif cellSpecs==2
            figFile = sprintf('%s_%s_ContourmapTouchModul-%s.fig',condition,strjoin(trialType,'&'),area_waveforms{ar});
        end
        
        fileIdx = find(strcmpi({stageDir.name},figFile));
        if isempty(fileIdx) && cellSpecs==1 % Un-filtered file search, if no file was found
            fileIdx = find(strcmpi({stageDir.name},sprintf('%s_%s_ContourmapBurst-%s.fig',condition,strjoin(trialType,'&'),area_names{ar})));
        elseif isempty(fileIdx) && cellSpecs==2 
            fileIdx = find(strcmpi({stageDir.name},sprintf('%s_%s_ContourmapTouchModul-%s.fig',condition,strjoin(trialType,'&'),area_names{ar})));
        end

        figName = fullfile(stageDir(fileIdx(1)).folder,stageDir(fileIdx(1)).name);
        fig = openfig(figName,'invisible');
        scat = findobj(fig,'Type','scatter');
        dataLabels = {scat.DataTipTemplate.DataTipRows.Label};
        
        % Save XData, as burst bias is encoded on the abscissa
        if endsWith(area_waveforms{ar},'RS')
            burstBias{ar,i} = scat.XData(scat.DataTipTemplate.DataTipRows(contains(dataLabels,'putExcit')).Value);
        elseif endsWith(area_waveforms{ar},'FS')
            burstBias{ar,i} = scat.XData(~scat.DataTipTemplate.DataTipRows(contains(dataLabels,'putExcit')).Value);
        else
            burstBias{ar,i} = scat.XData;
        end
        if cellSpecs==1
            close(fig)
        end

        % Retrieve tonic bias data
        if cellSpecs==1
            figFile = sprintf('%s_%s_ContourmapTonic-%s.fig',condition,strjoin(trialType,'&'),area_waveforms{ar});
            fileIdx = find(strcmpi({stageDir.name},figFile));
            if isempty(fileIdx) % Un-filtered file search, if no file was found
                fileIdx = find(strcmpi({stageDir.name},sprintf('%s_%s_ContourmapTonic-%s.fig',condition,strjoin(trialType,'&'),area_names{ar})));
            end
            
            figName = fullfile(stageDir(fileIdx(1)).folder,stageDir(fileIdx(1)).name);
            fig = openfig(figName,'invisible');
            scat = findobj(fig,'Type','scatter');
            dataLabels = {scat.DataTipTemplate.DataTipRows.Label};
        end
        
        % Save YData, as burst bias is encoded on the ordinate
        if endsWith(area_waveforms{ar},'RS')
            tonicBias{ar,i} = scat.YData(scat.DataTipTemplate.DataTipRows(contains(dataLabels,'putExcit')).Value);
        elseif endsWith(area_waveforms{ar},'FS')
            tonicBias{ar,i} = scat.YData(~scat.DataTipTemplate.DataTipRows(contains(dataLabels,'putExcit')).Value);
        else
            tonicBias{ar,i} = scat.YData;
        end
        close(fig)
    end
end

%% Plot data
% Generate tiled figure per area, in which overall number of significantly
% responding units are plotted, and units with respective weights.
for ar = 1:numel(area_waveforms)+1
    if ar==5 % All areas
        burstVals = cellfun(@(col) horzcat(burstBias{:, col}), num2cell(1:width(burstBias)), 'UniformOutput', false);
        tonicVals = cellfun(@(col) horzcat(tonicBias{:, col}), num2cell(1:width(tonicBias)), 'UniformOutput', false);
        plotColor = [0 0 0];
        if cellSpecs==1
            fig = figure('Name','CellBias_allAreas',...
                'Units','normalized','WindowState', 'maximized');
        elseif cellSpecs==2
            fig = figure('Name','CellBias_allAreas_allTouchModul',...
                'Units','normalized','WindowState', 'maximized');
        end
        tiledlayout(3,2)
        sgtitle('Aperture bias - All areas')
    else
        burstVals = burstBias(ar,:);
        tonicVals = tonicBias(ar,:);
        plotColor = hex2rgb(area_colors{ar});
        if cellSpecs==1
            fig = figure('Name',sprintf('CellBias_%s',area_waveforms{ar}),...
                'Units','normalized','WindowState', 'maximized');
        elseif cellSpecs==2
            fig = figure('Name',sprintf('CellBias_%s_allTouchModul',area_waveforms{ar}),...
                'Units','normalized','WindowState', 'maximized');
        end
        tiledlayout(3,2)
        sgtitle(sprintf('Aperture bias - %s',area_waveforms{ar}))
    end

    % Plot the number of cells with respective aperture bias
    % n_cells result from the number of signed values in each cell
    nexttile(1)
    n_cells = NaN(2,numel(stagePick));
    for i = 1:numel(stagePick)
        % Negative values = biased towards first condition
        n_cells(1,i) = -sum(burstVals{i}<0);
        % Positive values = biased towards second condition
        n_cells(2,i) = sum(burstVals{i}>0);
    end
    b = bar(n_cells', 'stacked','FaceColor',plotColor,'EdgeColor','None');
    b(1).FaceAlpha = 0.4;
    set(gca,'Layer', 'Top')
    title('Cells with burst bias')
    set(gca,'box','off')
    axis padded
    xticklabels(stageLabels)
    xtickangle(45)
    % Remove negative signs from y-ticks
    yticklabels(cellfun(@(x) strrep(x,'-',''), yticklabels,'UniformOutput',false))
    ylabel(sprintf('%s \\leftarrow Cells \\rightarrow %s',trialType{:}))

    nexttile(2)
    n_cells = NaN(2,numel(stagePick));
    for i = 1:numel(stagePick)
        % Negative values = biased towards first condition
        n_cells(1,i) = -sum(tonicVals{i}<0);
        % Positive values = biased towards second condition
        n_cells(2,i) = sum(tonicVals{i}>0);
    end
    b = bar(n_cells', 'stacked','FaceColor',plotColor,'EdgeColor','None');
    b(1).FaceAlpha = 0.4;
    set(gca,'Layer', 'Top')
    title('Cells with tonic bias')
    set(gca,'box','off')
    axis padded
    xticklabels(stageLabels)
    xtickangle(45)
    % Remove negative signs from y-ticks
    yticklabels(cellfun(@(x) strrep(x,'-',''), yticklabels,'UniformOutput',false))

    % Plot weighted number of cells with respective bias magnitude
    nexttile(3)
    weighted_cells = NaN(2,numel(stagePick));
    for i = 1:numel(stagePick)
        % Negative values = biased towards first condition
        cellIdx = burstVals{i}<0;
        weighted_cells(1,i) = sum(burstVals{i}(cellIdx));
        % Positive values = biased towards second condition
        cellIdx = burstVals{i}>0;
        weighted_cells(2,i) = sum(burstVals{i}(cellIdx));
    end
    b = bar(weighted_cells', 'stacked','FaceColor',plotColor,'EdgeColor','None');
    b(1).FaceAlpha = 0.4;
    set(gca,'Layer', 'Top')
    title('Weighted burst bias')
    set(gca,'box','off')
    axis padded
    xticklabels(stageLabels)
    xtickangle(45)
    % Remove negative signs from y-ticks
    yticklabels(cellfun(@(x) strrep(x,'-',''), yticklabels,'UniformOutput',false))
    ylabel(sprintf('%s \\leftarrow Cells*bias magnitude \\rightarrow %s',trialType{:}))

    nexttile(4)
    weighted_cells = NaN(2,numel(stagePick));
    for i = 1:numel(stagePick)
        % Negative values = biased towards first condition
        cellIdx = tonicVals{i}<0;
        weighted_cells(1,i) = sum(tonicVals{i}(cellIdx));
        % Positive values = biased towards second condition
        cellIdx = tonicVals{i}>0;
        weighted_cells(2,i) = sum(tonicVals{i}(cellIdx));
    end
    b = bar(weighted_cells', 'stacked','FaceColor',plotColor,'EdgeColor','None');
    b(1).FaceAlpha = 0.4;
    set(gca,'Layer', 'Top')
    title('Weighted tonic bias')
    set(gca,'box','off')
    axis padded
    xticklabels(stageLabels)
    xtickangle(45)
    % Remove negative signs from y-ticks
    yticklabels(cellfun(@(x) strrep(x,'-',''), yticklabels,'UniformOutput',false))

    % Plot bias values as violin plots
    nexttile(5)
    stageCats = cell2mat(cellfun(@(x,y) y*ones(numel(x),1), burstVals, num2cell(1:width(burstVals)), 'UniformOutput', false)');
    violinplot([burstVals{:}]',stageCats,'ShowMean',true,'ShowBox',...
        false,'ShowMedian',false,'ShowWhiskers',false,'ViolinAlpha',{[0.3,0.6]},...
        'MeanLineWidth',4,'ViolinColor',plotColor,'EdgeColor',plotColor);
    hold on
    line([0.5 width(burstVals)+0.5], [0 0], 'Color','k');
    xlim([0.5 width(burstVals)+0.5])
    ylims = ylim;
    pvals_burst = NaN(1,numel(stagePick));
    for i = 1:numel(stagePick)
        % Check if distribution differs significantly from 0
        if ~isempty(burstVals{i})
            pvals_burst(i) = signrank(burstVals{i});
        end

        if isempty(burstVals{i})
            text(i, 1.1*max(ylims), 'n/a','HorizontalAlignment','center','FontSize',14)
        elseif pvals_burst(i) < 0.001
            text(i, 1.1*max(ylims), '***','HorizontalAlignment','center','FontSize',14)
        elseif pvals_burst(i) < 0.01
            text(i, 1.1*max(ylims), '**','HorizontalAlignment','center','FontSize',14)
        elseif pvals_burst(i) < 0.05
            text(i, 1.1*max(ylims), '*','HorizontalAlignment','center','FontSize',14)
        else
            text(i, 1.1*max(ylims), 'n.s.','HorizontalAlignment','center','FontSize',14)
        end
    end
    ylim([ylims(1), 1.2*ylims(2)])
    title('Distribution of burst bias')
    set(gca,'box','off')
    xticks(1:numel(stagePick))
    xticklabels(stageLabels)
    xtickangle(45)
    % Remove negative signs from y-ticks
    yticklabels(cellfun(@(x) strrep(x,'-',''), yticklabels,'UniformOutput',false))
    ylabel(sprintf('%s \\leftarrow Bias \\rightarrow %s',trialType{:}))

    nexttile(6)
    stageCats = cell2mat(cellfun(@(x,y) y*ones(numel(x),1), tonicVals, num2cell(1:width(tonicVals)), 'UniformOutput', false)');
    violinplot([tonicVals{:}]',stageCats,'ShowMean',true,'ShowBox',...
        false,'ShowMedian',false,'ShowWhiskers',false,'ViolinAlpha',{[0.3,0.6]},...
        'MeanLineWidth',4,'ViolinColor',plotColor,'EdgeColor',plotColor);
    hold on
    l = line([0.5 width(tonicVals)+0.5], [0 0], 'Color','k');
    xlim([0.5 width(tonicVals)+0.5])
    ylims = ylim;
    pvals_tonic = NaN(1,numel(stagePick));
    for i = 1:numel(stagePick)
        % Check if distribution differs significantly from 0
        if ~isempty(tonicVals{i})
            pvals_tonic(i) = signrank(tonicVals{i});
        end

        if isempty(tonicVals{i})
            text(i, 1.1*max(ylims), 'n/a','HorizontalAlignment','center','FontSize',14)
        elseif pvals_tonic(i) < 0.001
            text(i, 1.1*max(ylims), '***','HorizontalAlignment','center','FontSize',14)
        elseif pvals_tonic(i) < 0.01
            text(i, 1.1*max(ylims), '**','HorizontalAlignment','center','FontSize',14)
        elseif pvals_tonic(i) < 0.05
            text(i, 1.1*max(ylims), '*','HorizontalAlignment','center','FontSize',14)
        else
            text(i, 1.1*max(ylims), 'n.s.','HorizontalAlignment','center','FontSize',14)
        end
    end
    ylim([ylims(1), 1.2*ylims(2)])
    title('Distribution of tonic bias')
    set(gca,'box','off')
    xticks(1:numel(stagePick))
    xticklabels(stageLabels)
    xtickangle(45)
    % Remove negative signs from y-ticks
    yticklabels(cellfun(@(x) strrep(x,'-',''), yticklabels,'UniformOutput',false))

    % Add information to figures
    fig.UserData = struct('fullStageNames', {analyzedStages(stagePick)},...
        'pvals_burst',pvals_burst, ...
        'pvals_tonic',pvals_tonic, ...
        'signTest','signrank');
end

%% Save figures
% Destination folder for matlab .fig files
definput = strjoin(analyzedStages(stagePick),'&');
stageDescription = inputdlg('Enter a description of the stages for saving the files:','Output Folder',[1 150],{definput});
destfold = fullfile(cohortPath,'Analysis-Figures','Burstiness-Bias',stageDescription{:});

if exist(destfold,"dir") == 0
    mkdir(destfold)
end

figHandles = handle(sort(double(findall(0, 'type', 'figure'))));

fprintf("Saving Figures.\n")

figure_suffix = [condition, '_', strjoin(trialType,'&')];
for i = 1:numel(figHandles)
    destfile = sprintf('%s\\%s_%s.fig', destfold, figure_suffix, figHandles(i).Name);
    savefig(figHandles(i), destfile);
end
fprintf("\nDone!\n\n")

%% Helper functions

% Function executed by the "doneButton" in stage ordering
function doneExeStage(~,~,f,hListbox,stageLabels,stagePick)

[~, changedOrder] = ismember(hListbox.String, stageLabels);
stagePick = stagePick(changedOrder);
stageLabels = stageLabels(changedOrder);

% Save variables
assignin('caller','stageLabels',stageLabels)
assignin('caller','stagePick',stagePick)

close(f)
end

function doneExeWave(~,~,f,ui_field,area_names)

area_waveforms = cell(1,numel(area_names));
for ar = 1:numel(area_names)
    switch ui_field(ar).Value
        case 2
            area_waveforms{ar} = [area_names{ar},'-RS'];
        case 3
            area_waveforms{ar} = [area_names{ar},'-FS'];
        otherwise
            area_waveforms{ar} = area_names{ar};
    end
end

% Overwrite variables
assignin('caller','area_waveforms',area_waveforms)

close(f)
end
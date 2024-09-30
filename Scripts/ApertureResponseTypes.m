%% Visualize individual responses upon stimulus presentation
close all; clearvars; clc;
% Make the user pick the two "ResponsePattern" mat files to compare
scriptFullPath = matlab.desktop.editor.getActiveFilename();
try load(regexprep(scriptFullPath, 'Scripts.*', 'Scripts\\userDataPath.mat'), 'cohortPath');
catch
    error('No userDataPath.mat file found. Run the userDataPath.m file first.')
end
startPath = fullfile(cohortPath,'Analysis-Figures\Burstiness-Analysis');

burstResponses = cell(1,2);
tonicResponses = cell(1,2);
trialTypeLabels = cell(1,2);
for trialType = 1:2
    [file,fileDir] = uigetfile('*.mat','Choose ResponsePattern mat file for each trial type.',startPath);
    startPath = fileDir;

    % Extracting respond window
    respondMatches = regexp(file, 'respond([-]?\d+\.\d+)-([-]?\d+\.\d+)', 'tokens');
    if trialType==1
        respWindow = str2double(respondMatches{1});
    else
        assert(isequal(respWindow,str2double(respondMatches{1})),'Respond windows of the chosen files do not match.')
    end
    % Extracting baseline window
    baselineMatches = regexp(file, 'baseline([-]?\d+\.\d+)-([-]?\d+\.\d+)', 'tokens');
    if trialType==1
        spontWindow = str2double(baselineMatches{1});
    else
        assert(isequal(spontWindow,str2double(baselineMatches{1})),'Baseline windows of the chosen files do not match.')
    end

    % Check if the ResponsePattern file is up to date
    vars = who('-file', fullfile(fileDir,file));
    if ismember('fileVersion', vars)
        % Load fileVersion if it exists
        load(fullfile(fileDir,file),'fileVersion')
        if fileVersion>=2
            load(fullfile(fileDir,file),'burstAPs_Info','tonicAPs_Info');
            burstResponses{trialType} = burstAPs_Info;
            tonicResponses{trialType} = tonicAPs_Info;
            clear('fileVersion')
        else
            fprintf(2,'\nResponsePattern file is not up to date. Re-run burstsUponTrigger.m script...\n')
            return
        end
    else
        fprintf(2,'\nResponsePattern file is not up to date. Re-run burstsUponTrigger.m script...\n')
        return
    end

    prompt = sprintf('Label %i: ',trialType);
    definput = regexp(file, '(?=narrow|wide)([^_]+)', 'match');
    trialTypeLabels{trialType} = cell2mat(inputdlg(prompt,'Trial type labels',[1 50],{strcat(upper(definput{:}(1)),lower(definput{:}(2:end)))}));

    condition = regexp(file, '_([a-zA-Z_]+)(?=narrow|wide)', 'match','once');
    condition = condition(2:end-1);

    stageDescript = cell2mat(regexp(fileparts(fileDir),'\\([^\\]*)$','tokens','once'));
end

% Destination folder for matlab .fig files
destfold = fullfile(cohortPath,'Analysis-Figures','Burstiness-Scatter',stageDescript);
if exist(destfold,"dir") == 0
    mkdir(destfold)
end
responseFileName = fullfile(destfold,sprintf('responseTypes_%s_%s.mat',condition,strjoin(trialTypeLabels,'&')));

% Save responseTypes.mat for all cells
runNsaveResponseTypes(burstResponses, tonicResponses, condition, spontWindow, respWindow, responseFileName)

%% Decide which areas to plot
area_names = {'BC','VPM','POm','ZIv'};
area_colors = {'#377eb8','#4daf4a','#984ea3','#ff7f00'};

[area_idx,tf] = listdlg('PromptString',{'Which area(s) do you want to analyze?'},...
    'SelectionMode','multiple','ListString',area_names,'InitialValue',1);
if tf
    for trialType = 1:2
        burstResponses{trialType} = burstResponses{trialType}(ismember(burstResponses{trialType}.Area,area_names(area_idx)),:);
        tonicResponses{trialType} = tonicResponses{trialType}(ismember(tonicResponses{trialType}.Area,area_names(area_idx)),:);
    end
else
    return
end

if numel(area_idx)==4
    area_descript = 'allAreas';
else
    area_descript = cell2mat(join(area_names(area_idx),'&'));
end

f = figure('Name','Waveform type');
set(gcf,'Position',[1000 600 420 300])
ui_field = gobjects(numel(area_idx),1);
ui_text = gobjects(numel(area_idx)+1,1);

ui_text(1) = uicontrol(f,'Style','text','Units','normalized',...
    'HorizontalAlignment','left','Position',[0.1 0.85 0.8 0.1],...
    'FontSize',10,'String','Choose the desired subpopulation for each area.');

for i = 1:numel(area_idx)
    ui_text(i+1) = uicontrol(f,'Style','text','Units','normalized',...
        'HorizontalAlignment','left','Position',[0.1 0.85-0.1*i 0.1 0.05],...
        'FontSize',10,'String',area_names{area_idx(i)});

    if isequal(area_names{area_idx(i)},'ZIv')
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
    'Callback',{@doneExe,f,ui_field,area_names,area_idx,...
    burstResponses,tonicResponses});

waitfor(doneButton)

animalNum = cellfun(@(x) str2double(regexp(x,'#(\d+)','tokens','once')), burstResponses{1}.SessionName, 'UniformOutput', false);
stageNum = cellfun(@(x) str2double(regexp(x,'P3\.(\d+\.\d+|\d+)','tokens','once')), burstResponses{1}.SessionName, 'UniformOutput', false);
sessionNum = cellfun(@(x) str2double(regexp(x,'ses\D*(\d+)','tokens','once')), burstResponses{1}.SessionName, 'UniformOutput', false);

%% Calculate relevant indices
% Calculate burst-index and burst-preference, e.g.,
% burstIndex [-1;1] = burstProportion_post - burstProportion_pre
% burstPreference [-1;1] = 0.5*(burstIndex_wide - burstIndex_narrow)
burstPreference = 0.5 * (cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)-x(:,1), burstResponses{2}.CellResponses, 'UniformOutput',false)) - ...
    cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)-x(:,1), burstResponses{1}.CellResponses, 'UniformOutput',false)));

% Perform a ranksum test on the results
burstSign = false(numel(burstPreference),1);
burstSign(~isnan(burstPreference)) = cellfun(@(y,z) ranksum(y,z), ...
    cellfun(@(x) x(:,2)-x(:,1), burstResponses{2}.CellResponses(~isnan(burstPreference)), 'UniformOutput',false),...
    cellfun(@(x) x(:,2)-x(:,1), burstResponses{1}.CellResponses(~isnan(burstPreference)), 'UniformOutput',false))<=0.05;

% Calculate tonic-index and tonic-preference, e.g.,
% tonicIndex [0;1] = tonicSpikes_post / (tonicSpikes_post + tonicSpikes_pre)
% tonicPreference [-1;1] = tonicIndex_wide - tonicIndex_narrow
tonicPreference = cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)./(x(:,2)+x(:,1)), tonicResponses{2}.CellResponses, 'UniformOutput',false)) - ...
    cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)./(x(:,2)+x(:,1)), tonicResponses{1}.CellResponses, 'UniformOutput',false));

% Perform a ranksum test on the results
tonicSign = false(numel(tonicPreference),1);
tonicSign(~isnan(tonicPreference)) = cellfun(@(y,z) ranksum(y,z), ...
    cellfun(@(x) x(:,2)-x(:,1), tonicResponses{2}.CellResponses(~isnan(tonicPreference)), 'UniformOutput',false),...
    cellfun(@(x) x(:,2)-x(:,1), tonicResponses{1}.CellResponses(~isnan(tonicPreference)), 'UniformOutput',false))<=0.05;

% For retrieving the touch-modulated cells
spontWindow = [-600,-400]; % spontaneous window in msec (default: [-600,-400])
respWindow = [0,200]; % response window in msec (default: [0,200])
singleCellReactivityName = ['SingleCellReactivity_', condition, '_withTargetEstimate_', num2str(respWindow/1000,'respond%.3f-%.3f_'), num2str(spontWindow/1000,'baseline%.3f-%.3f.mat')];

fileSelection = unique(burstResponses{1}.SessionName);
touchModulated = false(numel(burstPreference),1);
for file = 1:height(fileSelection)
    load(fullfile(fileSelection{file},singleCellReactivityName),'SingleCellResponse')
    % Match the selected units from the tonic and burst responses with
    % the respective ones from the SingleCellResponse file
    for i = find(ismember(burstResponses{1}.SessionName,fileSelection{file}))'
        idx = find(SingleCellResponse.cluster_id==str2double(burstResponses{1}.UnitID{i}));
        touchModulated(i) = SingleCellResponse.respWide(idx) | SingleCellResponse.respNarrow(idx) | SingleCellResponse.respInterm(idx);
    end
end

%% Scatter results
fig = figure('Name',sprintf('%s_BurstTonicScatter',area_descript));
ax1 = gca;
hold on
% grouping = categorical(burstResponses{1}.Area);
% scatterhist(burstPreference,tonicPreference,'Group',grouping,'NBins',40,'Color',hex2rgb(area_colors(area_idx)));

scat0 = gobjects(1,numel(area_idx));
scat1 = gobjects(1,numel(area_idx));
scat2 = gobjects(1,numel(area_idx));
scat3 = gobjects(1,numel(area_idx));
signProportions = NaN(numel(area_idx),4);
for i = 1:numel(area_idx)
    idx = ismember(burstResponses{1}.Area,area_names{area_idx(i)});
    % Non-responders
    scat0(i) = scatter(ax1, burstPreference(idx & ~burstSign & ~tonicSign),tonicPreference(idx & ~burstSign & ~tonicSign),'.','MarkerEdgeColor',area_colors{i});
    signProportions(i,1) = sum(idx & ~burstSign & ~tonicSign);
    scat0(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unitID',burstResponses{1}.UnitID(idx & ~burstSign & ~tonicSign));
    scat0(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('animalNum',animalNum(idx & ~burstSign & ~tonicSign));
    scat0(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('stageNum',stageNum(idx & ~burstSign & ~tonicSign));
    scat0(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('sessionNum',sessionNum(idx & ~burstSign & ~tonicSign));
    scat0(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('putExcit',burstResponses{1}.PutExcitatory(idx & ~burstSign & ~tonicSign));

    % Only significant burst preference
    scat1(i) = scatter(ax1, burstPreference(idx & burstSign & ~tonicSign),tonicPreference(idx & burstSign & ~tonicSign),'o','filled','MarkerFaceColor',area_colors{i});
    signProportions(i,2) = sum(idx & burstSign & ~tonicSign);
    scat1(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unitID',burstResponses{1}.UnitID(idx & burstSign & ~tonicSign));
    scat1(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('animalNum',animalNum(idx & burstSign & ~tonicSign));
    scat1(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('stageNum',stageNum(idx & burstSign & ~tonicSign));
    scat1(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('sessionNum',sessionNum(idx & burstSign & ~tonicSign));
    scat1(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('putExcit',burstResponses{1}.PutExcitatory(idx & burstSign & ~tonicSign));

    % Only significant tonic preference
    scat2(i) = scatter(ax1, burstPreference(idx & ~burstSign & tonicSign),tonicPreference(idx & ~burstSign & tonicSign),'square','filled','MarkerFaceColor',area_colors{i});
    signProportions(i,3) = sum(idx & ~burstSign & tonicSign);
    scat2(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unitID',burstResponses{1}.UnitID(idx & ~burstSign & tonicSign));
    scat2(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('animalNum',animalNum(idx & ~burstSign & tonicSign));
    scat2(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('stageNum',stageNum(idx & ~burstSign & tonicSign));
    scat2(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('sessionNum',sessionNum(idx & ~burstSign & tonicSign));
    scat2(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('putExcit',burstResponses{1}.PutExcitatory(idx & ~burstSign & tonicSign));

    % Both significant burst and tonic preference
    scat3(i) = scatter(ax1, burstPreference(idx & burstSign & tonicSign),tonicPreference(idx & burstSign & tonicSign),'pentagram','filled','MarkerFaceColor',area_colors{i});
    signProportions(i,4) = sum(idx & burstSign & tonicSign);
    scat3(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unitID',burstResponses{1}.UnitID(idx & burstSign & tonicSign));
    scat3(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('animalNum',animalNum(idx & burstSign & tonicSign));
    scat3(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('stageNum',stageNum(idx & burstSign & tonicSign));
    scat3(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('sessionNum',sessionNum(idx & burstSign & tonicSign));
    scat3(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('putExcit',burstResponses{1}.PutExcitatory(idx & burstSign & tonicSign));
end

xlims = xlim;
ylims = ylim;

l = line([0 0], [-1.1 1.1], 'Color','k');
uistack(l,'bottom')

l = line([-1.1 1.1],[0 0], 'Color','k');
uistack(l,'bottom')

xlabel(sprintf('%s \\leftarrow Burst preference \\rightarrow %s',trialTypeLabels{:}))
ylabel(sprintf('%s \\leftarrow Tonic preference \\rightarrow %s',trialTypeLabels{:}))

if numel(area_idx)>1
    % Copy markers for the first legend
    scatCopy = copyobj(scat1, ax1);
    set(scatCopy,'XData', NaN', 'YData', NaN)
    legend1 = legend(ax1,scatCopy,area_names(area_idx),'Location','northeast');
end

% Cumulative cell response proportions
axes('Position',[.15 .15 .15 .15]), box off
piePlot = pie(sum(signProportions,1),{'none','burst','tonic','both'});
patchHand = findobj(piePlot, 'Type', 'Patch');
patchHand(1).FaceColor = 0.8*[1 1 1];
patchHand(2).FaceColor = 0.6*[1 1 1];
patchHand(3).FaceColor = 0.4*[1 1 1];
patchHand(4).FaceColor = 0.2*[1 1 1];

% Copy the axes and plot the second legend
ax2 = axes('position',get(ax1,'position'),'visible','off');
legend2 = legend(ax2, [scat0(1), scat1(1),scat2(1),scat3(1)],{'none','sign burst','sign tonic','sign burst & tonic'},'Location','northwest');
linkaxes([ax1 ax2],'xy')

xlim([-max(abs(xlims))*1.1 max(abs(xlims))*1.1])
ylim([-max(abs(ylims))*1.1 max(abs(ylims))*1.1])

definput = sprintf('Burst and tonic response profile of %s',area_descript);
userTitle = inputdlg('Enter a figure title:','Figure title',[1 100],{definput});
title(ax1,userTitle)

fig.UserData = struct('Areas',{area_names(area_idx)}, ...
    'waveformTypes', {waveformTypes},...
    'groupCat', {{'none','burst','tonic','both'}},...
    'signProportions',sum(signProportions,1), ...
    'cellNum',sum(signProportions,2));

% Highlight FS and RS in different shades

% Define subgroups, either as their position in the diagram (i.e., the quadrants),
% or by a grouping algorithm (e.g., k-means clustering)

%% Plot heat map of each area

spikeDescript = {'burst responders','tonic responders', 'touch-modul cells'};
if numel(area_idx)==4 % If all areas, plot also a pooled analysis
    loopRuns = 5;
else
    loopRuns = numel(area_idx);
end
for i = 1:loopRuns
    for spikeType = 1:3
        if i==5 % Plot pooled data
            idx = true(height(burstResponses{1}.Area),1);
        else
            idx = ismember(burstResponses{1}.Area,area_names{area_idx(i)});
        end
        % Compute the 2D histogram
        stepSize = 0.02;
        Xedges = (-1:stepSize:1);
        Yedges = (-1:stepSize:1);
        if spikeType==1
            x = burstPreference(idx & burstSign);
            y = tonicPreference(idx & burstSign);
        elseif spikeType==2
            x = burstPreference(idx & tonicSign);
            y = tonicPreference(idx & tonicSign);
        elseif spikeType==3
            x = burstPreference(idx & touchModulated);
            y = tonicPreference(idx & touchModulated);
        end
        histogram2D = histcounts2(x,y, Xedges, Yedges);

        % Define the structuring element for dilation
        se = strel('disk', 1);  % Adjust the disk size as needed

        % Smooth the matrix over avg_bins x avg_bins
        avg_bins = 5;
        avg_matrix = ones(avg_bins,avg_bins*2)/avg_bins^2*2;
        histogram2D = conv2(histogram2D,avg_matrix,'same');

        % Plot the dilated image as a heat map
        if i==5 % Plot pooled data
            figName = 'allAreas';
            areaDescript = [sprintf('%s, ',area_names_waveforms{1:end-1}),sprintf('%s',area_names_waveforms{end})];
        else
            figName = area_names_waveforms{area_idx(i)};
            areaDescript = area_names_waveforms{area_idx(i)};
        end

        if spikeType==1
            fig = figure('Name',sprintf('HeatmapBurst-%s',figName));
            fig.UserData = struct('avg_bins', avg_bins,...
                'cellNum',sum(idx & burstSign), ...
                'areaDescript', areaDescript);
        elseif spikeType==2
            fig = figure('Name',sprintf('HeatmapTonic-%s',figName));
            fig.UserData = struct('avg_bins', avg_bins,...
                'cellNum',sum(idx & tonicSign), ...
                'areaDescript', areaDescript);
        elseif spikeType==3
            fig = figure('Name',sprintf('HeatmapTouchModul-%s',figName));
            fig.UserData = struct('avg_bins', avg_bins,...
                'cellNum',sum(idx & touchModulated), ...
                'areaDescript', areaDescript);            
        end
        hold on
        %     s = scatter(x,y,'o','MarkerEdgeColor',area_colors{area_idx(i)});
        imvals_X = (-1+stepSize/2:stepSize:1-stepSize/2);
        imvals_Y = (-1+stepSize/2:stepSize:1-stepSize/2);
        im = imagesc(imvals_X, imvals_Y, histogram2D');

        if i==5 % Plot pooled data
            cm = brewermap([],"Greys");
        else
            switch area_names{area_idx(i)}
                case 'BC'
                    cm = brewermap([],"Blues");
                case 'VPM'
                    cm = brewermap([],"BuGn");
                case 'POm'
                    cm = brewermap([],"BuPu");
                case 'ZIv'
                    cm = brewermap([],"Oranges");
            end
        end
        % Switch first value to one to make background appear white
        cm(1,:) = [1 1 1];
        colormap(cm)

        % Add grid lines
        line([0 0], [-1.1 1.1], 'Color','k');
        line([-1.1 1.1],[0 0], 'Color','k');

        xlim([-max(abs(xlims))*1.1 max(abs(xlims))*1.1])
        ylim([-max(abs(ylims))*1.1 max(abs(ylims))*1.1])

        set(gca,'Layer', 'Top')
        xlabel(sprintf('%s \\leftarrow Burst preference \\rightarrow %s',trialTypeLabels{:}))
        ylabel(sprintf('%s \\leftarrow Tonic preference \\rightarrow %s',trialTypeLabels{:}))
        if i==5 % Plot pooled data
            title(sprintf('Heat map of %s - All areas',spikeDescript{spikeType}));
        else
            title(sprintf('Heat map of %s - %s %s',spikeDescript{spikeType},area_names{area_idx(i)},waveformTypes{i}));
        end
    end
end

%% Plot first distribution contours of each area

for i = 1:loopRuns
    for spikeType = 1:3
        if i==5 % Plot pooled data
            idx = true(height(burstResponses{1}.Area),1);
        else
            idx = ismember(burstResponses{1}.Area,area_names{area_idx(i)});
        end
        % Compute the 2D histogram
        stepSize = 0.02;
        Xedges = (-1:stepSize:1);
        Yedges = (-1:stepSize:1);
        if spikeType==1
            x = burstPreference(idx & burstSign);
            y = tonicPreference(idx & burstSign);
        elseif spikeType==2
            x = burstPreference(idx & tonicSign);
            y = tonicPreference(idx & tonicSign);
        elseif spikeType==3
            x = burstPreference(idx & touchModulated);
            y = tonicPreference(idx & touchModulated);
        end

        numLevels = 3;

        if i==5 % Plot pooled data
            figName = 'allAreas';
            areaDescript = [sprintf('%s, ',area_names_waveforms{1:end-1}),sprintf('%s',area_names_waveforms{end})];
        else
            figName = area_names_waveforms{area_idx(i)};
            areaDescript = area_names_waveforms{area_idx(i)};
        end

        if spikeType==1
            fig = figure('Name',sprintf('ContourmapBurst-%s',figName));
            fig.UserData = struct('numLevels', numLevels,...
                'cellNum',sum(idx & burstSign), ...
                'areaDescript', areaDescript);
        elseif spikeType==2
            fig = figure('Name',sprintf('ContourmapTonic-%s',figName));
            fig.UserData = struct('numLevels', numLevels,...
                'cellNum',sum(idx & tonicSign), ...
                'areaDescript', areaDescript);
        elseif spikeType==3
            fig = figure('Name',sprintf('ContourmapTouchModul-%s',figName));
            fig.UserData = struct('numLevels', numLevels,...
                'cellNum',sum(idx & touchModulated), ...
                'areaDescript', areaDescript);
        end

        if i==5 % Plot pooled data
            plotColor = [0 0 0];
        else
            plotColor = hex2rgb(area_colors{area_idx(i)});
        end

        if numel(x)<=1
            scat = scatter(x,y, 'MarkerFaceColor', plotColor, 'MarkerEdgeColor','None');
        else
            scat = kscontour([x,y], 'Color', plotColor, 'ShowPoints', true, 'Nlevels', numLevels);
        end

        % Add additional info to data points
        if spikeType==1
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unitID',burstResponses{1}.UnitID(idx & burstSign));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('animalNum',animalNum(idx & burstSign));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('stageNum',stageNum(idx & burstSign));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('sessionNum',sessionNum(idx & burstSign));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('putExcit',burstResponses{1}.PutExcitatory(idx & burstSign));
        elseif spikeType==2
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unitID',burstResponses{1}.UnitID(idx & tonicSign));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('animalNum',animalNum(idx & tonicSign));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('stageNum',stageNum(idx & tonicSign));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('sessionNum',sessionNum(idx & tonicSign));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('putExcit',burstResponses{1}.PutExcitatory(idx & tonicSign));
        elseif spikeType==3
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('unitID',burstResponses{1}.UnitID(idx & touchModulated));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('animalNum',animalNum(idx & touchModulated));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('stageNum',stageNum(idx & touchModulated));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('sessionNum',sessionNum(idx & touchModulated));
            scat.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('putExcit',burstResponses{1}.PutExcitatory(idx & touchModulated));
        end

        % Add grid lines
        line([0 0], [-1.1 1.1], 'Color','k');
        line([-1.1 1.1],[0 0], 'Color','k');

        xlim([-max(abs(xlims))*1.1 max(abs(xlims))*1.1])
        ylim([-max(abs(ylims))*1.1 max(abs(ylims))*1.1])

        set(gca,'Layer', 'Top')
        xlabel(sprintf('%s \\leftarrow Burst preference \\rightarrow %s',trialTypeLabels{:}))
        ylabel(sprintf('%s \\leftarrow Tonic preference \\rightarrow %s',trialTypeLabels{:}))
        if i==5 % Plot pooled data
            title(sprintf('Contour lines of %s - All areas',spikeDescript{spikeType}));
        else
            title(sprintf('Contour lines of %s - %s %s',spikeDescript{spikeType},area_names{area_idx(i)},waveformTypes{i}))
        end

    end
end

%% Save figures
figHandles = handle(sort(double(findall(0, 'type', 'figure'))));

fprintf("Saving Figures.\n")

figure_suffix = [condition, '_', strjoin(trialTypeLabels,'&')];
for i = 1:numel(figHandles)
    destfile = sprintf('%s\\%s_%s.fig', destfold, figure_suffix, figHandles(i).Name);
    savefig(figHandles(i), destfile);
end
fprintf("\nDone!\n\n")


%% Helper functions
% Function executed by the "doneButton" in subpopulation refinement
function doneExe(~,~,f,ui_field,area_names,area_idx,...
    burstResponses,tonicResponses)

idx = false(height(burstResponses{1}),1);
waveformTypes = cell(1,numel(area_idx));
area_names_waveforms = cell(1,numel(area_idx));
for i = 1:numel(area_idx)
    switch ui_field(i).Value
        case 1 % 'all'
            waveformTypes{i} = 'all units';
            area_names_waveforms{i} = area_names{i};
            idx = idx | contains(burstResponses{1}.Area,area_names{area_idx(i)});
        case 2 % 'RS units'
            waveformTypes{i} = 'RS units';
            area_names_waveforms{i} = [area_names{i},'-RS'];
            idx = idx | (contains(burstResponses{1}.Area,area_names{area_idx(i)}) & ...
                burstResponses{1}.PutExcitatory);
        case 3 % 'FS units'
            waveformTypes{i} = 'FS units';
            area_names_waveforms{i} = [area_names{i},'-FS'];
            idx = idx | (contains(burstResponses{1}.Area,area_names{area_idx(i)}) & ...
                ~burstResponses{1}.PutExcitatory);
    end
end

for trialType = 1:2
    burstResponses{trialType} = burstResponses{trialType}(idx,:);
    tonicResponses{trialType} = tonicResponses{trialType}(idx,:);
end

% Overwrite variables
assignin('caller','burstResponses',burstResponses)
assignin('caller','tonicResponses',tonicResponses)
assignin('caller','waveformTypes',waveformTypes)
assignin('caller','area_names_waveforms',area_names_waveforms)

close(f)
end

function runNsaveResponseTypes(burstResponses, tonicResponses, condition, spontWindow, respWindow, responseFileName)
% Calculate burst-index and burst-preference, e.g.,
% burstIndex [-1;1] = burstProportion_post - burstProportion_pre
% burstPreference [-1;1] = 0.5*(burstIndex_wide - burstIndex_narrow)
burstPreference = 0.5 * (cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)-x(:,1), burstResponses{2}.CellResponses, 'UniformOutput',false)) - ...
    cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)-x(:,1), burstResponses{1}.CellResponses, 'UniformOutput',false)));

% Perform a ranksum test on the results
burstSign = false(numel(burstPreference),1);
burstSign(~isnan(burstPreference)) = cellfun(@(y,z) ranksum(y,z), ...
    cellfun(@(x) x(:,2)-x(:,1), burstResponses{2}.CellResponses(~isnan(burstPreference)), 'UniformOutput',false),...
    cellfun(@(x) x(:,2)-x(:,1), burstResponses{1}.CellResponses(~isnan(burstPreference)), 'UniformOutput',false))<=0.05;

% Calculate tonic-index and tonic-preference, e.g.,
% tonicIndex [0;1] = tonicSpikes_post / (tonicSpikes_post + tonicSpikes_pre)
% tonicPreference [-1;1] = tonicIndex_wide - tonicIndex_narrow
tonicPreference = cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)./(x(:,2)+x(:,1)), tonicResponses{2}.CellResponses, 'UniformOutput',false)) - ...
    cellfun(@(y) mean(y(isfinite(y))), cellfun(@(x) x(:,2)./(x(:,2)+x(:,1)), tonicResponses{1}.CellResponses, 'UniformOutput',false));

% Perform a ranksum test on the results
tonicSign = false(numel(tonicPreference),1);
tonicSign(~isnan(tonicPreference)) = cellfun(@(y,z) ranksum(y,z), ...
    cellfun(@(x) x(:,2)-x(:,1), tonicResponses{2}.CellResponses(~isnan(tonicPreference)), 'UniformOutput',false),...
    cellfun(@(x) x(:,2)-x(:,1), tonicResponses{1}.CellResponses(~isnan(tonicPreference)), 'UniformOutput',false))<=0.05;

% For retrieving the touch-modulated cells
singleCellReactivityName = ['SingleCellReactivity_', condition, '_withTargetEstimate_', num2str(respWindow/1000,'respond%.3f-%.3f_'), num2str(spontWindow/1000,'baseline%.3f-%.3f.mat')];

fileSelection = unique(burstResponses{1}.SessionName);
touchModulated = false(numel(burstPreference),1);
for file = 1:height(fileSelection)
    load(fullfile(fileSelection{file},singleCellReactivityName),'SingleCellResponse')
    % Match the selected units from the tonic and burst responses with
    % the respective ones from the SingleCellResponse file
    for i = find(ismember(burstResponses{1}.SessionName,fileSelection{file}))'
        idx = find(SingleCellResponse.cluster_id==str2double(burstResponses{1}.UnitID{i}));
        touchModulated(i) = SingleCellResponse.respWide(idx) | SingleCellResponse.respNarrow(idx) | SingleCellResponse.respInterm(idx);
    end
end

animalNum = cellfun(@(x) str2double(regexp(x,'#(\d+)','tokens','once')), burstResponses{1}.SessionName, 'UniformOutput', false);
stageNum = cellfun(@(x) str2double(regexp(x,'P3\.(\d+\.\d+|\d+)','tokens','once')), burstResponses{1}.SessionName, 'UniformOutput', false);
sessionNum = cellfun(@(x) str2double(regexp(x,'ses\D*(\d+)','tokens','once')), burstResponses{1}.SessionName, 'UniformOutput', false);

% Save the response types as an individual variable.
responseTypes = table(burstResponses{1}.SessionName, burstResponses{1}.Area, animalNum, stageNum, sessionNum, ...
    burstResponses{1}.UnitID, burstResponses{1}.PutExcitatory,burstSign,tonicSign,touchModulated, ...
    'VariableNames',{'SessionName','Area','Animal','Stage','Session','UnitID','PutExcitatory','BurstSign','TonicSign','TouchModul'});
save(responseFileName,'responseTypes');

end
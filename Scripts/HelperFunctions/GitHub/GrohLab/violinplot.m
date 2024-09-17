%% Firing Rates Violin Plots - Separate Plots per Group

mechResponse = clInfo.Mech_Control_15mW_R;
ruNames =   {'VPL'}; % [{'Group 1'}, {'Group 2'}, {'L6'}];


red = [0.75, 0, 0];
green = [0, 0.75, 0];
blue = [0.25, 0.5, 1];
purple = [0.5,0,0.5];
% colours = [green; red; blue];
colours = [0,0,0];

%  idxLaterTbl = ismember(clInfo.id, Later);
%  idxLatestTbl = ismember(clInfo.id, Latest);
%  idxL6Tbl = ismember(clInfo.id, L6);
%  idxMat = [idxLaterTbl, idxLatestTbl, idxL6Tbl]  & mechResponse;
idxMat = mechResponse; % & clInfo.Mech_Laser_1Hz_15mW_Rate_Evoked < 60;
matDims = size(idxMat);
nGroups = matDims(2);

% yMax = [max(clInfo.Laser_5sec_pulse_15mW_Block_Rate_Spont), max(clInfo.Laser_5sec_pulse_15mW_Block_Rate_Evoked)];
yMax = round(max(yMax)+5, -1);
yMin = 0;

for unitGroup = 1:nGroups
    
    
    unitIds = idxMat(:,unitGroup)==true;
    nUnits = sum(unitIds);
    rts = NaN(nUnits, 2);
    rts(:,1) = clInfo.Mech_Control_15mW_Rate_Evoked(unitIds);
    rts(:,2) = clInfo.Mech_Laser_1Hz_15mW_Rate_Evoked(unitIds);
    [p, sig] = ranksum(rts(:,1), rts(:,2));
    fprintf([ruNames{unitGroup}, ' p value: ', num2str(p), '\n']);
    figure('Color', 'White', 'Name', [ruNames{unitGroup}, '_Evoked_Rates_ViolinPlot']);
    hold on
    
    
    [yCtr, xCtr] = ksdensity(rts(:,1),'Bandwidth',0.7);
    patch(yCtr * -1,xCtr,colours(unitGroup,:),'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
    [yCtr, xCtr] = ksdensity(rts(:,2),'Bandwidth',0.7);
    patch(yCtr,xCtr,colours(unitGroup,:),'EdgeAlpha',0.1,'FaceAlpha',0.5);
    
    
    
%     plot(zeros(nUnits,1) - 0.0025, rts(:,1),...
%         'LineStyle','none', 'Marker', '*', 'Color', colour(unitGroup,:), 'MarkerSize', 2);
%     
%     plot(zeros(nUnits,1) + 0.0025, rts(:,2),...
%         'LineStyle','none', 'Marker', '*', 'Color', [0,0,1], 'MarkerSize', 2);
    
    medControl = median(rts(:,1));
    medLaser = median(rts(:,2));
       
    
%     leg = legend('Mech Control','Mech Laser', 'Location','northwest');
%     leg.Box = 'off';
%     leg.Location = 'northeast';
    
    
    
    ax = gca;
    ax.Title.String = ruNames{unitGroup};
    ax.FontName = 'Arial';
    ax.FontSize = 10;
    xscale = max(abs(ax.XAxis.Limits));
    yMax = ax.YAxis.Limits(2);
    yMax = round(max(yMax)+5, -1);
    
    
    plot([-xscale/5,0], [medControl, medControl], 'LineStyle', '-', 'LineWidth', 1, 'Color', colours(unitGroup,:));
    plot([0,xscale/5], [medLaser, medLaser], 'LineStyle', '-', 'LineWidth', 1, 'Color', colours(unitGroup,:));
    
    plot(zeros(nUnits,1) - xscale/50, rts(:,1),...
        'LineStyle','none', 'Marker', '*', 'Color', colours(unitGroup,:), 'MarkerSize', 2);
    
    plot(zeros(nUnits,1) + xscale/50, rts(:,2),...
        'LineStyle','none', 'Marker', '*', 'Color', colours(unitGroup,:), 'MarkerSize', 2);
    
    
    %leg.String = [{'Mech Control'}, {'Mech + Laser'}, {'Population Median'}];
    ax.XAxis.Limits = [-xscale, xscale];
    ax.YAxis.Limits = [yMin, yMax];
    ax.YLabel.String = 'Firing Rate [Hz]';
    ax.XTick = [];
    %     ax.Visible = 'off';
end

%% Firing Rates Violin Plots - One Plot

% MechResponse = clInfo.Mech_Control_15mW_MR;
ruNames = [{'Group 1'}, {'Group 2'}, {'L6'}];


red = [0.75, 0, 0];
green = [0, 0.75, 0];
blue = [0.25, 0.5, 1];
purple = [0.5,0,0.5];
colours = [red; green; blue];


idxLaterTbl = ismember(clInfo.id, Later);
idxLatestTbl = ismember(clInfo.id, Latest);
idxL6Tbl = ismember(clInfo.id, L6);
idxMat = [idxLaterTbl, idxLatestTbl, idxL6Tbl]; % & MechResponse;
matDims = size(idxMat);
nGroups = matDims(2);

yMax = [max(clInfo.Mech_Control_4mW_Rate_Evoked), max(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked)];
yMax = round(max(yMax)+5, -1);
yMin = 0;
 figure('Color', 'White', 'Name','Evoked_Rates_ViolinPlot');
for unitGroup = 1:nGroups
    
    
    unitIds = idxMat(:,unitGroup)==true;
    nUnits = sum(unitIds);
    rts = NaN(nUnits, 2);
    rts(:,1) = clInfo.Laser_5sec_pulse_4mW_Block_Rate_Spont(unitIds);
    rts(:,2) = clInfo.Laser_5sec_pulse_4mW_Block_Rate_Evoked(unitIds);
    
    subplot(nGroups,1,unitGroup)
    hold on
    
    
    [yCtr, xCtr] = ksdensity(rts(:,1),'Bandwidth',0.7);
    patch(yCtr * -1,xCtr,colour(unitGroup,:),'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
    [yCtr, xCtr] = ksdensity(rts(:,2),'Bandwidth',0.7);
    patch(yCtr,xCtr,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
    
    
%     plot(zeros(nUnits,1) - 0.0025, rts(:,1),...
%         'LineStyle','none', 'Marker', '*', 'Color', colour(unitGroup,:), 'MarkerSize', 2);
%     
%     plot(zeros(nUnits,1) + 0.0025, rts(:,2),...
%         'LineStyle','none', 'Marker', '*', 'Color', [0,0,1], 'MarkerSize', 2);
    
    medControl = median(rts(:,1));
    medLaser = median(rts(:,2));
       
    
    leg = legend('Mech Control','Mech Laser', 'Location','northwest');
    leg.Box = 'off';
    leg.Location = 'northeast';
    
    
    
    ax = gca;
    ax.Title.String = ruNames{unitGroup};
    ax.FontName = 'Arial';
    ax.FontSize = 10;
    xscale = max(abs(ax.XAxis.Limits));
    
    
    
    plot([-xscale/10,0], [medControl, medControl], 'LineStyle', '-', 'LineWidth', 1, 'Color', colour(unitGroup,:));
    plot([0,xscale/10], [medLaser, medLaser], 'LineStyle', '-', 'LineWidth', 1, 'Color', [0,0,1]);
    
    plot(zeros(nUnits,1) - xscale/50, rts(:,1),...
        'LineStyle','none', 'Marker', '*', 'Color', colour(unitGroup,:), 'MarkerSize', 2);
    
    plot(zeros(nUnits,1) + xscale/50, rts(:,2),...
        'LineStyle','none', 'Marker', '*', 'Color', [0,0,1], 'MarkerSize', 2);
    
    
    leg.String = [{'Mech Control'}, {'Mech + Laser'}, {'Population Median'}];
    ax.XAxis.Limits = [-xscale, xscale];
    %ax.YAxis.Limits = [yMin, yMax];
    ax.YAxis.Limits(1) = 0;
    ax.YLabel.String = 'Firing Rate [Hz]';
    ax.XTick = [];
    %     ax.Visible = 'off';
end

%% Relative Rates Violin Plots

mechResponse = clInfo.Mech_Control_4mW_R;
ids = [{'L6'}, {'S1'}, {'Vpl'}];
idxMat = [idxTagged, idxNonTagged]; % & MechResponse;
matDims = size(idxMat);
figure('Name', 'Relative Firing Rate Violin Plots', 'Color', 'White');
nGroups = matDims(2);

yMax = [max(clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Control_4mW_Rate_Spont_1_to_2s),...
    max(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_1_to_2s)];
yMax = max(yMax);
yMin = [min(clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Control_4mW_Rate_Spont_1_to_2s),...
    min(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_1_to_2s)];
yMin = min(yMin);


for unitGroup = 1:nGroups
    subplot(1, nGroups, unitGroup)
    
    unitIds = idxMat(:,unitGroup);
    nUnits = sum(unitIds);
    rts = NaN(nUnits, 2);
    rts(:,1) = clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s(unitIds) - clInfo.Mech_Control_4mW_Rate_Spont_1_to_2s(unitIds);
    rts(:,2) = clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s(unitIds) - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_1_to_2s(unitIds);
    
    subplot(1,2,unitGroup)
    hold on
    
    
    [yCtr, xCtr] = ksdensity(rts(:,1),'Bandwidth',0.7);
    patch(yCtr * -1,xCtr,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
    [yCtr, xCtr] = ksdensity(rts(:,2),'Bandwidth',0.7);
    patch(yCtr,xCtr,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
    
    
    plot(zeros(nUnits,1) - 0.01, rts(:,1),...
        'LineStyle','none', 'Marker', '*', 'Color', [1,0,0], 'MarkerSize', 5);
    
    plot(zeros(nUnits,1) + 0.01, rts(:,2),...
        'LineStyle','none', 'Marker', '*', 'Color', [0,0,1], 'MarkerSize', 5);
    
    
    
    legend('Mech Control','Mech Laser', 'Location','northwest');
    
    ax = gca;
    ax.Title.String = ids{unitGroup};
    ax.FontName = 'Arial';
    ax.FontSize = 25;
    ax.YAxis.Limits = [yMin, yMax];
    
end
%% SNR Violin Plots


mechResponse = clInfo.Mech_Control_4mW_R;
ids = [{'L6'}, {'S1'}, {'Vpl'}];
idxMat = [idxTagged, idxNonTagged]; % & MechResponse;
matDims = size(idxMat);
figure('Name', 'SNR Violin Plots', 'Color', 'White');
nGroups = matDims(2);

yMax = [max(clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Control_4mW_Rate_Spont_0_to_1s),...
    max(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_0_to_1s)];
yMax = max(yMax);
yMin = [min(clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Control_4mW_Rate_Spont_0_to_1s),...
    min(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_0_to_1s)];
yMin = min(yMin);

for unitGroup = 1:nGroups
    subplot(1, nGroups, unitGroup)
    
    unitIds = idxMat(:,unitGroup);
    nUnits = sum(unitIds);
    rts = NaN(nUnits, 2);
    rts(:,1) = clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s(unitIds) - clInfo.Mech_Control_4mW_Rate_Spont_0_to_1s(unitIds);
    rts(:,2) = clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s(unitIds) - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_0_to_1s(unitIds);
    
    subplot(1,2,unitGroup)
    hold on
    
    
    [yCtr, xCtr] = ksdensity(rts(:,1),'Bandwidth',0.7);
    patch(yCtr * -1,xCtr,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
    [yCtr, xCtr] = ksdensity(rts(:,2),'Bandwidth',0.7);
    patch(yCtr,xCtr,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
    
    
    plot(zeros(nUnits,1) - 0.01, rts(:,1),...
        'LineStyle','none', 'Marker', '*', 'Color', [1,0,0], 'MarkerSize', 5);
    
    plot(zeros(nUnits,1) + 0.01, rts(:,2),...
        'LineStyle','none', 'Marker', '*', 'Color', [0,0,1], 'MarkerSize', 5);
    
    
    
    
    legend('Mech Control','Mech Laser', 'Location','northwest');
    
    ax = gca;
    ax.Title.String = ids{unitGroup};
    ax.FontName = 'Arial';
    ax.FontSize = 25;
    ax.YAxis.Limits = [yMin, yMax];
    
end
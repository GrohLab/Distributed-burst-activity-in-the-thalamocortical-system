function fig =...
    plotClusterReactivity(PSTH,trig,sweeps,timeLapse,binSz,IDe,expName,stims,IDs)
%PLOTCLUSTERREACTIVITY gets the individual PSTH from the getPSTH function
%and plots it for each cluster. It also plots the population PSTH on the
%bottom of the figure and draws a line for the beguinning and end of the
%stimuli.
%   fig = PlotClusterReactivity(PSTH, trig, sweeps, timeLapse, binSz,...
%               IDe, expName, stims, IDs)
%
%       INPUTS:
%           PSTH - a MxN matrix containing the counts for the M individual
%           clusters along the binned N time scale.
%           trig - a 1xT array containing the stimulus probability without
%           binning.
%           sweeps - number of trials for the given experiment.
%           timeLapse - peri-stimulus window. 1x2 array having the time
%           before and after the stimulus, both in seconds.
%           binSz - Bin size in seconds
%           IDe - (1+M)x1 cell array contining the names for the stimulus
%           and for the clusters involved.
%           expName - String variable containing the text for the figure
%           name.
%
%       OUTPUT:
%           fig - Figure object containing all the plots created.
%
%  Emilio Isaias-Camacho @ GrohLab 2019

% If the trigger is the whisker stimulating device, then the color of the
% trigger will be green. Otherwise, it would be light blue for the laser.
clr = defineColorForStimuli(IDe);

fthAxFlag = false;
if ~exist('stims','var') || isempty(stims)
    totlX = 4;
elseif ~isempty(stims)
    fthAxFlag = true;
    totlX = 5;
end

% Number of clusters and number of bins in the PSTH
[Ncl, Npt] = size(PSTH);
% Normalizing the PSTH to the maximal value on each cluster.
PSTHn = PSTH./max(PSTH,[],2);
fig = figure('Name',expName,'Color',[1,1,1]);
ax1 = subplot(totlX,1,1:3,'Parent',fig);
psthTX = linspace(timeLapse(1),timeLapse(2),Npt);
clrmp = parula;
% clrmp = defineWhYellRedColormap;
% clrmp = defineBlueRedColormap();
imagesc(ax1,'XData',psthTX,'CData',PSTHn);
colormap(ax1,clrmp)

trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));
trigTXc = linspace(timeLapse(1),timeLapse(2),size(stims,2));

% Plotting lines for depicting the on- and offset of the trigger
tObj = StepWaveform(trig,1/mean(diff(trigTX)));
tSubs = tObj.subTriggers;
if ~isempty(tSubs)
    for ct = 1:size(tSubs,1)
        for nf = 1:size(tSubs,2)
            line(ax1,'XData',repmat(trigTX(tSubs(ct,nf)),1,2),'YData',[0.5,Ncl+0.5],...
                'LineStyle',':','Color',clr,'LineWidth',2)
        end
    end
else
    trig = logical(trig/max(trig(:)));
    line(ax1,'XData',trigTX,'YData',(Ncl+1.5)*trig - 0.5,...
        'LineStyle',':','Color',clr,'LineWidth',2)
end

% Formatting the heatmap
ax1.YLim = [0.5,size(PSTH,1)+0.5];
ax1.XLim = [timeLapse(1)-binSz/2, timeLapse(2)+binSz/2];
ax1.YTick = 1:Ncl;
ax1.YTickLabel = IDe(2:end);
ax1.YAxis.Label.String = sprintf('Cluster ID_{%d #clusters}', size(PSTH,1));
ax1.XAxis.Visible = 'off';

% Plotting the population PSTH together with the trigger probability
ax2 = subplot(totlX,1,4,'Parent',fig);
popPSTH = sum(PSTH,1,'omitnan')/(Ncl * sweeps);
plot(ax2,psthTX,popPSTH,'Color',[0.8,0.8,0.8],'DisplayName','Population PSTH')
axLabel = 'Population activity';
yyaxis(ax2,'right')
plot(ax2,trigTX,trig,'LineWidth',1.5,'Color',clr,'DisplayName',IDe{1},...
    'LineStyle',':')

% Formatting the population PSTH plot
ax2.YAxis(1).Label.String = axLabel;
% ax2.YAxis(1).Limits = [0,1.01];
ax2.YAxis(2).Limits = [0,1.01];
ax2.YAxis(2).Color = 'k';
ax2.YAxis(2).Label.String = 'Stimulus probability';


ax2.XLabel.String = sprintf('Time_{%.2f ms} [s]',binSz*1e3);
ax2.XLim = [timeLapse(1), timeLapse(2)];
ax2.Box = 'off';
ax2.ClippingStyle = 'rectangle';

lgnd = legend(ax2,'show'); set(lgnd, 'Location','best', 'Box', 'off');
linkaxes([ax1,ax2],'x')
title(ax1,[expName,sprintf(' %d trials',sweeps)])

if fthAxFlag
    ax3 = subplot(totlX,1,5,'Parent',fig);
    [r,c] = size(stims);
    if r < c
        stims = stims';
    end
    
    for cs = 1:min(r,c)
        if exist('IDs','var')
            plot(ax3,trigTXc,stims(:,cs),'LineStyle','-.','LineWidth',0.5,...
                'DisplayName', IDs{cs})
        else
            plot(ax3,trigTXc,stims(:,cs),'LineStyle','-.','LineWidth',0.5)
        end
        
        ax3.Children(1).Color = defineColorForStimuli(IDs(cs));
        
        if cs == 1
            ax3.NextPlot = 'add';
        end
    end
    lgnd = legend(ax3,'show'); set(lgnd, 'Location', 'best', 'Box', 'off');
    ax3.Box = 'off';
    ax3.XLim = [timeLapse(1), timeLapse(2)];
    ax3.XAxis.Visible = 'off';
    ax3.YAxis.Visible = 'off';
    linkaxes([ax1,ax2,ax3],'x')
end
end

function clr = defineColorForStimuli(IDe)
if contains(IDe{1},'piezo','IgnoreCase',true) ||...
        contains(IDe{1},'whisker','IgnoreCase',true)
    clr = [0, 0.8, 0];
elseif contains(IDe{1},'light','IgnoreCase',true) ||...
        contains(IDe{1},'laser','IgnoreCase',true)
    clr = [80, 187, 211]/255;
else
    clr = [165, 70, 87]/255; 
end
end

function clrmp = defineBlueRedColormap()
clrmp = [linspace(0,1,256)', zeros(256,1), linspace(1,0,256)'];
end

function clrmp = defineWhYellRedColormap()
%Helper function to store the color map
clrmp = [1.0000    1.0000    1.0000
    1.0000    0.9909    0.9727
    1.0000    0.9818    0.9455
    1.0000    0.9727    0.9182
    1.0000    0.9636    0.8909
    1.0000    0.9545    0.8636
    1.0000    0.9455    0.8364
    1.0000    0.9364    0.8091
    1.0000    0.9273    0.7818
    1.0000    0.9182    0.7545
    1.0000    0.9091    0.7273
    1.0000    0.9000    0.7000
    1.0000    0.8909    0.6727
    1.0000    0.8818    0.6455
    1.0000    0.8727    0.6182
    1.0000    0.8636    0.5909
    1.0000    0.8545    0.5636
    1.0000    0.8455    0.5364
    1.0000    0.8364    0.5091
    1.0000    0.8273    0.4818
    1.0000    0.8182    0.4545
    1.0000    0.8091    0.4273
    1.0000    0.8000    0.4000
    0.9905    0.7746    0.3810
    0.9810    0.7492    0.3619
    0.9714    0.7238    0.3429
    0.9619    0.6984    0.3238
    0.9524    0.6730    0.3048
    0.9429    0.6476    0.2857
    0.9333    0.6222    0.2667
    0.9238    0.5968    0.2476
    0.9143    0.5714    0.2286
    0.9048    0.5460    0.2095
    0.8952    0.5206    0.1905
    0.8857    0.4952    0.1714
    0.8762    0.4698    0.1524
    0.8667    0.4444    0.1333
    0.8571    0.4190    0.1143
    0.8476    0.3937    0.0952
    0.8381    0.3683    0.0762
    0.8286    0.3429    0.0571
    0.8190    0.3175    0.0381
    0.8095    0.2921    0.0190
    0.8000    0.2667         0
    0.7800    0.2633         0
    0.7600    0.2600         0
    0.7400    0.2567         0
    0.7200    0.2533         0
    0.7000    0.2500         0
    0.6800    0.2467         0
    0.6600    0.2433         0
    0.6400    0.2400         0
    0.6200    0.2367         0
    0.6000    0.2333         0
    0.5800    0.2300         0
    0.5600    0.2267         0
    0.5400    0.2233         0
    0.5200    0.2200         0
    0.5000    0.2167         0
    0.4800    0.2133         0
    0.4600    0.2100         0
    0.4400    0.2067         0
    0.4200    0.2033         0
    0.4000    0.2000         0];
end
function [Figs, Results] =...
    scatterSignificance(Results, Counts, CondNames, Delta_t, clID)
%SCATTERSIGNIFICANCE takes the statistical results from the statTests and
%plots them against each other (conditions)


meanfr = cellfun(@(x) mean(x,2)/Delta_t, Counts, 'UniformOutput', false);
mxfr = cellfun(@max, meanfr);
mxfr = round(mxfr*1.15, -1);
% mxfr = max(mxfr);
if ~min(mxfr(:))
    mxfr(mxfr==0) = 10;
end
Nr = numel(Results);
Figs = gobjects(Nr,1);
ax = gobjects(2,1);
n = getHesseLineForm([1,0]);
axLabels = {'Spontaneous_{', 'Evoked_{'};
Ncond = numel(CondNames);
Ne = size(Results(1).Activity(1).Pvalues,1);
Hc = false(Ne, Nr*2 - Ncond);
hCount = 1;
for cr = 1:Nr
    combCell = textscan(Results(cr).Combination,'%d %d\t%s');    
    cond1 = double(combCell{1}); cond2 = double(combCell{2});
    Figs(cr) = figure('Color',[1,1,1],'Visible','off','Units','normalized');
    actvty = Results(cr).Activity(1).Type;
    if contains(actvty,'condition')
        figType = 1;
    else
        figType = 2;
    end
    csp = 1;
    while csp <= figType
        actvty = Results(cr).Activity(csp).Type;
        H = Results(cr).Activity(csp).Pvalues < 0.05;
        Hc(:,hCount) = H;
        hCount = hCount + 1;        
        ttle = sprintf('%s: %s vs. %s',actvty, CondNames{cond1},...
            CondNames{cond2});
        ax(csp) = subplot(1,figType,csp,'Parent',Figs(cr));
        switch actvty
            case 'Spontaneous'
                aslSubX = 1; aslSubY = 1;
            case 'Evoked'
                aslSubX = 2; aslSubY = 2;
            otherwise
                aslSubX = 1; aslSubY = 2;
        end
        xaxis = meanfr{cond1, aslSubX}; yaxis = meanfr{cond2, aslSubY};
        % Directionality of change
        dPts = [xaxis, yaxis] * n;
        Results(cr).Activity(csp).Direction = dPts;
        scatter(ax(csp), xaxis, yaxis, 'DisplayName', 'Clusters');
        grid(ax(csp), 'on');
        text(ax(csp), xaxis, yaxis, clID)
        grid(ax(csp), 'minor'); axis('square'); 
        axMx = max(mxfr(cond1,aslSubX),mxfr(cond2,aslSubY));
        % axis(ax(csp), [0, mxfr(cond1,aslSubX), 0, mxfr(cond1,aslSubY)]);
        axis(ax(csp), [0, axMx, 0, axMx]);
        % axis(ax(csp), [0, mxfr(cond1), 0, mxfr(cond1)]);
        ax(csp).NextPlot = 'add';
        % axMx = max(mxfr(cond1),mxfr(cond2));
        line(ax(csp), 'XData', [0, axMx],...
            'YData', [0, axMx],...
            'Color', [0.8, 0.8, 0.8], 'LineStyle', '--',...
            'DisplayName', 'y = x');
        title(ax(csp),ttle); 
        xlabel(ax(csp), [axLabels{aslSubX},CondNames{cond1},'} [Hz]']); 
        ylabel(ax(csp), [axLabels{aslSubY},CondNames{cond2},'} [Hz]']);
        scatter(ax(csp),xaxis(H), yaxis(H),'SizeData',30,...
            'DisplayName', 'Significant', 'Marker', 'pentagram')
        [mdl,~,rsq] = fit_poly(xaxis, yaxis, 1);
        line(ax(csp),'XData',[0, axMx],...
            'YData', [0*mdl(1) + mdl(2), axMx*mdl(1) + mdl(2)],...
            'DisplayName', sprintf('Population trend %.3f (%.1fx %+.1f)',...
            rsq, mdl(1), mdl(2)), 'LineStyle', ':', 'LineWidth', 0.5)
        csp = csp + 1;
        if aslSubX == 1 && aslSubY == 2
            H2 = Results(cr).Activity(csp).Pvalues < 0.05;
            scatter(ax(csp-1), xaxis(H2), yaxis(H2), '+', 'DisplayName', 'Shuffled')
        end
        if csp == 3
            Figs(cr).OuterPosition =...
                [Figs(cr).OuterPosition(1:2), 0.5344, 0.4275];
        end
        lgnd = legend(ax(csp-1),'show'); lgnd.Box = 'off'; 
        lgnd.Location = 'best';
    end
end
set(Figs,'Visible','on')
end


function [timePts, xslope] = getWaveformCriticalPoints(avWaves, fs)
% GETWAVEFORMCRITICALPOINTS returns the time points for which the first and
% second derivatives of the signal equals zero for the minima and/or maxima
% and the inflection points for a matrix NxM where N is the numbe of
% samples and M the number of signals.
%        [timePts, xslope] = getWaveformCriticalPoints(avWaves, fs)

% Zero-crossing variables
% [dw, zcIdx, zcSubs, zcN, zcSel] = computeDerivativePoints(avWaves, 1);
[dw, ~, zcSubs, ~, zcSel] = computeDerivativePoints(avWaves, 1);
% Fiducial point variables
% [ddw, fpIdx, fpSubs, fpN, fpSel] = computeDerivativePoints(avWaves, 2);
[ddw, ~, fpSubs, ~, fpSel] = computeDerivativePoints(avWaves, 2);
dt = 1/fs;
[Nt, Ncl] = size(avWaves);
tx = (0:Nt-1)/fs;
timePts = cell(Ncl,2);
xslope = timePts;
for ccl = 1:Ncl
    %{
%     lnIdx = arrayfun(@(x) x:x+1, zcSubs(zcSel(ccl)+1:zcSel(ccl+1)),...
%         funOpts{:});
%     dMdl = cellfun(@(x) fit_poly(tx(x), dw(x,ccl), 1), lnIdx, funOpts{:});
%     xz = cellfun(@(x) (-x(2)/x(1)) + dt/2, dMdl);
    %}
    [xz, xzs] = computeZeroCrossings(dt, ccl, zcSel, zcSubs, dw, tx, 1);
    [fp, fps] = computeZeroCrossings(dt, ccl, fpSel, fpSubs, ddw, tx, 2);
%{
%     fp = zeros(fpN(ccl),1);
%     for cfp = 1:numel(fpSubs)
%         lnIdx = fpSubs(cfp):fpSubs(cfp)+1;
%         mdl = fit_poly(tx(lnIdx), ddw(lnIdx, ccl), 1);
%         fp(cfp) = -mdl(2)/mdl(1) + dt;
%     end
%     lmt = abs(median(devs(ccl,:)));
%     xzImp = abs(interpolateTimeSeries(tx, devs(ccl,:), xz)) > lmt;
%     fpImp = abs(interpolateTimeSeries(tx, devs(ccl,:), fp)) > lmt;
%}
    timePts(ccl,:) = {xz, fp};
    xslope(ccl,:) = {xzs, fps};
end


end

function [dw, dIdx, dSubs, dN, dSel] = computeDerivativePoints(...
    signalMatrix, dOrder)
dw = diff(signalMatrix, dOrder, 1); dIdx = diff(sign(dw),1,1) ~= 0;
[dSubs,~] = find(dIdx); dN = sum(dIdx,1); dSel = [0,cumsum(dN)];
end

function [zCrossings, zxSlopeDir] =...
    computeZeroCrossings(dt, ci, dSel, dSubs, dw, tx, dOrder)
funOpts = {'UniformOutput', 0};
lnIdx = arrayfun(@(x) x:x+1, dSubs(dSel(ci)+1:dSel(ci+1)), funOpts{:});
dMdl = cellfun(@(x) fit_poly(tx(x), dw(x,ci), 1), lnIdx, funOpts{:});
dBase = mod(dOrder,2) + 1;
zCrossings = cellfun(@(x) (-x(2)/x(1)) + dt/dBase, dMdl);
zxSlopeDir = cellfun(@(x) sign(x(1)), dMdl);
end
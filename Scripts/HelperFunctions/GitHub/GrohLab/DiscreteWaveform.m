classdef (Abstract) DiscreteWaveform < GeneralWaveform
    %DISCRETEWAVEFORM This class would be abstract class to have
    %   Detailed explanation goes here
    
    properties (Abstract, SetAccess = 'private') 
        Triggers 
    end
    
    properties (SetAccess = 'private')
        FirstOfTrain
        LastOfTrain
        Count
        Delta
    end
    
    properties (Access = 'public')
        MinIEI (1,1) = NaN; % Time in seconds
    end
    
    methods
        function obj = DiscreteWaveform(data, samplingFreq, units, title)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@GeneralWaveform(data, samplingFreq, units, title);
            
        end 
        function h = plot(obj,varargin)
            h = plot@GeneralWaveform(obj,varargin{:});
            ax = get(h,'Parent');
            set(ax,'NextPlot','add');
            TG = obj.Triggers;
            tUp = find(TG(:,1));
            tDn = find(TG(:,2));
            spkOpls = diff([tUp,tDn],1,2);
            ttx = tUp / obj.SamplingFreq;
            plot(ax,ttx,obj.Data(obj.Triggers(:,1)),...
                'LineStyle','none','Marker','.');
            if sum(abs(spkOpls)) > numel(spkOpls)
                ttd = tDn / obj.SamplingFreq;
                plot(ax,ttd,obj.Data(obj.Triggers(:,2)),...
                    'LineStyle','none', 'Marker','.');
            end            
        end
        
        % Get the first event of a train
        function fot = get.FirstOfTrain(obj)
            evntTms = obj.loadTimeData();
            if ~isnan(obj.MinIEI)
                fot = DiscreteWaveform.firstOfTrain(evntTms,obj.MinIEI);
            else
                fprintf('The minimum inter-event interval MinIEI is not set!\n')
                fprintf('The function will use the mean MinIEI...\n')
                [fot, obj.MinIEI] = DiscreteWaveform.firstOfTrain(evntTms);
            end
        end
        
        % Get the last event of a train
        function lot = get.LastOfTrain(obj)
            evntTms = flip(obj.loadTimeData());
            if ~isempty(obj.MinIEI)
                lot = flip(DiscreteWaveform.firstOfTrain(evntTms,obj.MinIEI));
            else
                fprintf('The minimum inter-event interval MinIEI is not set!\n')
                fprintf('The function will use the mean MinIEI...\n')
                [lot, obj.MinIEI] = flip(DiscreteWaveform.firstOfTrain(evntTms));
            end
        end
        
        % Get the number of events in a train
        function c = get.Count(obj)
            fot = obj.FirstOfTrain;
            lot = obj.LastOfTrain;
            feSubs = find(fot);
            leSubs = find(lot);
            c = (leSubs - feSubs) + 1;
        end
        
        % Get the duration of the train in seconds
        function d = get.Delta(obj)
            fot = obj.FirstOfTrain;
            lot = obj.LastOfTrain;
            feSubs = find(fot);
            leSubs = find(lot);
            d = abs(feSubs-leSubs)./obj.SamplingFreq;
        end
        
        function tdTms = loadTimeData(obj)
            fs = obj.SamplingFreq;
            try
                % In case of stepWaveform
                trigs = obj.subTriggers;
            catch
                % In case of spikeWaveform
                trigs = obj.Data(:);
            end
            tdTms = trigs(:,1)./fs;
        end
    end
    
    %% Static public methods
    methods (Static, Access = 'public')
        % Recreate the signal with logic values.
        function logicalTrace = subs2idx(subs,N)
            if size(subs,2) == 2
                logicalTrace = false(1,N);
                for cmt = 1:size(subs,1)
                    logicalTrace(subs(cmt,1):subs(cmt,2)) = true;
                end
            else
                [Nr, Nc] = size(subs);
                logicalTrace = false(Nr * (Nr < Nc) + Nc * (Nc < Nr),...
                    N);
                subsClass = class(subs);
                switch subsClass
                    case 'cell'
                        for cmt = 1:size(subs,1)
                            logicalTrace(cmt,subs{cmt}) = true;
                        end
                    case {'single','double'}
                        logicalTrace(1,subs) = true;
                    otherwise
                        fprintf('Case not yet implemented...\n')
                end
            end
        end
        
        % Get a semi-logic trigger signal
        function semiLogicSignal = SCBrownMotion(RaF)
            [R,C] = size(RaF);
            if ~any([R == 2,C == 2])
                disp('What kind of rising and falling edges were given?')
                semiLogicSignal = RaF;
                return
            end
            if R > C
                RaF = RaF';
            end
            semiLogicSignal = cumsum(RaF(1,:) - RaF(2,:));
        end
        
        % Get the first true value of a logic pulse
        function [frstSpks, minIpi, Pks, Sps] = firstOfTrain(spkTimes, minIpi)
            % OUTPUTS a logical index for the edges which are the first in
            % time for the step pulse train.
            Ipi = abs(diff(spkTimes));
            if ~exist('minIpi','var')
                % minIpi = mean(Ipi);
                minIpi = DiscreteWaveform.computeIpiThresh(Ipi);
                if isempty(minIpi)
                    minIpi = 1;
                end
            end
            Pks = Ipi < minIpi;
            Sps = DiscreteWaveform.addFst(~Pks,true);
            frstSpks = DiscreteWaveform.addLst(Sps(1:end-1) & Pks,false);
        end 
        
        function [frstSpks, minIpi, Pks, Sps] = lastOfTrain(spkTimes, minIpi)
            % OUTPUTS a logical index for the edges which are the first in
            % time for the step pulse train.
            Ipi = abs(diff(spkTimes));
            if ~exist('minIpi','var')
                % minIpi = mean(Ipi);
                minIpi = DiscreteWaveform.computeIpiThresh(Ipi);
                if isempty(minIpi)
                    minIpi = 1;
                end
            end
            Pks = flip(Ipi) < minIpi;
            Sps = DiscreteWaveform.addFst(~Pks,true);
            frstSpks = flip(DiscreteWaveform.addLst(Sps(1:end-1) & Pks,false));
        end 
        
        function ipiThresh = computeIpiThresh(Ipi)
            funOpts = {'UniformOutput', false};
            [binCenters, binEdges, lData, ts] = prepareLogBinEdges(Ipi, 128);
            pCount = histcounts(lData, binEdges);
            mu = mean(pCount); sig = std(pCount);
            populatedBins = pCount ~= 0;
            if sum(populatedBins)/numel(pCount) < 0.5
                % Discrete intervals
                fprintf(1, 'Discrete ')
                [~, bigGap] = max(diff(binEdges(populatedBins)));
                populatedSubs = find(populatedBins);
                ipiMdl = round(populatedSubs([bigGap, bigGap+1])*[0.33;0.67]);
                % ipiThresh = 10.^max(binEdges([0,(pCount.*binCenters)] < 0));
                ipiThresh = 10.^mean([binCenters(ipiMdl),...
                    max(binEdges([(pCount.*binCenters), 0] < 0))]);
            else
                % Continuous (spikes)
                fprintf(1, 'Continuous ')
                % Fitting a spline to smooth the sampling artifacts.
                spl = fitSpline(binCenters, pCount, 2,...
                    range(binCenters)/9, 1/3);
                % Getting the critical points of the smoothed PDF
                [crPts, crPtSl] = getWaveformCriticalPoints(spl(:), 1/ts);
                crPts = cellfun(@(x) x + binCenters(1), crPts, funOpts{:});
                % Interpolating the 'y' values of the critical points
                yCrPt = cellfun(@(x) interp1(binCenters, spl, x), crPts,...
                    funOpts{:});
                % Sorting the peaks amplitude to get the tallest crossing
                % at least 1 standard deviation
                [lmm, peakOrd] = sort(crPtSl{1}.*yCrPt{1}, 'ascend');
                pdfPeak = abs(lmm(lmm < 0) - mu)/sig > 1;
                Ng = sum(pdfPeak);
                % if there are more than 1 peak, chances are that the
                % distribution is bimodal.
                if Ng > 1
                    % Search for the valley between the peaks
                    
                else
                    % Search for an ISI value or return 1 millisecond
                    ipiThresh = 1e-3;
                end
                
            end
            fprintf(1, 'distribution of intervals\n')
%             ipiThresh = mean(Ipi);
%             uIpi = uniquetol(round(Ipi,2));
%             uIpi_p = diff(uIpi);
%             zUip = abs(zscore(uIpi_p));
%             uIpi_d = uIpi_p/max(uIpi_p(zUip < 2));
%             if numel(uIpi) >= numel(Ipi)*0.1
%                 % Continuous intervals
%                 % To be implemented
%                 fprintf(1, 'Continuous ')
%             else
%                 % Discrete intervals
%                 fprintf(1, 'Discrete ')
%                 Ipi_p = abs(diff(Ipi));
%                 tatGap = find(diff(uIpi) > 1,1);
%                 mxIpi = max(Ipi(Ipi_p < 2/3e4));
%                 if isempty(mxIpi)
%                     mxIpi = 0;
%                 end
%                 ipiThresh = [mxIpi, uIpi(tatGap)] * ...
%                     [1-uIpi_d(tatGap); uIpi_d(tatGap)];
%                 if any(ismember(uIpi, ipiThresh))
%                     ipiThresh = ipiThresh * 1.2;
%                 end
%             end
%             fprintf(1, 'distribution of intervals\n')
        end
    end
    
    %% Static private methods
    methods (Static, Access = 'private')
        function new_array = addFst(array,element)
            new_array = cat(find(size(array)~=1), element, array);
        end
        
        function new_array = addLst(array,element)
            new_array = cat(find(size(array)~=1), array, element);
        end
    end
end

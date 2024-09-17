function [discreteStack, continuouStack] =...
    getStacks(spT,alignP,ONOFF,timeSpan,fs,fsLFP,consEvents,varargin)
% GETSTACK returns a stack of spikes aligned to a certain event ''alignT''
% considering the events in the cell array ''consEvents''. The alignment
% can be done with the on-set or the off-set of the triggers using the
% ONOFF option. The default is 'on'. The timeSpan is given in a 1x2 array
% with the time before the trigger as the first element and the time after
% the trigger as the second element. The sampling frequency is always
% important and it is given in Hz. The LFP and whisker movement are
% examples of triggered averages from continuous signals. Finally, the
% consEvents are 
discreteStack = NaN;
continuouStack = NaN;
%% Computing the size of the stack
if isa(alignP,'logical')
    alWf = StepWaveform(alignP,fs,'on-off','Align triggers');
    alignP = alWf.Triggers;
end
switch ONOFF
    case 'on'
        disp('Considering onset of the triggers')
    case 'off'
        disp('Considering offset of the triggers')
    otherwise
        disp('Unrecognized trigger selection. Considering onsets')
        ONOFF = 'on';
end

[auxR,auxC] = size(alignP);
% if auxR < auxC
%     alignP = alignP';
% end
Na = auxR;
if auxC > 2 || auxC < 1
    fprintf(['Warning! The alignment matrix is expected to have ',...
        'either only rising or rising and falling edges time indices.\n'])
    return;
end

%% Considered events arrangement.
Ne = 0;
if exist('consEvents','var') && ~isempty(consEvents)
    typ = whos('consEvents');
    switch typ.class
        case {'double', 'single'}
            [rws, cols] = size(consEvents);
            if rws < cols
                Ne = rws;
                consEvents = consEvents';
            else
                Ne = cols;
            end
        case 'cell'
            Ne = length(consEvents);
            consEvents2 = consEvents;
            evntTrain = cellfun(@islogical,consEvents);
            % Converting the logical event trains into indices
            for ce = 1:Ne
                if ~isempty(consEvents2{ce})
                    if evntTrain(ce)
                        % Logical event
                        stWv = StepWaveform(consEvents{ce},fs);
                        consEvents2{ce} = stWv.Triggers;
                    elseif ~sum(round(consEvents{ce}) - consEvents{ce})
                        % Subscript event
                        consEvents2(ce) = consEvents(ce);
                    else
                        % 'Raw' signal
                        fprintf(1,'Transforming the considered event %d to',ce)
                        fprintf(1,' logical\n');
                        aux = abs(consEvents2{ce});
                        consEvents2{ce} = aux > mean(aux);
                    end
                else
                    % Empty variable will be deleted
                    fprintf(1,'Input event %d is empty and will be deleted\n',...
                        ce)
                    consEvents2(ce) = [];
                end
            end
            consEvents = consEvents2;
        otherwise
            fprintf('The events to consider are not in recognized format.\n')
    end
end
%% Preallocation of the discrete stack:
prevSamples = ceil(abs(timeSpan(1)) * fs);
postSamples = ceil(timeSpan(2) * fs);
Nt = prevSamples + postSamples + 1;
try
    discreteStack = false(2+Ne,Nt,Na);
catch
    fprintf(1,'The requested array size exceeds the free RAM memory.\n')
    fprintf(1,'Please, make note of where this message appeared.\n')
    discreteStack = [];
    return;
end
% Creation of the logical spike train
if isnumeric(spT) && ~sum(round(spT) - spT) && nnz(spT) == max(size(spT))
    mxS = spT(end) + Nt;
    spTemp = false(1,mxS);
    spTemp(spT) = true;
    spT = spTemp;
elseif ~sum(round(spT) - spT) && nnz(spT) == max(size(spT))
    spT = round(spT * fs);
    mxS = spT(end) + Nt;
    spT = StepWaveform.subs2idx(spT,mxS);
end
%% Preallocation of the continuous stack:
if ~exist('fsLFP','var')
    fsLFP = fs;
end 
fsConv = fsLFP/fs;
% Signal validation
Ns = numel(varargin);
sparseFlag = false;
if Ns
    if Ns == 1
        [Nrow, Ncol] = size(varargin{1});
        if (Nrow ~= 1 && Ncol > Nrow) || (Ncol ~= 1 && Nrow > Ncol)
            Ns = Nrow * (Nrow < Ncol) + Ncol * (Ncol < Nrow);
        end
    end
    signalCheck = cellfun(@isnumeric, varargin);
    signalCheck2 = cellfun(@length, varargin);
    signalCheck3 = cellfun(@iscell, varargin);
    signalCheck4 = cellfun(@issparse, varargin);
    if any(signalCheck3)
        varargin = varargin{1};
        Ns = numel(varargin);
        lngths = cellfun(@length,varargin);
        lnthCheck = std(lngths);
        if any(lnthCheck./lngths > 1)
            outlier = zscore(lngths).^2 > 0.5;
            fprintf(1,'Warning! One of the considered continuous signals')
            fprintf(1,' length deviates considerably from the rest.\n')
            fprintf(1,'This signal(s) is (are) going to be deleted!\n')
            varargin(outlier) = [];
            Ns = Ns - sum(outlier);
        end
        MAX_CONT_SAMP = min(cellfun(@length,varargin));
        % Convert all signals into a row vector
        RoC = cellfun(@isrow,varargin);
        if ~all(RoC)
            varargin(~RoC) = cellfun(@transpose,varargin(~RoC),...
                'UniformOutput',false);
        end
        % Convert all signals to double
        clss = cellfun(@class, varargin, 'UniformOutput', false);
        if numel(unique(clss)) > 1
            diffTypeFlag = ...
                ~cellfun(@contains, clss, repmat({'double'},numel(clss),1));
            varargin(diffTypeFlag) =...
                cellfun(@double, varargin(diffTypeFlag), 'UniformOutput', false);
        end
    end
    if sum(signalCheck) ~= Ns && ~any(signalCheck3)
        fprintf('Discarding those inputs which are not numeric...\n')
        disp(varargin(~signalCheck))
        varargin(~signalCheck) = [];
        Ns = Ns - sum(~signalCheck);
        signalCheck2(~signalCheck) = [];
    end
    if any(signalCheck)
        if std(signalCheck2) ~= 0
            fprintf('The signals have not the same length...\n')
            fprintf('Considering the smallest: %d\n',min(signalCheck2))
            MAX_CONT_SAMP = min(signalCheck2);
        else
            MAX_CONT_SAMP = signalCheck2(1);
        end
    end
    % prevSamplesLFP = ceil(timeSpan(1) * fsLFP);
    % postSamplesLFP = ceil(timeSpan(2) * fsLFP);
    % NtLFP = prevSamplesLFP + postSamplesLFP + 1;
    % continuouStack = zeros(Ns,NtLFP,Na,'single');
    if any(signalCheck4)
        sparseFlag = true;
    end
end
prevSamplesLFP = ceil(abs(timeSpan(1)) * fsLFP);
postSamplesLFP = ceil(timeSpan(2) * fsLFP);
NtLFP = prevSamplesLFP + postSamplesLFP + 1;
continuouStack = zeros(Ns,NtLFP,Na,'single');
if sparseFlag
    fprintf(1,'Ideal case: sparse output\n');
end
%% Cutting the events into the desired segments.
for cap = 1:Na
    % Considering the rising or the falling edge of the step function.
    if strcmp(ONOFF, 'on')
        segmIdxs = [alignP(cap,1)-prevSamples,alignP(cap,1)+postSamples];
        segmIdxsLFP = round([(alignP(cap,1)*fsConv)-prevSamplesLFP,...
            (alignP(cap,1)*fsConv)+postSamplesLFP]);
    elseif strcmp(ONOFF, 'off')
        segmIdxs = [alignP(cap,2)-prevSamples,alignP(cap,2)+postSamples];
        segmIdxsLFP = round([(alignP(cap,2)*fsConv)-prevSamplesLFP,...
            (alignP(cap,2)*fsConv)+postSamplesLFP]);
    end
    % The segments should be in the range of the spike train.
    % Validations for the discrete stack
    spSeg = false(1,Nt);
    if numel(spT) ~= 1
        if segmIdxs(1) >= 1 && segmIdxs(2) <= length(spT)
            spSeg = spT(segmIdxs(1):segmIdxs(2));
        elseif segmIdxs(1) <= 0
            fprintf('The viewing window is out of the signal range.\n')
            fprintf('%d samples required before the signal. Omitting...\n',...
                segmIdxs(1))
            segmIdxs(segmIdxs <= 0) = 1;
            try
                spSeg(Nt - segmIdxs(2) + 1:Nt) = spT(segmIdxs(1):segmIdxs(2));
            catch
                fprintf(1, 'Not enough samples given to create the trial %d\n',...
                    cap)
            end
        else
            fprintf('The viewing window is out of the signal range.\n')
            fprintf('%d samples required after the signal. Omitting...\n',...
                segmIdxs(2))
            segmIdxs(segmIdxs > length(spT)) = length(spT);
            spSeg(1:diff(segmIdxs)+1) = spT(segmIdxs(1):segmIdxs(2));
        end
    end
    discreteStack(2,:,cap) = spSeg;
    % Find 'overlapping' events in time of interest
    alignPeriod = getEventPeriod(alignP, {alignP}, ONOFF, cap,...
        prevSamples, postSamples);
    discreteStack(1,:,cap) = alignPeriod;
    % If there are continuous events given
    if Ne
        if isa(consEvents,'cell')
        discreteStack(3:2+Ne,:,cap) =...
            getEventPeriod(alignP, consEvents, ONOFF, cap,...
            prevSamples, postSamples);
        elseif isnumeric(consEvents)
            discreteStack(3:2+Ne,:,cap) =...
            getEventPeriod(alignP, {consEvents}, ONOFF, cap,...
            prevSamples, postSamples);
        else
        end
    end
    % Getting the continuous segments if the time window is in range of the
    % signal domain and validation for the continuous stack
    if Ns
        % segmIdxsLFP(1) >= 1 && segmIdxsLFP(2) <= MAX_CONT_SAMP
        
        % continuouStack(:,:,cap) = single(cell2mat(signalSegments'));
        continuouStack(:,:,cap) =...
            getContinuousSignalSegment(varargin, segmIdxsLFP, MAX_CONT_SAMP, NtLFP);
    end
end
fprintf(1,'Stacks constructed!\n')
end

% Aligning the events according to the considered time point. The inputs
% are time indices called Tdx as in Time inDeX.
function evntOn = getEventPeriod(alignTdx, evntTdx, ONOFF, cap, prev, post)
% Assuming that the considered events are always cells.
if isempty(evntTdx)
    evntOn = [];
    return;
else
    evntOn = false(numel(evntTdx),prev+post+1);
    for ce = 1:length(evntTdx)
        % Change of implementation. Considering the on set OR the offset of
        % a step pulse.
        if strcmp(ONOFF,'on')
            relTdx = evntTdx{ce} - alignTdx(cap,1);
        else
            relTdx = evntTdx{ce} - alignTdx(cap,2);
        end
        
        % Indexes for encountered rising _onIdx_, falling _offIdx_ edges,
        % or a bigger pulse envolving the viewing window _invIdx_.
        onIdx = relTdx(:,1) >= -prev & relTdx(:,1) < post; 
        if size(relTdx, 2) == 2
            offIdx = relTdx(:,2) >= -prev & relTdx(:,2) < post;
            invIdx = relTdx(:,1) < -prev & relTdx(:,2) > post;
        else
            offIdx = onIdx;
            invIdx = onIdx;
        end
        % Indexes considering the partial or complete step that falls into
        % the considered segment: _onIdx_ & _offIdx_
        % Indices considering the inclusion of the considered window into a
        % lengthy pulse of the considered event: _invIdx_
        allIdx = find(onIdx | offIdx | invIdx);
        initStep = relTdx(onIdx,1) + prev + 1;
        if size(relTdx, 2) == 2
            fnalStep = relTdx(offIdx,2) + prev + 1;
        else
            fnalStep = initStep;
        end
        % Indices considering the inclusion of the considered window into a
        % lengthy pulse of the considered event.
        %
        % HERE
        %
        
        % Event rising, falling or both into the alignment window as a
        % first condition. 
        stCount = 1;
        ndCount = 1;
        for cstp = 1:numel(allIdx)
            % If there is no rising edge found, then the step is considered
            % to start before the window and it will be true since the
            % start of the window. Otherwise, the obvious index is taken.
            if onIdx(allIdx(cstp))
                strt = initStep(stCount);
                stCount = stCount + 1;
            else
                strt = 1;
            end
            % Same case. If the falling edge is not found, the step is
            % considered to end in a later time point than the considered
            % window and it will be true from the rising edge to the end.
            if offIdx(allIdx(cstp))
                fnal = fnalStep(ndCount);
                ndCount = ndCount + 1;
            else
                fnal = post + prev + 1;
            end
            evntOn(ce,strt:fnal) = true;
        end
    end
end
end

function contSigSeg = getContinuousSignalSegment(signalCell, Idxs, N, Nt)
Ns = numel(signalCell);
if Ns == 1
    [Nrow, Ncol] = size(signalCell{1});
    if (Nrow ~= 1 && Ncol > Nrow) || (Ncol ~= 1 && Nrow > Ncol)
        Ns = Nrow * (Nrow < Ncol) + Ncol * (Ncol < Nrow);
    end
else
    
end
contSigSeg = zeros(Ns,Nt,'single');

noSigFlag = cellfun(@(x) isempty(x), signalCell);
if all(noSigFlag)
    return
else
    signalCell(noSigFlag) = [];
end

Subs = Idxs;
if Idxs(1) >= 1 && Idxs(2) <= N
    SegSubs = 1:Nt;
elseif Idxs(1) <= 0
    Subs(Idxs <= 0) = 1;
    SegSubs = Nt - Subs(2) + 1:Nt;
else
    Subs(Idxs > N) = N;
    SegSubs = 1:diff(Subs)+1;
end
%signalSegments = getSignalSegments(signalCell, Subs);
signalSegments = getSignalSegments();
signalMat = cell2mat(signalSegments);
if Ns == size(signalMat,1)
    if ~issparse(signalMat)
        contSigSeg(:,SegSubs) = signalMat;
    else
        contSigSeg(:,SegSubs) = full(signalMat);
    end
else
    contSigSeg(:,SegSubs) = cell2mat(signalSegments');
end

    %function sigSeg = getSignalSegments(signalCell,Subs)
    function sigSeg = getSignalSegments()
        
        % Transpose column vectors
        transpSign = cellfun(@iscolumn, signalCell);
        
        if any(transpSign)
            try
                signalCell(transpSign) =...
                    {signalCell{transpSign}'};
            catch
                auxSegm = cellfun(@transpose,signalCell(transpSign),...
                    'UniformOutput',false);
                signalCell(transpSign) = auxSegm;
            end
        end
        
        % Equalize data type (position)

        try
            sigSeg =...
                cellfun(...
                @(x) x(:,round((Subs(1):Subs(2)))),...
                signalCell, 'UniformOutput', false);
        catch ME
            fprintf('The is an issue with the continuous cell array.\n')
            fprintf('Worth debugging!\n')
            disp(ME.getReport)
            dbstop in getStacks at 378
        end
   
    end

        %{
% % Equalize data type
% clss = cellfun(@class, signalCell, 'UniformOutput', false);
% if numel(unique(clss)) > 1
%     diffTypeFlag = ...
%         ~cellfun(@contains, clss, repmat({'double'},numel(clss),1));
%     signalCell(diffTypeFlag) =...
%         cellfun(@double, signalCell(diffTypeFlag), 'UniformOutput', false);
% end
        %}

end


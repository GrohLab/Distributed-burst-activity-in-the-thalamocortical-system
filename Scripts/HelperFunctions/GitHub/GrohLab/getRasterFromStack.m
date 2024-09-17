function [relativeSpikeTimes, tx] =...
    getRasterFromStack(...
    discreteStack,... Discrete stack containing the alignment points
    kIdx,... Boolean array indicating which trials should be ignored
    koIdx,... Boolean array indicating which events should be kept.
    timeLapse,... 2 element array contaning the time before the trigger and the time after the trigger in seconds
    fs,... Original sampling frequency
    tmsORsubsFlag,... Subscripts or times boolean flag
    ERASE_kIDX... Boolean flag indicating if
    )
%GETRASTERFROMSTACK converts the logical positions of the spikes into
%spike times relative to the trigger time. It returns a cell array ExT,
%where E is the number of events and T is the number of triggers. The user
%has the 'flexibility' to choose which triggers to take out of the results
%with the koIdx input variable. This should be a vector of Tx1 elements.
%Similarly, the user can also select which events to include in the
%result by using the input variable koIdx, which should be Ex1 boolean
%vector.
%   [relativeSpikeTimes, tx] = getRasterFromStack(discreteStack, kIdx,
%                               koIdx, timeLapse, fs, ERASE_kIDX)
%       INPUTS:
%           discreteStack - boolean matrix sized ExSxT, where E is number
%           of events to consider, S is the number of samples, and T is the
%           number of triggers. This matrix is obtained by calling the
%           function getStacks. 
%           kIdx - Tx1 boolean vector which indicate which triggers should
%           be kicked out by having a true value.
%           koIdx - Ex1 boolean vector indicating which events should be
%           ignored for the results.
%           timeLapse - 2x1 numerical vector indicating the time in seconds
%           before and after the trigger
%           fs - sampling frequency
%           tmsORsubsFlag - boolean number to indicate the output to be
%           either subscripts or time stamps.
%           ERASE_kIDX - boolean scalar to indicate if the kicked out
%           events should be erased from the results. default false.
%       OUTPUTS:
%           relativeSpikeTimes - (E-(#koIdx))x(T-(#kIdx)) cell array
%           containing the relative time spikes in seconds for each event
%           and for each trigger.
%           tx - time axis for the time span given.
%   Emilio Isaias-Camacho @ GrohLab 2019
if ~exist('tmsORsubsFlag','var')
    tmsORsubsFlag = true;
end
[Ne, Nt, Na] = size(discreteStack);
% There will be for each neuron (or event) Na - !(kIdx) number of trials
relativeSpikeTimes = cell(Ne-(sum(~koIdx) + 1), Na);
iE = find(koIdx);
if isempty(iE)
    iE = 2;
else
    iE = reshape(iE,1,numel(iE));
    iE = [2, 2 + iE];
end
% Time axis
spIdx = 1;
tx = linspace(timeLapse(1),timeLapse(2),Nt);
for cse = iE
    % For each spike train
    for cap = Na:-1:1
        % For each alignment point
        if ~kIdx(cap)
            % If the alignment point should be ignored
            % upEdge = diff(discreteStack(cse,:,cap)) > 0;
            % if sum(upEdge) ~= 0
                % If there are events in this alignment point
                % downEdge = diff(discreteStack(cse,:,cap)) < 0;
                % isSpike = upEdge(1:end-1) - downEdge(2:end);
                % if sum(isSpike) == 0
                    % If the event contains spikes
                    spikeTimes = tx(squeeze(discreteStack(cse,:,cap)));
                    if tmsORsubsFlag
                        % Time
                        relativeSpikeTimes(spIdx,cap) = {spikeTimes};
                    else
                        % Subscripts
                        relativeSpikeTimes(spIdx,cap) =...
                            {round(fs*spikeTimes)};
                    end
                    %lvl = (cse - 2)*Na + cap;
                % end
            % end
        end
    end
    spIdx = spIdx + 1;
    
end
tx = seconds(tx);
if exist('ERASE_kIDX','var') && ERASE_kIDX
    % If the user chose to erase the kick out alignment trials
    relativeSpikeTimes(:,cellfun(@isempty,relativeSpikeTimes(1,:))) = [];
end




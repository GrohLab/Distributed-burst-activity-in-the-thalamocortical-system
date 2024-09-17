function [TypeOfSpike,BurstStarts,BurstEnds,BurstDurations,NumSpikesInBursts,numberBursts] = BurstDetect(spikeTimes_msec, BurstThreshold_msec, TailThreshold_msec, minNoSpikes)
TypeOfSpike = zeros(1,length(spikeTimes_msec));

if spikeTimes_msec(2) - spikeTimes_msec(1) < BurstThreshold_msec
    TypeOfSpike(1) = 1; % 1 for bursts, 2 for tails , 0 for individuals
end

for j = 2 : length(spikeTimes_msec)-1
    if spikeTimes_msec(j) - spikeTimes_msec(j-1) < BurstThreshold_msec || spikeTimes_msec(j+1) - spikeTimes_msec(j)< BurstThreshold_msec
        TypeOfSpike(j) = 1;
    end
end

if spikeTimes_msec(length(spikeTimes_msec)) - spikeTimes_msec(length(spikeTimes_msec)-1) < BurstThreshold_msec
    TypeOfSpike(length(spikeTimes_msec)) = 1;
end

%% omitting bursts having fewer spikes then desired
BreakPointAfter=[];
BreakPointBefore=[];
burstspikes = find(TypeOfSpike==1);
if ~isempty(burstspikes)
    for j = 1:length(burstspikes)-1
        if spikeTimes_msec(burstspikes(j+1))- spikeTimes_msec(burstspikes(j)) > BurstThreshold_msec
            BreakPointAfter = [BreakPointAfter burstspikes(j)];
        end
    end
    BreakPointAfter = [BreakPointAfter burstspikes(end)];
    
    for j = 2:length(burstspikes)
        if spikeTimes_msec(burstspikes(j))- spikeTimes_msec(burstspikes(j-1)) > BurstThreshold_msec
            BreakPointBefore = [BreakPointBefore burstspikes(j)];
        end
    end
    BreakPointBefore = [burstspikes(1) BreakPointBefore];
    
    ZeroIndexes = find(BreakPointAfter-BreakPointBefore+ 1<minNoSpikes); % spike indexes to be nullified
    for l=1:length(ZeroIndexes)
        TypeOfSpike(BreakPointBefore(ZeroIndexes(l)):BreakPointAfter(ZeroIndexes(l))) = 0;
    end
end

%% TAIL CALCULATION

if spikeTimes_msec(2) - spikeTimes_msec(1) <= TailThreshold_msec && ...
        TypeOfSpike(1) ~= 1
    TypeOfSpike(1) = 2; % 1 for bursts, 2 for tails , 0 for individuals
end

for j = 2 : length(spikeTimes_msec)-1
    if (spikeTimes_msec(j) - spikeTimes_msec(j-1) <= TailThreshold_msec)...
            && TypeOfSpike(j) ~= 1
        TypeOfSpike(j) = 2;
    end
end

if (spikeTimes_msec(length(spikeTimes_msec)) - spikeTimes_msec(length(spikeTimes_msec)-1) <= TailThreshold_msec) ...
        && TypeOfSpike(length(spikeTimes_msec)) ~= 1
    TypeOfSpike(length(spikeTimes_msec)) = 2;
end


%% omitting tails away from bursts and merging rest with the bursts

tailspikes = find(TypeOfSpike==2);
controlForLoop = TypeOfSpike;
resultOfmerging=zeros(1,length(TypeOfSpike));
%/loop till all the burst related tails merged and counted as bursts
while ~isequal(controlForLoop,resultOfmerging)
    controlForLoop = TypeOfSpike;
    for m = 1 : length(tailspikes)
        if (tailspikes(m)== 1)...
                && (spikeTimes_msec(tailspikes(m) + 1) - spikeTimes_msec(tailspikes(m)) <= TailThreshold_msec...
                && TypeOfSpike(tailspikes(m) + 1)==1)
            TypeOfSpike(tailspikes(m)) = 1;
        elseif tailspikes(m)== length(TypeOfSpike)...
                && (spikeTimes_msec(tailspikes(m))- spikeTimes_msec(tailspikes(m)-1) <= TailThreshold_msec...
                && TypeOfSpike(tailspikes(m) - 1)==1)
            TypeOfSpike(tailspikes(m)) = 1;
        elseif (tailspikes(m)~= 1) && tailspikes(m)~= length(TypeOfSpike)
            if(spikeTimes_msec(tailspikes(m) + 1) - spikeTimes_msec(tailspikes(m)) <= TailThreshold_msec...
                    && TypeOfSpike(tailspikes(m) + 1)==1) || (spikeTimes_msec(tailspikes(m))...
                    - spikeTimes_msec(tailspikes(m)-1) <= TailThreshold_msec...
                    && TypeOfSpike(tailspikes(m) - 1)==1)
                TypeOfSpike(tailspikes(m)) = 1;
            end
        end
    end
    resultOfmerging = TypeOfSpike;
    tailspikes = find(TypeOfSpike==2);
end

TypeOfSpike(TypeOfSpike ~= 1)=0;

%%Burst Starts  and Burst Ends
burstspikes = find(TypeOfSpike==1);
if ~isempty(burstspikes)
    BreakPoint=[];
    for j = 1:length(burstspikes)-1
        if spikeTimes_msec(burstspikes(j+1))- spikeTimes_msec(burstspikes(j)) > TailThreshold_msec
            BreakPoint = [BreakPoint burstspikes(j)];
        end
    end
    BreakPointEnds = [BreakPoint burstspikes(end)];
    BurstEnds = spikeTimes_msec(BreakPointEnds);
    
    BreakPoint=[];
    for j = 2:length(burstspikes)
        if spikeTimes_msec(burstspikes(j))- spikeTimes_msec(burstspikes(j-1)) > TailThreshold_msec
            BreakPoint = [BreakPoint burstspikes(j)];
        end
    end
    BreakPointStarts = [burstspikes(1) BreakPoint];
    BurstStarts = spikeTimes_msec(BreakPointStarts);
    NumSpikesInBursts = BreakPointEnds - BreakPointStarts +1;
    BurstDurations = BurstEnds - BurstStarts;
    numberBursts = length(BurstStarts);
else
    BurstStarts = NaN;
    BurstEnds= NaN;
    BurstDurations= NaN;
    NumSpikesInBursts= NaN;
    numberBursts= 0;
end

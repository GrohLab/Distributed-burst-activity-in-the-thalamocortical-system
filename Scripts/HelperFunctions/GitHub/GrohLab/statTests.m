function [Results, Counts] = statTests(dStack, condFlag, timeFlags, varargin)
%STATTESTS counts the spikes in the given time flags and compares them
%within and against other conditions indicated in condFlag. 
%   [Results, Counts] = statTests(dStack, condFlag, timeFlags, statTest)
%       INPUTS
%           - dStack - a NxTxA logical matrix which is created from the
%           getStacks function. N is number of events plus the considered
%           trigger, T is the time around the trigger, and A is the total
%           number of triggers.
%           - condFlag - a CxO logical matrix which indicates the condition
%           separation among the same trigger. C is the number of clusters
%           and O is the number of conditions available.
%           - timeFlags - a 2xT logical matrix which indicates the two time
%           windows to be compared.
%           - [OPTIONAL] testType - is a character string indicating which
%           test should be applied to the considered observations (i.e.
%           kstest, mcnemar, chi^2, binomial, wilcoxon).
%       OUTPUTS
%           - Results - a (O*(O+1))/2 structure array containing 2 fields;
%           Combination and Activity. Activity is in turn a structure array
%           of 2 elements containing Type and Pvalues. Type can either be
%           Spontaneous or Evoked, or Unaltered condition or Shuffled
%           condition.
%           - Counts - a Ox2 cell array containing the total count of
%           spikes in the given conditions and during the given time
%           windows.
% Emilio Isaias-Camacho @ GrohLab 2019
%% Validating the inputs
p = inputParser;

defaultTest = 'kstest';
validTests = {'kstest','mcnemar','chi2','binomial','wilcoxon'};
checkTest = @(x) any(validatestring(x,validTests));

checkStack = @(x) any([islogical(x), size(x) >= 1, numel(size(x)) == 3]);

[~, Nt, NTa] = size(dStack);

checkCondFlag = @(x) all([size(x,1) == NTa, islogical(x)]);
checkTimeFlags = @(x) all([size(x,1) == 2,size(x,2) == Nt, islogical(x)]);


addRequired(p, 'dStack', checkStack)
addRequired(p, 'condFlag', checkCondFlag)
addRequired(p, 'timeFlags', checkTimeFlags)
addOptional(p, 'test', defaultTest, checkTest)

p.KeepUnmatched = true;

parse(p,dStack, condFlag, timeFlags, varargin{:})
testType = p.Results.test;

%% Preparatory variables
% Number of alignment points per conditions
Na = sum(condFlag, 1);

% Test selector
switch testType
    case validTests{1} % KStest
        statFun = @kstest2;
    case validTests{2} % McNemar
        fprintf(1,'McNemar not yet implemented. Selecting kstest\n')
        statFun = @kstest2;
        testType = validTests{1};
    case validTests{3} % Chi^2
        fprintf(1,'Chi^2 not yet implemented. Selecting kstest\n')
        statFun = @kstest2;
        testType = validTests{1};
    case validTests{4} % Binomial
        fprintf(1,'Binomial not yet implemented. Selecting kstest\n')
        statFun = @kstest2;
        testType = validTests{1};
    case validTests{5} % Wilcoxon
        statFun = @signrank;
end

% Number of conditions to cycle through
Ncond = size(condFlag, 2);

% Output variable
possComb = (Ncond*(Ncond+1))/2;
Counts = cell(Ncond,2);
Results = struct('Combination',{},'Activity',struct('Type',{},'Pvalues',{}));
Results = repmat(Results,possComb,1);

%% Main loop
cr = 1;
for ccond = 1:Ncond
    sponCountA = ...
        squeeze(sum(dStack(2:end,timeFlags(1,:),condFlag(:,ccond)),2));
    evokCountA = ...
        squeeze(sum(dStack(2:end,timeFlags(2,:),condFlag(:,ccond)),2));
    Counts(ccond,:) = [{sponCountA},{evokCountA}];
    cc2 = ccond + 1;
    while cc2 <= Ncond
        % Comparing condition A versus condition B; spontaneous and evoked
        sponCountB = ...
            squeeze(sum(dStack(2:end,timeFlags(1,:),condFlag(:,cc2)),2));
        evokCountB = ...
            squeeze(sum(dStack(2:end,timeFlags(2,:),condFlag(:,cc2)),2));
        [~, P] = runStatTest(statFun, sponCountA, sponCountB);
        Results(cr).Combination = sprintf('%d %d\t%s',ccond,cc2,testType);
        Results(cr).Activity(1).Type = 'Spontaneous';
        Results(cr).Activity(1).Pvalues = P;
        [~, P] = runStatTest(statFun, evokCountA, evokCountB);
        Results(cr).Activity(2).Type = 'Evoked';
        Results(cr).Activity(2).Pvalues = P;
        cc2 = cc2 + 1;
        cr = cr + 1;
    end
    % Spontaneous vs evoked in the same condition. 
    [~, P] = runStatTest(@signrank, sponCountA, evokCountA);
    Results(cr).Combination = sprintf('%d %d\t%s',ccond,ccond,'signrank');
    Results(cr).Activity(1).Type = 'Unaltered condition';
    Results(cr).Activity(1).Pvalues = P;
    % Spontaneous vs evoked in the same condition after shuffling.
    shufSubsA = randperm(Na(ccond),Na(ccond));
    shufSubsB = randperm(Na(ccond),Na(ccond));
    [~, P] = runStatTest(@signrank,...
        sponCountA(:,shufSubsA), evokCountA(:,shufSubsB));
    Results(cr).Activity(2).Type = 'Shuffled condition';
    Results(cr).Activity(2).Pvalues = P;
    cr = cr + 1;
end

end

function [H, P] = runStatTest(statFun, obs1, obs2)
Ncl = size(obs1,1);
Out1 = zeros(Ncl,1);
Out2 = Out1;
for ccl = 1:Ncl
   [Out1(ccl), Out2(ccl)] = statFun(obs1(ccl,:), obs2(ccl,:));
end
if sum(round(Out1) - Out1)
    H = logical(Out2);
    P = Out1;
else
    H = logical(Out1);
    P = Out2;
end
end
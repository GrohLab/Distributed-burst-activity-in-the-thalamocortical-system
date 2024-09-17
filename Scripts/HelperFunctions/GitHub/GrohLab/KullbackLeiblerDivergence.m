function kbd = KullbackLeiblerDivergence(P1,P2)
% KULLBACKLEIBLERDIVERGENCE returns the similarity measure between two
% probability distributions
if abs(sum(P1) - 1.0) > 1e-5 || sum(~P1)
    P1(P1 == 0) = range(P1)*1e-6;
    P1 = P1/sum(P1);
end
if length(P1) == length(P2)
    N = length(P1);
    if abs(sum(P2) - 1) > 1e-3 || sum(~P2)
        P2(P2 == 0) = range(P2)*1e-6;
        P2 = P2/sum(P2);
    end
else
    kbd = NaN;
    fprintf('! The distributions have different resolutions (different length)!\n')
    fprintf('No divergence calculated\n')
    return;
end
aux = rdivide(P1,P2);
kbd = dot(P1,log2(aux))/log2(N);
end
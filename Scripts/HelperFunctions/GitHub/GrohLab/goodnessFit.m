function [ rsq ] = goodnessFit( pts, mdl )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Observed points
y = pts(:,2);
x = pts(:,1);
% Estimated linear model
Yhat = x*mdl(1) + mdl(2);

ss = sumsqr(y - nanmean(y));    % Sum of squares w.r.t. the original mean
ssr = sumsqr(Yhat - y);         % Sum of error
rsq = 1-(ssr/ss);               % Squared areas. R^2

end


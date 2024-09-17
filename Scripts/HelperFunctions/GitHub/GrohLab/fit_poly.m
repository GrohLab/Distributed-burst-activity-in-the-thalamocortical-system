function [ mdl, Yhat, rsq ] = fit_poly(x,y,m)
% fit_poly returns the fitted model, the curve estimation and the R^2 to
% the data. It accepts polynomials from first order up until computer
% memory and precision allows.
% [MODEL, Y_ESTIMATE, R^2] = fit_poly( X, Y, POLY_ORDER)

[x, y] = check_orientation(x,y);
if isnumeric(m)
    if size(x,2) >= m+1 && size(y,2) >= m+1
        [Yhat, mdl] = detrend_profile(m,x,y);
        rsq = goodnessFit([x',y'],mdl);
    else
        mdl = [NaN;NaN];
        Yhat = -1;
        rsq = -1;
        % display(['There must be at least ',num2str(m+1),' points to ',...
        %    'estimate a ',num2str(m),' order curve'])
    end
elseif ischar(m)
    switch m
        case 'e'
            y(y==0) = y(y==0) + eps;
            [Yhat, mdl] = detrend_profile(1,x,log(y));
            rsq = goodnessFit([x',log(y)'],mdl);
            auxA = exp(mdl(2));
        case 'log'
            [Yhat, mdl] = detrend_profile(1,x,exp(y));
            rsq = goodnessFit([x',exp(y)'],mdl);
            auxA = log(mdl(2));
        case {'sin','cos'}
            [Yhat, mdl] = detrend_profile(1,x,acos(y));
            rsq = goodnessFit([x',acos(y)'],mdl);
            auxA = cos(mdl(2));
        otherwise
            error(...
                'The option is either not yet implemented or there is a typo: %s\nOnly ''e'', ''log'', ''cos'', or  ''sin''',m)
    end
    auxB = mdl(1);
    mdl = [auxA;auxB];
end
end

function [x_out,y_out]=check_orientation(x,y)
[xx,~] = size(x);
[yx,~] = size(y);
x_out = x;
y_out = y;
if xx > 1
    x_out = x';
end
if yx > 1
    y_out = y';
end
end
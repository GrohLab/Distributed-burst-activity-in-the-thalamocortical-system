function [polft,ft] = detrend_profile(m,x,y)
%DETREND_PROFILE Returns the polynomial fit and the fit parameters for the
%given $x$ and $y$
% Uses the pseudo-inverse matrix to compute the fit.

% Avoiding undeterminations (NaN) in the computations.
x(y==-Inf)=[];
x(y==Inf)=[];
y(y==-Inf)=[];
y(y==Inf)=[];
if length(x) == length(y)
    % Polynomial with positive exponent
    
    if m >= 1
        if isrow(x)
            x = x';
        end
        M = x.^(m:-1:0);
        ft = pinv(M) * y';
        polft = zeros(size(x));
        for ex = m:-1:1
            polft = polft + ft(m-ex+1)*x.^ex;
        end
        polft = polft + ft(m+1);
    else
        % Either a root, an inverse function or a constant.
        if m < 0
            % Inverse function
            warning('Attempting a 1/x^n fit.')
            M = (x.^m)';
            ft = pinv(M) * y';
            polft = x.^ft;
        elseif(m>0)
            % Root
            warning('Fitting a root!')
            M = (x.^m)';
            ft = pinv(M) * y';
            polft = x.^ft;
        else
            % Polynomial order 0
            warning(...
                'The polynomial order 0 is the mean of the ''y'' values')
            ft = nanmean(y);
        end
    end
else
    error('The input vectors must have equal length')
end
end
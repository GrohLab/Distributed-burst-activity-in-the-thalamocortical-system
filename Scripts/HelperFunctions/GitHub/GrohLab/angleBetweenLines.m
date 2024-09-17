function angulo = angleBetweenLines(m1,m2,radDeg)
%ANGLEBETWEENLINES returns the angle in degrees between the lines given
%their slopes.
%   The function accepts only two arguments: the slope for each line and
%   returns the angle between them according to the following formula: 
%                   $\alpha = \atan{\frac{m2-m1}{1+m1\times m2}}$

tan_arg = (m2-m1)./(1+m1.*m2);
if nargin == 2
    radDeg = 'deg';
end
angulo = atan(tan_arg); 
if strcmp(radDeg,'deg') || nargin == 2
    angulo = (180/pi)*angulo;
end

end
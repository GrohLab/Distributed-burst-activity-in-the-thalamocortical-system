function [n, d] = getHesseLineForm(mdl)
%GETHESSELINEFORM computes the unique orthogonal unitary vector
%corresponding to a given line equation.
%   [n, d] = getHesseLineForm(mdl)
% INPUTS
%       mdl - 2x1 array containing the slope and y-intersection of a linear
%             equation in the form y = m·x + b
% OUTPUT
%       n - 2x1 vector coordinates
%       d - distance from the origin
% Emilio Isaias-Camacho @Neurophotonics Lab 2016
angl = angleBetweenLines(0,mdl(1),'rad');
n = [-sin(angl);...
    cos(angl)];
x = 0;
y = mdl(2);
d = [x y]*n;
end


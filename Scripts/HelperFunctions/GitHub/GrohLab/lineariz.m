function [m,b]=lineariz(average,top,buttom)
m = (top-buttom)/(max(average)-min(average));
b = top - m*max(average);
end
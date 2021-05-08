function [perc] = returnPercentile(scores,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t = 1:0.05:100;
perc = 0;
for i = 1:length(t)
    y = prctile(scores,t(i));
    if y >= x
       perc = t(i);
       break;
    end
end
end


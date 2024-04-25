function [vals,S] = powerfit(x,y)
%POWERFIT
[b,a,~,~,S]=logfit(x,y,'loglog');
a  = 10^a;
vals = [a b];
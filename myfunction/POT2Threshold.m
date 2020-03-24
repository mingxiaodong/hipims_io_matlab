function [TH_out,RT_out]=POT2Threshold(X,RT)
% find the threhold value for daily records under a given return period(year)
r = tiedrank(X);
n = numel(X);
EP = (r-0.44)/(n+0.12);%EP Gringorten (WMO, 1983)
RT_day = 1./(1-EP);
RT_year = RT_day/365.25;
[~,ind] = min(abs(RT_year-RT));
TH_out = X(ind);
RT_out = RT_year(ind);
end

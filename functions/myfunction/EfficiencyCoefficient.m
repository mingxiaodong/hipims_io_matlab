function [NSE,RMSE] = EfficiencyCoefficient(varargin)
% EfficiencyCoefficient calculate the Nash–Sutcliffe efficiency(NSE) 
% and root-mean-square error (RMSE)
%   Qsimu and Qobv are vectors withe the same length
if length(varargin) == 2
    Qobv = varargin{1};
    Qsimu = varargin{2};
elseif length(varargin) == 4
    Qobv = varargin(1:2);
    Qsimu = varargin(3:4);
else
    error('wrong number of arguments')
end
if size(Qobv,2)==2
    T_simu = Qsimu{1};
    q_simu = Qsimu{2};
    T_obv = Qobv{1};
    q_obv = Qobv{2};
    minT = max(min(T_obv(:)),min(T_simu(:)));
    maxT = min(max(T_obv(:)),max(T_simu(:)));
    ind1 = T_obv>=minT&T_obv<=maxT;
    T = T_obv(ind1);
    y_obv = q_obv(ind1);
    if isdatetime(T(1))
        T = datenum(T);
        T_simu = datenum(T_simu);        
    end
    y_simu = interp1(T_simu,q_simu,T,'linear','extrap');
else
    y_simu = Qsimu;
    y_obv = Qobv;
end

NSE = 1-sum((y_simu-y_obv).^2)/sum((y_obv-mean(y_obv)).^2);
RMSE = sqrt(sum((y_simu-y_obv).^2)/numel(y_obv));
end


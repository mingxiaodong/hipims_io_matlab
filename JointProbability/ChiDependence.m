function [chiu,chibaru,u] = ChiDependence(data,varargin)
%% [chiu,chibaru,u] = ChiDependence(data,u/nq,conf)
% data: two column matrix
nq = 100; %The number of quantiles at which the measures are evaluated
conf = 0.95;
if isempty(varargin)
    u = [];
else
    u = varargin{1};
    if length(varargin)==2
        conf = varargin{2};
    end
end


% remove nan value in pair
ind = isnan(data(:,1))|isnan(data(:,1));
data(ind,:)=[];
n = length(data);
% uniform data to [0 1]
xrank = tiedrank(data(:,1));
yrank = tiedrank(data(:,2));
data = [xrank/(n+1),yrank/(n+1)];
rowmax = max(data,[],2);
rowmin = min(data,[],2);

% qlim2 = c(min(rowmax) + eps, max(rowmin) - eps);
%   if (~isempty(qlim))
%     if (qlim(1) < qlim2(1)) 
%       error("lower quantile limit is too low")
%     elseif (qlim(2) > qlim2(2)) 
%       error("upper quantile limit is too high")
%     elseif (qlim(1) > qlim(2)) 
%       error("lower quantile limit is less than upper quantile limit")
%     else
%       qlim = qlim2;
%     end
%   end


if isempty(u)
    u = (linspace(min(rowmax)+0.00000001,max(rowmin)-0.00000001,nq))';
elseif numel(u) == 1&&u>1
    nq = u;
    u = (linspace(min(rowmax)+0.00000001,max(rowmin)-0.00000001,nq))';
else
    nq = numel(u);
    if size(u,1)==1
        u=u';
    end
end
cu = u*0;
cbaru = cu;
for i=1:nq
    cu(i) = mean(rowmax<u(i));  %P{F(X)?u,G(Y)?u}
    cbaru(i) = mean(rowmin>u(i)); %P{F(X)>u,G(Y)>u}
end
chiu = 2-log(cu)./log(u);
chibaru = 2*log(1-u)./log(cbaru)-1;
cnst = norminv((1 + conf)/2);
varchi = ((1./log(u).^2 * 1)./cu.^2 .* cu .* (1 - cu))./n;
varchi = cnst * sqrt(varchi);
varchibar = (((4 * log(1 - u).^2)./(log(cbaru).^4 .* cbaru.^2)) .* cbaru .* (1 - cbaru))./n;
varchibar = cnst * sqrt(varchibar);
chiu = [chiu-varchi,chiu,chiu+varchi];
chibaru = [chibaru-varchibar,chibaru,chibaru+varchibar];

chiulb = 2 - log(max(2 * u - 1, 0))./log(u);
chibarulb = 2 * log(1 - u)./log(1 - 2 * u + max(2 * u -1, 0)) - 1;
chiu(chiu>1)=1;
chibaru(chibaru>1)=1;
chiu = max(chiu,chiulb);
chibaru = max(chibaru,chibarulb);
end

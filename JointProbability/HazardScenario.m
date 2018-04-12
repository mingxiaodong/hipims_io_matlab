clear,clc
n = 1000;
X = wblrnd(1,0.5,[n,1]); [~,XI] = sort(X);
Y = gamrnd(1,2,[n,1]); [~,YI] = sort(Y);
Z = gamrnd(1,2,[n,1]);
%%
addpath '/Users/b4042552/Documents/MATLAB/hline_vline'
group = ones(size(X));
group(X<mean(X)&Y<mean(Y)&Z<mean(Z)) = 0;
group(X>=mean(X)&Y>=mean(Y)&Z>=mean(Z)) = 2;
gscatter(X,Y,group,'kbr','***')
legend('No hazard','Single hazard', 'Multi-hazard')
ax_h = gca;
Xrange = [0 mean(X)*2];
Yrange = [0 mean(Y)*2];
Zrange = [0 mean(Z)*2];
ax_h.XTick=Xrange;
ax_h.XTickLabel = {'',''};
ax_h.XLim = Xrange;
ax_h.YTick=Yrange;
ax_h.YTickLabel = {'',''};
ax_h.YLim = Yrange;
ax_h.ZTick=Zrange;
ax_h.ZTickLabel = {'',''};
ax_h.ZLim = Zrange;
hline(mean(Y),'r')
vline(mean(X),'r')
% print('HazardScenario','-djpeg','-r300')

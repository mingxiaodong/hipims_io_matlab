%% Example for how to plot two maps in one figure with different colormap and adjust the order of map
% created on 5 May 2017 by Xiaodong
%% basic data
[x,y,z] = peaks(100);
% data for the background map
z1 = z;
% data for the showing map
z2 = z; z2(z<=0)=nan; z2 =z2+100;
%% Method 1
figure
h1 = surface(x,y,z1);
ax = gca; 
ax.NextPlot = 'add';
ax.SortMethod = 'childOrder';
h2 = surface(x,y,z2);
numColor = 32;
cmap = [gray(numColor);parula(numColor)];

valueZ = z1;
A = (numColor-1)*(valueZ-min(valueZ(:)))/range(valueZ(:));
A = round(A);
% A = min(numColor,A);
set(h1,'CData',A)
valueZ = z2;
B = (numColor-1)*(valueZ-min(valueZ(:)))/range(valueZ(:));
B = round(B)+numColor;
set(h2,'CData',B)
colormap(ax,cmap);
caxis([min(A(:)),max(B(:))]);
cb = colorbar;
cb.Limits = [numColor+1, numColor*2];
cbTicks = str2double(cb.TickLabels)/numColor;
cbTicks = cbTicks*range(valueZ(:))+min(valueZ(:));
cb.TickLabels = num2cell(round(cbTicks));
%% Method 2
figure
% ax.NextPlot = 'replace';
ax1 = axes;
h1 = surface(x,y,z1);
ax2 = axes;
h2 = surface(x,y,z2);
linkaxes([ax1 ax2])
% ax1.Visible = 'off';
ax2.Visible = 'off';
colormap(ax1,'gray');
colormap(ax2,'parula');

cb1 = colorbar(ax1);
cb1.Location = 'westoutside';
cb2 = colorbar(ax2);
cb2.Location = 'eastoutside';
set([ax1,ax2],'Position',[0.15 0.15 0.8 0.75])
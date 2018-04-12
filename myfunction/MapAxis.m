function MapAxis(ax)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
XTickgap = 10000;
YTickgap = 10000;
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XTick = ax.XLim(1):XTickgap:ax.XLim(2);
ax.YTick = ax.YLim(1):YTickgap:ax.YLim(2);
ax.XTickLabel = string((ax.XTick-ax.XLim(1))/1000);
ax.YTickLabel = string((ax.YTick-ax.YLim(1))/1000);
end


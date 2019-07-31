function MapAxis(ax,varargin)
%MapAxis(gca,mapUnit,Tickgap,FontSize)  Summary of this function goes here
%   mapUnit: 1(m),1000(km)
%   Tickgap: scalar(m) or [XTickgap YTickgap]
%   FontSize: scalar
% Example
% MapAxis(gca,1000)
% MapAxis(gca,1000,[2000, 1000])
% MapAxis(gca,1000,[2000, 1000],9)
Tickgap = 10000;
FontSizeNew = [];
if isempty(varargin)
    mapUnit=1;
else
    mapUnit = varargin{1};
    if length(varargin)>1
        Tickgap = varargin{2};
        if length(varargin)>2
            FontSizeNew = varargin{3};
        end
    end
end
if numel(Tickgap)==1
    XTickgap = Tickgap;
    YTickgap = Tickgap;
else
    XTickgap = Tickgap(1);
    YTickgap = Tickgap(2);
end
ax.Box = 'on';
% ax.XGrid = 'on';
% ax.YGrid = 'on';
ax.XTick = ax.XLim(1):XTickgap:ax.XLim(2);
ax.YTick = ax.YLim(1):YTickgap:ax.YLim(2);
ax.XTickLabel = string((ax.XTick-ax.XLim(1))/mapUnit);
ax.YTickLabel = string((ax.YTick-ax.YLim(1))/mapUnit);
ax.YTickLabelRotation = 90;
% ax.XLimMode = 'manual';
% ax.YLimMode = 'manual';
% ax.ZLimMode = 'manual';
% mapUnitStr = 'm';
% if mapUnit==1000
%     mapUnitStr = 'km';
% end
% dim = [ax.Position(1) ax.Position(2) 0.1 0.05];
% 
% annotation('textbox',dim,'String',['Map unit: ' mapUnitStr],'FitBoxToText','on')
% txContent = ['map unit: ' mapUnitStr];
% tx = text(ax.XLim(end),ax.YLim(1),txContent);
% tx.Units = 'characters';
% tx.Position(1) = tx.Position(1)-length(txContent)-4;
% tx.Position(2) = tx.Position(2)+1;
if ~isempty(FontSizeNew)
    tx.FontSize = FontSizeNew;
    ax.FontSize = FontSizeNew;
end

end


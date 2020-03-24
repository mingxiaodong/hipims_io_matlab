function RotateColorbarTickLabel(hc,al)
% Rotate the tick label of the vertical colorbar
% RotateColorbarTickLabel(hc,al)
% hc is the handle of color bar, al is the rotation angle
% Example
% Z = peaks(10)*100;
% x11 = 1; y11 = 10; dx = 1;
% R = makerefmat(x11,y11,dx,-dx);
% mapshow(Z,R,'DisplayType','Surface')
% hc = colorbar;
% RotateColorbarTickLabel(hc,90)

% Created by Xiaodong Ming on 14 Aug 2018

hcUnits = hc.Units;
hc.Units = 'inches';
tickValues = hc.Ticks;
numLabel = numel(tickValues);
yPos = (tickValues-hc.Limits(1))/range(hc.Limits)*hc.Position(4);
TickLabels = hc.TickLabels;
% hc.TickDirection = 'out';
hc.TickLabels = [];
ax = gca;
axUnits = ax.Units;
ax.Units = 'inches';
x = hc.Position(1)-ax.Position(1)+hc.Position(3);
y = hc.Position(2)-ax.Position(2)+yPos;
for i=1:numLabel
    text(x,y(i),TickLabels(i),'Rotation',al,'Units','inches');
end
ax.Units = axUnits;
hc.Units = hcUnits;
end
%%
% ax = gca;
% axUnits = ax.Units;
% ax.Units = 'points';
% hcUnits = hc.Units;
% hc.Units = 'points';
% 
% 
% tickValues = hc.Ticks;
% numLabel = numel(tickValues);
% yPos = (tickValues-hc.Limits(1))/range(hc.Limits)*hc.Position(4);
% TickLabels = hc.TickLabels;
% % hc.TickDirection = 'out';
% hc.TickLabels = [];
% 
% x = hc.Position(1)-ax.Position(1) + hc.Position(3);%+ax.FontSize/2;%
% y = yPos;
% % hy = text(x,0,'abc','Units','pixels');
% for i=1:numLabel
%     text(x,y(i),TickLabels(i),'Rotation',90,'Units',ax.Units,'HorizontalAlignment','center');
% end
% %-ax.FontSize*0*length(TickLabels(i))
% % ax.Units = axUnits;
% % hc.Units = hcUnits;
% %

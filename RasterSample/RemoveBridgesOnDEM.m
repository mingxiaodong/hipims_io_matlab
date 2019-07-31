%% read a river DEM and create a bridge
clear,clc
[Z,R] = ArcgridreadM('fakeriver.asc');
bridgeValues = Z(:,130:135);
bridgeValues(~isnan(bridgeValues)) = max(bridgeValues(:));
Z(:,130:135) = bridgeValues;
figure('units','normalized','outerposition',[0 0 1 1])
mapshow(Z,R,'DisplayType','mesh')
view([23 26])
box on
grid on
ax = gca;
ax.BoxStyle = 'full';
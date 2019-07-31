clc,clear
addpath /Users/b4042552/Dropbox/Matlab/myfunction
load('CarlisleRiverChannelData.mat')
cellsize = 2;
errorValue = 1;
tic
bankLine1 = [channelBoundaryLine(1).X', channelBoundaryLine(1).Y'];
bankLine2 = [channelBoundaryLine(2).X', channelBoundaryLine(2).Y'];
[Z,R] = GenerateBedElevationRaster(bankLine1,bankLine2,cellsize,crossSectionLine,10);
toc
%%
figure
mapshow(Z,R,'DisplayType','Surface')
ax=gca;
ax.DataAspectRatio = [1,1,0.1];
% axis image
%%
figure;surf(flipud(Z)','EdgeColor','none')

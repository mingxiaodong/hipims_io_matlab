clc
tic
bankLine1 = [channelBoundaryLine(1).X', channelBoundaryLine(1).Y'];
bankLine2 = [channelBoundaryLine(2).X', channelBoundaryLine(2).Y'];
[Z,R] = GenerateBedElevationRaster(bankLine1,bankLine2,5,crossSectionLine,10);
toc
%%
figure;mapshow(Z,R,'DisplayType','Surface')
axis image
%%
figure;surf(flipud(Z)','EdgeColor','none')

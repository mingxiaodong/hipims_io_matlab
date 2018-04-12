% MainRoad = shaperead('MainRoad_Local.shp');
clear,clc
[Z,R] = arcgridread('h_144001.asc');
Z(Z<=0) = nan;

%%
figure
mapshow(Z,R,'DisplayType','surface')
mapshow('MainRoad_Local.shp','FaceColor','none','LineWidth',0.1);%'LineStyle','--',
mapshow('WaterSystem_Local.shp','FaceColor',[0.5 0.5 0.5]);
mapshow('Buildings_Local.shp','FaceColor',[0.1 0.1 0.1],'EdgeColor','none');
mapshow('AreaMask_Local.shp','FaceColor','none','LineWidth',2);
%% generate files for multi-gpu model
OrginalInputLocation = 'G:/cudaSWEsSolver/LondonFrame4';
NewMultiInputLocation = 'G:/cudaSWEsSolver/LondonFrame4/MultiGPU/';
NumSec = 8;
DomainDecomposite(OrginalInputLocation,NewMultiInputLocation,NumSec)

%% load normalInundatedArea.asc
[z_waterSurface,~]= arcgridread('normalInundatedArea.asc');
%%
[z_result,r]=arcgridread('results/h_max_event9.asc');
Z_plot = z_result;
Z_plot(z_waterSurface>0)=nan;
figure
mapshow(Z_plot,r,'DisplayType','Surface')
axis image
%% unzip file
clear,clc
gridShp = shaperead('C:\Users\b4042552\Google Drive\MyResearch\London\T10kmGrid_LC4.shp');
outputdir = 'F:\Data\Environment Agency (National 2m)\LC4\';
zipFileFolder = 'F:\Data\Environment Agency (National 2m)\DSM\';
for i=1:length(gridShp)
    filename = [zipFileFolder 'LIDAR-DSM-2M-' gridShp(i).TILE_NAME '.zip'];
    unzip(filename,outputdir)
end
%%
[z,r] = arcgridread('G:\cudaSWEsSolver\LondonFrame4\MultiGPU\7\input\mesh\DEM.txt');
figure
mapshow(z,r,'DisplayType','Surface')
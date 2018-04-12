%% rectify DEM data based on shapefile polygon
%% read file
clear,clc
demFileName = 'G:\cudaSWEsSolver\LondonFlood1953\dsm.asc';
shpFileName = 'G:\Data\London\DykeArea.shp';
[z_dem, r_dem] = arcgridread(demFileName);
changeArea_shp = shaperead(shpFileName);

%% plot
figure
fig_h1 = mapshow(z_dem, r_dem, 'DisplayType','Surface'); 
zdatam(fig_h1,z_dem-1000)
% demcmap(z_dem)
axis manual %off
mapshow(changeArea_shp,'FaceColor','none','EdgeColor','red') % show catchment boundary
axis image 
caxis([-20,100])


%% find the points to be changed
z_dem_RemoveDyke = z_dem;
for i = 1:length(changeArea_shp)
    pix_RowCol = map2pix(r_dem,[changeArea_shp(i).X]', [changeArea_shp(i).Y]');
    pix_RowCol(isnan(pix_RowCol(:,1)),:)=[];
    mask = poly2mask(pix_RowCol(:,2),pix_RowCol(:,1),size(z_dem,1),size(z_dem,2));
    innerPointsValue = z_dem(mask);
    innerPointsValue(innerPointsValue>median(innerPointsValue)) = median(innerPointsValue);
    z_dem_RemoveDyke(mask) = innerPointsValue;
end
%%
figure
fig_h2 = mapshow(z_dem_RemoveDyke, r_dem, 'DisplayType','Surface'); 
zdatam(fig_h2,z_dem_RemoveDyke-1000)
% demcmap(z_dem_RemoveDyke)
axis manual %off
mapshow(changeArea_shp,'FaceColor','none','EdgeColor','red') % show catchment boundary
axis image
caxis([-20,100])
%% write Arc ASCII file for new DEM
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
Arcgridwrite('London1953_RemoveDyke.asc',z_dem_RemoveDyke,r_dem)
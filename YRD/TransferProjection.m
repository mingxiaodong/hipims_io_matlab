%% Transfer Projection Coords to Geographic Coords
clear,clc
% input the county boundary of YRD
CountyBound = shaperead('H:\Data\YRD bound\YRD_county.shp');
% input the projection of YRD
proj = geotiffinfo('H:\Data\YRD bound\DEM_YRD1.tif');
% input km grid of YRD
load YRD_KmsGrid
% Transfer XY coords of county boundary to latitude and longitude
% projinv: X&Y to Lat&Lon
% projfwd: Lat&Lon to X&Y
for i = 1:140
    [CountyBound(i).Lat,CountyBound(i).Lon] =...
        projinv(proj,CountyBound(i).X,CountyBound(i).Y);    
end
% Transfer km grid of YRD to lat and lon coords
[Grid_Lat,Grid_Lon] = projinv(proj,YRD_KmsGrid(:,1),YRD_KmsGrid(:,2));
geoshow(Grid_Lat,Grid_Lon,'DisplayType','point')
geoshow(CountyBound,'DisplayType','polygon','FaceColor','none');
YRD_GeoGrid = [Grid_Lat,Grid_Lon, YRD_KmsGrid(:,3)];
[lat, lon] = projinv(proj,XkmG(:),YkmG(:));
LatG = reshape(lat,size(XkmG));
LonG = reshape(lon,size(YkmG));
save YRD_GeoGrid YRD_GeoGrid CountyBound LatG LonG CkmG


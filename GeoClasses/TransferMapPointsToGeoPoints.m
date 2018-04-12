% read dem file
[Z1,R1] = arcgridread('DEM.txt');
%%obtain XY coordinates of all raster points
R_map = R1;
dx = R_map(2,1); dy = R_map(1,2);
x11 = R_map(3,1); 
y11 = R_map(3,2);
x_se = x11+dx*size(Z1,2); %x coor in southeast point
y_se = y11+dy*size(Z1,1); %y coor in southeast point
lon11 = input('Geo lon of x11:'); %find the value in arcmap
lat11 = input('Geo lat of y11:');
lon_se = input('Geo lon of x_se:');
lat_se = input('Geo lat of y_se:');
%%
%define Rgeo

dlon = (lon_se - lon11)/size(Z1,2);
dlat = (lat_se - lat11)/size(Z1,1);
R_geo = makerefmat(lon11, lat11, dlon, dlat);
XYV = dlmread('h_172800.dat');
x = XYV(:,1); y = XYV(:,2);
[row,col] = map2pix(R_map,x,y);
[lat, lon] = pix2latlon(R_geo,row,col);
Fuzhou30mGeoCoor = table(lon,lat);
writetable(Fuzhou30mGeoCoor,'Fuzhou30mGeoCoor.txt')
%%
lon_se = lon11+dlon*size(Z2,2);
lat_se = lat11+dlat*size(Z2,1);
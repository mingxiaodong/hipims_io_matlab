function [rain_source,t_value_s,z_RainMask] = Netcdf2RainInput(ncName,z_dem,r_dem)
%NETCDF2RAININPUT generate netcdf data from UKV model to rainfall mask and
%rainfall source data of the HiPIMS model domain.
% ncName is the full name of netcdf file
% z_dem, r_dem are the elevation and reference matrix of model domain DEM
% rain_source: matrix of time series(s, 1st column) and rainfall rate(m/s, 2nd~end columns)
% z_RainMask: matrix of rainfall mask, same
% Created by Xiaodong Ming on 2017-11-8
%%
size_z_dem = size(z_dem);
source = ncName;
ncid = netcdf.open(source);
[~,nvars,~,~] = netcdf.inq(ncid);
varNames = cell(nvars,1);
varValues = cell(nvars,1);
for i=1:nvars
    varNames{i} = netcdf.inqVar(ncid,i-1);
    varValues{i} = netcdf.getVar(ncid,i-1);
%     eval([varNames{i} '=  varValues{i};'])
end
netcdf.close(ncid);
clear ncid nvars source
ind = strcmp('rlon',varNames);
rlon = double(varValues{ind});
ind = strcmp('rlat',varNames);
rlat = double(varValues{ind});
ind = strcmp('lsrain',varNames);
rainData = double(varValues{ind});
ind = strcmp('t',varNames);
t_value = double(varValues{ind});

%% convert coordinates
[LAT,LON] = meshgrid(rlat, rlon);
SP_coor = [177.5-180 -37.5];
grid_in = [LON(:),LAT(:)];
[grid_out] = rotated_grid_transform(grid_in, 2, SP_coor);
lon = grid_out(:,1);
lat = grid_out(:,2);
lon = reshape(lon,size(LON));
lat = reshape(lat,size(LAT));
[E, N] = ll2os(lat, lon);
%% create rainfall raster
[x11,y11] = pix2map(r_dem,1,1);
dx = 1500; dy = -dx;
r_rain = makerefmat(x11,y11,dx,dy);
nrows_rain = ceil(size_z_dem(1)/dx*r_dem(2));
ncols_rain = ceil(size_z_dem(2)/dx*r_dem(2));
z_rain = nan(nrows_rain,ncols_rain);
[X_rain, Y_rain] = Raster2FeaturePoints(z_rain,r_rain);

%% rain source
[t_value_u,ia,~]= unique(t_value);
% figure(1)
% t0 = datetime(2015,12,3,0,0,0);
% t_DT = t0+double(t_value_u);
t_value_s = t_value_u*24*3600; % time unit converted from day to second
rain_mat = nan(nrows_rain,ncols_rain,numel(t_value_s));
% R = makerefmat(0,0,1500,-1500);
for i=1:length(ia)
    Zdata = rainData(:,:,:,ia(i));% unit convert from (kg m-2 s-1) to (m s-1)
    F = scatteredInterpolant(E(:),N(:),Zdata(:),'nearest');
    z_rain = F(X_rain,Y_rain);
    rain_mat(:,:,i) = z_rain/1000;
%     mapshow(z_rain,r_rain,'DisplayType','Surface')
%     title(datestr(t_DT(i)))
%     pause(0.1)
end
rain_source = rain_mat;
% rain_source = [t_value_s rain_mat];
%% rain mask
[X_dem, Y_dem] = Raster2FeaturePoints(z_dem,r_dem);
Zdata = 0:numel(z_rain)-1;
F = scatteredInterpolant(X_rain(:),Y_rain(:),Zdata(:),'nearest');
z_RainMask = F(X_dem,Y_dem);

end


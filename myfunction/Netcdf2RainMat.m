function [rainMat,timeSeries] = Netcdf2RainMat(ncName,z_mask,r_mask)
%Netcdf2RainMat generate netcdf data from UKV model to rainfall matrix and
%rainfall grid reference mat, time series.
% ncName is the full name of netcdf file
% rainMat: 3D matrix of rainfall rate(m/s, 2nd~end columns)
% timeSeries: time series(s, one column vector)
% Created by Xiaodong Ming on 2017-11-14
%%
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
[X_rain, Y_rain] = Raster2FeaturePoints(z_mask,r_mask);

%% rain source
[t_value_u,ia,~]= unique(t_value);
t_value_s = t_value_u*24*3600; % time unit converted from day to second
rainMat = nan(size(z_mask,1),size(z_mask,2),numel(t_value_s));
for i=1:length(ia)
    Zdata = rainData(:,:,:,ia(i));
    F = scatteredInterpolant(E(:),N(:),Zdata(:),'nearest');
    z_rain = F(X_rain,Y_rain);
    rainMat(:,:,i) = z_rain/1000;% unit convert from (kg m-2 s-1) to (m s-1)
end
timeSeries = t_value_s;
end

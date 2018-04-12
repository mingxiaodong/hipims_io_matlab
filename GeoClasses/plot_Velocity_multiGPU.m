%% velocity files
clear,clc
t = 21600;
z_U = 0;
z_h = 0;
%%
for n=7:-1:0
%*****file names
DEMName = ['G:\cudaSWEsSolver\Eden\MultiGPU\' num2str(n) '\input\mesh\DEM.txt'];
outputFolder = ['G:\cudaSWEsSolver\Eden\MultiGPU\' num2str(n) '\output\'];
cd(outputFolder)
fileName_h  = ['h_'  num2str(t) '.dat'];
fileName_hU = ['hU_' num2str(t) '.dat'];
%*****read files
[z_dem,r_dem] = arcgridread(DEMName);
    h_raster = nan(size(z_dem)); 
    Ux_raster = nan(size(z_dem)); 
    Uy_raster = nan(size(z_dem));
data = dlmread(fileName_h); 
    x = data(:,1);
    y = data(:,2);    
    h = data(:,3);
data = dlmread(fileName_hU); 
    Ux = data(:,3); Ux = Ux./(h+1e-6); 
    Uy = data(:,4); Uy = Uy./(h+1e-6);  clear data  
%%*****convert feature points to raster
    [rows,cols] = map2pix(r_dem,x,y); rows=round(rows); cols=round(cols);
    rows(rows>size(z_dem,1)) = size(z_dem,1); rows(rows<1) = 1;
    cols(cols>size(z_dem,2)) = size(z_dem,2); cols(cols<1) = 1;
    ind = sub2ind(size(h_raster),rows,cols); clear rows cols
h_raster(ind) = h;
Ux_raster(ind) = Ux; 
Uy_raster(ind) = Uy; clear ind
U = (Ux_raster.^2+Uy_raster.^2).^0.5;
if n==7
    z_U=U; r_U = r_dem;
    z_h=h_raster; 
else
    z_temp = [z_U;U(3:end,:)];
    z_U = z_temp;
    z_temp = [z_h;h_raster(3:end,:)];
    z_h = z_temp;
end
end
%%
% [M,I] = max(z_h(:));
I = find(z_h(:)>8);
[a,b] = ind2sub(size(z_U),I);
[x,y] = pix2map(r_U,a,b);
%%
z_new = z_h;
% z_new(z_new<0.1)=nan;
figure
hold on
scatter(x,y,'r*')
mapshow(z_new-100,r_U,'DisplayType','Surface');
axis image
hold off
%% *****plot velocity map
figure
quiver(x,y,Ux,Uy,3)
axis equal
title(['Velocity:' num2str(t)])

% figure;mapshow(result_h,R,'DisplayType','Surface')
clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\RainfallUK
cd('G:\Data\UK_EU\UK_Rainfall_Radar')
%% decompress the *gz.tar files
fileList = dir('*201512*gz.tar');
outputdir = 'G:\Data\UK_EU\UK_Rainfall_Radar\201512';
for i = 1:length(fileList)
    untar(fileList(i).name,outputdir)
end
%% decompress the *gz files
cd(outputdir)
files = '*dat.gz';
gunzip(files,outputdir)
delete(files)
%% get radar rain files to be extracted
mylist = dir('*1km-composite.dat');
%mylist = mylist(1:288);
%% load the catchment frame
% basinPoly = shaperead('G:\Data\London\LondonCatchment2.shp');
basinPoly = shaperead('G:\Data\Cumbria\EdenCatchment.shp');
x_west = basinPoly.BoundingBox(1);
x_east = basinPoly.BoundingBox(2);
y_south = basinPoly.BoundingBox(3);
y_north = basinPoly.BoundingBox(4);
% m = ceil(range(basinPoly.Y)/1000); %number of rows
% n = ceil(range(basinPoly.X)/1000); %number of columns
% x = basinPoly.X(1:end-1);
% y = basinPoly.Y(1:end-1);
% x = (x-basinPoly.BoundingBox(1,1))/1000;
% y = (y-basinPoly.BoundingBox(1,2))/1000;
% basinMask = poly2mask(x,y,m,n);
%% 
for i = 1:length(mylist)
    % read the binary file
    fileName = mylist(i).name;
    [int_gen_hd, rl_gen_hd, rl_datsp_hd, char_hd, int_datsp_hd, ...
        rr_dat_mat] = rdnim1km(fileName);
    % the Z data of the whole UK
    Z = rr_dat_mat;
    Z(Z<=0)=0;
    Z = Z/32;
    % make reference
    x_upleftcor = rl_gen_hd(5);
    y_upleftcor = rl_gen_hd(3);
    dx = rl_gen_hd(6);
    dy = -rl_gen_hd(4);
    x11 = x_upleftcor+0.5*dx; %half dx east of the upper left corner
    y11 = y_upleftcor+0.5*dy; %half dy south of the upper left corner
    R = makerefmat(x11, y11, dx, dy);
    % extract Z data based on catchment frame
    [row,col] = map2pix(R,[x_west,x_east],[y_south y_north]);
    row = round(row); row_vec = min(row)-1:max(row)+1;
    col = round(col); col_vec = min(col)-1:max(col)+1;
    [col_mat,row_mat] = meshgrid(col_vec,row_vec);
    extractInd = sub2ind(size(rr_dat_mat),row_mat(:),col_mat(:));
    Z_extracted = Z(extractInd);
    Z_extracted = reshape(Z_extracted,size(row_mat));
    % write asc file
    x0_arc = R(3,1)+R(2,1)/2 + min(col_vec)*dx;
    y0_arc = R(3,2)+R(1,2)*size(Z,1)+R(1,2)/2 - (size(Z,1)-max(row_vec))*dy;
    cellsize_arc = abs(dx);
    noDATA_value = -9999;
    Z_extracted(isnan(Z_extracted)) = noDATA_value;
    [nrows, ncols] = size(Z_extracted);
    writeFileName = [fileName(32:43) '.asc'];
    fileID = fopen(writeFileName,'w');
    fprintf(fileID,'ncols    %d\n', ncols);
    fprintf(fileID,'nrows    %d\n', nrows);
    fprintf(fileID,'xllcorner    %.2f\n', x0_arc);
    fprintf(fileID,'yllcorner    %.2f\n', y0_arc);
    fprintf(fileID,'cellsize    %.2f\n', abs(R(1,2)));
    fprintf(fileID,'NODATA_value    %d\n', noDATA_value);
    fclose(fileID);
    dlmwrite(writeFileName,Z_extracted,'-append','delimiter','\t')
end
%% make a rainfall mask file based on 1km cell
Z_mask = 1:numel(Z_extracted);
Z_mask = reshape(Z_mask,size(Z_extracted));
writeFileName = 'rainfall_mask.asc';
fileID = fopen(writeFileName,'w');
fprintf(fileID,'ncols    %d\n', ncols);
fprintf(fileID,'nrows    %d\n', nrows);
fprintf(fileID,'xllcorner    %.2f\n', x0_arc);
fprintf(fileID,'yllcorner    %.2f\n', y0_arc);
fprintf(fileID,'cellsize    %.2f\n', abs(R(1,2)));
fprintf(fileID,'NODATA_value    %d\n', noDATA_value);
fclose(fileID);
dlmwrite(writeFileName,Z_mask,'-append','delimiter','\t')
%% make animation for radar rainfall
clear,clc
mylist = dir('201412*asc');
mylist = mylist(289:576);
basinPoly = shaperead('G:\Data\London\LondonCatchment2.shp');
%%
for i=1:length(mylist)
h_fig = figure;
[Z,R] = arcgridread(mylist(i).name);
Z(Z==0) = nan;
mapshow(Z,R,'DisplayType','surface')
hold on
mapshow(basinPoly,'FaceColor','none')
hold off
colorbar;
caxis([0 15])
axis equal
title(mylist(i).name(1:12))
print(h_fig,'-djpeg','-r100',[mylist(i).name(1:12) 'h.jpg']);
close(h_fig);
end
%% make animations
clear,clc
cd('G:\Data\Cumbria\UK_Rainfall_Radar\Composite_1km')
mylist = dir('20151204*.jpg');
animationName = '20151204.gif';
delayTime = 0.3;
N = size(mylist);
for i = 1:N
    img = imread(mylist(i).name);
    imshow(img);
    frame=getframe(gcf);
    im=frame2im(frame);
    [A,map]=rgb2ind(img,256);
    if i == 1
		imwrite(A,map,animationName,'gif','LoopCount',Inf,'DelayTime',delayTime);
	else
		imwrite(A,map,animationName,'gif','WriteMode','append','DelayTime',delayTime);
    end
end
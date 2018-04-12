clear,clc
addpath 'C:/Users/b4042552/Documents/MATLAB/freezeColors'
CaseName = 'Fuzhou';
curFolder = ['F:\PhD\ResultsAndVisualization\' CaseName '\'];
output_path = [curFolder 'output\'];
cd(output_path)
%% read files of dem, h, hU
dT = 900;
[Z_DEM,R] = arcgridread([curFolder 'input\mesh\dem.txt']);% R colume1: [0,dx,x11], colume2: [dy,0,y11]
filename_z = [output_path 'z_0.dat']; result_z = dlmread(filename_z);
x = result_z(:,1); y = result_z(:,2); z = result_z(:,3);
% Gauges = load([curFolder 'input\field\gauges_pos.dat']);
%%transfer XY coordinates to index according to matrix Z
nrows = size(Z_DEM,1); ncols = size(Z_DEM,2);
[row,col] = map2pix(R,x,y); 
row = ceil(row); col = ceil(col); row(row==0)=1; col(col==0)=1;
row(row>nrows)= nrows; col(col>ncols) = ncols;
index = round(sub2ind(size(Z_DEM),row,col));
%%
for i = 11:40
T = (i-1)*dT;
%filename_hU = [output_path 'hU_' num2str(T) '.dat']; result_hU = dlmread(filename_hU);
%filename_eta = [output_path 'eta_' num2str(T) '.dat']; result_eta = dlmread(filename_eta);
filename_h = [output_path 'h_' num2str(T) '.dat']; result_h = dlmread(filename_h);
h = result_h(:,3); 
h_grid = nan(size(Z_DEM)); 
h_grid(index) = h;
h_grid(h_grid <= 0.01) = nan;
%%*plot result map
h_fig = figure;
h0 = mapshow(Z_DEM,R,'DisplayType','surface');demcmap(Z_DEM); freezeColors; 
zdatam(h0,Z_DEM-800)
hold on
Mapping_Result = h_grid; MapTitle = 'Water Depth'; %Mapping_Result(Mapping_Result==0) = nan;
h1 = mapshow(Mapping_Result,R,'DisplayType','surface');
colormap(jet); colorbar; 
%caxis([min(Mapping_Result(:)) max(Mapping_Result(:))]); 
title([MapTitle ' (Time = ' num2str(T/3600) 'h)']); xlabel('Meter towards East'); ylabel('Meter towards North');
% h2 = plot(Gauges(:,1),Gauges(:,2),'r*');zdatam(h2,max(Mapping_Result(:)))
set(gca,'layer','top')
hold off;
axis equal;
caxis([0 2]);
set(gca,'FontSize',10); 
set(gcf, 'PaperPositionMode', 'auto');
print(h_fig,'-djpeg','-r600',['floodmap' num2str(i) '_t' num2str(T/3600) 'h.jpg']);
pause(0.5);
close(h_fig);
end
%% make animation
cd('C:\Users\b4042552\Dropbox\ShareWithXilin\fuzhou_animation')
list = dir('floodmap*.jpg');
N = size(list);
for i = 1:N
    img = imread(list(i).name);
    imshow(img);
    frame=getframe(gcf);
    im=frame2im(frame);
    [A,map]=rgb2ind(img,256);
    if i == 1
		imwrite(A,map,'animation.gif','gif','LoopCount',Inf,'DelayTime',0.5);
	else
		imwrite(A,map,'animation.gif','gif','WriteMode','append','DelayTime',0.5);
    end
end
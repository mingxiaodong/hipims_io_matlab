clear,clc
addpath 'C:/Users/b4042552/Documents/MATLAB/freezeColors'
CaseName = 'Fuzhou30m';
curFolder = 'G:\cudaSWEsSolver\LondonSeabasin\20mNoDyke';%['C:/Users/b4042552/Google Drive/Data/' CaseName '/']; %current folder
output_path = 'G:\cudaSWEsSolver\LondonSeabasin\20mNoDyke';%[curFolder 'output\'];
cd(output_path)
[Z_DEM,R] = arcgridread([curFolder '\input\mesh\dem.txt']);
%%
Z_DEM = Z;
result_z = dlmread('z_0.dat');
x = result_z(:,1); y = result_z(:,2); z = result_z(:,3); clear result_z
% gauge_t_h = load('h_gauges.dat');
% gauge_t_hU = load('hU_gauges.dat');
%%*transfer XY coordinates to index according to matrix Z
nrows = size(Z_DEM,1); ncols = size(Z_DEM,2);
[row,col] = map2pix(R,x,y); 
row = round(row); col = round(col);
row(row==0) = 1; col(col==0) = 1;
row(row>nrows)= nrows; col(col>ncols) = ncols; clear nrows ncols
index = round(sub2ind(size(Z_DEM),row,col)); clear row col
%% read files of h, hU
T = 1200;
h_Or_eta = 'h';
filename = ['hU_' num2str(T) '.dat']; result_hU = dlmread(filename); clear filename
hu = result_hU(:,3); hv = result_hU(:,4); clear result_hU
filename = [h_Or_eta '_' num2str(T) '.dat'];
result_h_eta = dlmread(filename);clear filename
h_eta = result_h_eta(:,3); clear result_h_eta
if h_Or_eta=='h'
    h = h_eta; eta = h + z; clear h_eta
else
    eta = h_eta; h = h_eta-z; clear h_eta
end
u = hu./(h + 0.000000001);
v = hv./(h + 0.000000001);
vel = sqrt(u.*u + v.*v); clear u v
%%Data for plotting
% eta_grid = nan(size(Z_DEM)); eta_grid(index) = eta;
h_grid = nan(size(Z_DEM));h_grid(index) = h;
vel_grid = nan(size(Z_DEM)); vel_grid(index) = vel;

%% plot result map
% T = 262800;
% [h_grid,R] = arcgridread(['h_max_' num2str(T) '.asc']);
figure;
%h0 = mapshow(Z_DEM,R,'DisplayType','surface');demcmap(Z_DEM); freezeColors; zdatam(h0,Z_DEM-100)
hold on
% Mapping_Result = h_grid; MapTitle = 'Water Depth'; Mapping_Result(Mapping_Result<=0.1) = nan;
% Mapping_Result = h_grid+Z_DEM; MapTitle = 'Water Elevation'; 
Mapping_Result = vel_grid; MapTitle = 'Velocity';
h1 = mapshow(Mapping_Result,R,'DisplayType','surface');
colormap(parula); colorbar; 
caxis([min(Mapping_Result(:)) max(Mapping_Result(:))]); 
% caxis([0 20])
title([MapTitle ' (Time = ' num2str(T) 's)']); xlabel('Km towards East'); ylabel('Km towards North');
% h2 = plot(Gauges(:,1),Gauges(:,2),'r*');zdatam(h2,max(Mapping_Result(:)))
set(gca,'layer','top')
hold off
axis image;
box on; grid on
%% show the largest value
% [M,I] = max(vel);
I = vel(:)>30;
zdatam(h1,Mapping_Result-800)
hold on, plot(x(I),y(I),'*r'), hold off; 
%caxis([0 50]);
%% show gauges position
Gauges = load([curFolder 'input\field\gauges_pos.dat']);
for i = 1:length(Gauges)
    text(Gauges(i,1)+50,Gauges(i,2),max(Mapping_Result(:)),num2str(i));
end
set(gca,'FontSize',20,'XlimMode','auto','YlimMode','auto'); 
%% plot 'h_gauges.dat' SHOULD RUN THE FIRST SECTION IN ADVANCE!
gauge_t_h = load('h_gauges.dat');
t0 = datetime('2016-9-27 15:00','InputFormat','yyyy-MM-dd HH:mm','Format','yyyy-MM-dd HH:mm');
% TimeHis = gauge_t_h(:,1)/3600/24+t0; x_unit = 'h';
TimeHis = gauge_t_h(:,1)/3600; x_unit = 'h';
gauge_h = gauge_t_h(:,2:end); y_unit = 'm';
%% gauge_h = gauge_h; y_unit='cm';
ShowPoints = 35;%1:size(gauge_h,2);
figure
legendName = cell(1,size(gauge_h,2));
for i=1:size(gauge_h,2)
    legendName{i} = ['Gauge ' num2str(i)];
end
set(groot,'defaultAxesColorOrder',[1 0 0; 0 0.8 0; 0 0 1; 0.8 0.8 0; 1 0 0.8],...
    'defaultAxesLineStyleOrder','-|--|-.|:')
plot(TimeHis,gauge_h(:,ShowPoints),'LineWidth',1.2)
set(groot,'defaultAxesLineStyleOrder','remove');set(groot,'defaultAxesColorOrder','remove')
legend(legendName{ShowPoints},'Location','best');
xlabel(['Time (' x_unit ')']); ylabel(['Water depth (' y_unit ')'])
%% plot 'hu_gauges.dat' SHOULD RUN THE FIRST SECTION IN ADVANCE!
gauge_t_hU = load('hU_gauges.dat');
gauge_hu = gauge_t_hU(:,2:2:end-1);
gauge_hv = gauge_t_hU(:,3:2:end);
gauge_hU = (gauge_hu.^2 + gauge_hv.^2).^0.5;

figure
set(groot,'defaultAxesColorOrder',[1 0 0;0 1 0;0 0 1;0 1 1; 1 0 1],'defaultAxesLineStyleOrder','-|--|-.|:');
plot(TimeHis(1:length(gauge_hU)),gauge_hU(:,ShowPoints),'LineWidth',1.2)
legend(legendName{ShowPoints});
xlabel(['Time (' x_unit ')']); ylabel('Velocity (m/s)')
set(groot,'defaultAxesColorOrder','remove');set(groot,'defaultAxesLineStyleOrder','remove');
%% plot timestep data
timestep = dlmread('timestep_log.txt');
startind = 1;
figure, plot(timestep(startind:end,1),timestep(startind:end,2))
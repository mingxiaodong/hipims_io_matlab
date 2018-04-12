%% import result files and save as mat variable
%% *load basic files
clear,clc
addpath 'C:/Users/b4042552/Documents/MATLAB/freezeColors'
CaseName = 'Fuzhou';
curFolder = ['G:\cudaSWEsSolver\' CaseName '\']; %current folder
output_path = [curFolder 'output\'];
cd(output_path)
[Z_DEM,R] = arcgridread([curFolder '\input\mesh\dem.txt']);
result_z = dlmread('z_0.dat');
x = result_z(:,1); y = result_z(:,2); z = result_z(:,3); clear result_z
gauge_t_h = load('h_gauges.dat');
gauge_t_hU = load('hU_gauges.dat');
%%*transfer XY coordinates to index according to matrix Z
nrows = size(Z_DEM,1); ncols = size(Z_DEM,2);
[row,col] = map2pix(R,x,y); 
row = ceil(row); col = ceil(col);
row(row==0) = 1; col(col==0) = 1;
row(row>nrows)= nrows; col(col>ncols) = ncols; clear nrows ncols
index = round(sub2ind(size(Z_DEM),row,col)); clear row col
%% *get the file names and time series
h_Or_eta = 'h';
list_h_Eta = dir([h_Or_eta '_*0.dat']);
TimeSeries = zeros(length(list_h_Eta),1);
for i = 1:length(TimeSeries)
    TimeSeries(i) = str2double(list_h_Eta(i).name(length(h_Or_eta)+2:end-4));
end
[TimeSeries, sortID] = sort(TimeSeries);
list_h_Eta = list_h_Eta(sortID);
% list_hU = dir('hU_*0.dat'); 
% list_hU = list_hU(sortID);

%% load result files(h only)
ResultCell_H = cell(length(TimeSeries),1);
h_Or_eta = 'h';
tic
for i = 1:length(TimeSeries)
    T = TimeSeries(i);
    filename = [h_Or_eta '_' num2str(T) '.dat'];
    result_h_eta = dlmread(filename);clear filename
    h_eta = result_h_eta(:,3); clear result_h_eta
    if h_Or_eta=='h'
        h = h_eta; eta = h + z; clear h_eta
    else
        eta = h_eta; h = h_eta-z; clear h_eta
    end    
    h_grid = nan(size(Z_DEM));h_grid(index) = h;    
    ResultCell_H{i} = h_grid;
    disp(T)
end
toc
savePath = 'C:/Users/b4042552/Google Drive/data';%C:\Users\b4042552\Google Drive\Data
save([savePath '/ResultCell_H72'],'ResultCell_H','TimeSeries', 'R', 'Z_DEM','-v7.3')
%% load result files(h and hU)
ResultCell_H = cell(length(TimeSeries),1);
ResultCell_vel = ResultCell_H;
tic
for i = 1:length(TimeSeries)    
    result_hU = dlmread(list_hU.name);
    hu = result_hU(:,3); hv = result_hU(:,4); clear result_hU    
    
    result_h_eta = dlmread(list_h_Eta.name);
    h_eta = result_h_eta(:,3); clear result_h_eta
    if h_Or_eta=='h'
        h = h_eta; eta = h + z; clear h_eta
    else
        eta = h_eta; h = h_eta-z; clear h_eta
    end    
    u = hu./(h + 0.000000001);
    v = hv./(h + 0.000000001);
    vel = sqrt(u.*u + v.*v);    
    %%Data for plotting
%     eta_grid = nan(size(Z_DEM)); eta_grid(index) = eta;
    h_grid = nan(size(Z_DEM));h_grid(index) = h;
    vel_grid = nan(size(Z_DEM)); vel_grid(index) = vel;        
    ResultCell_H{i} = h_grid;
    ResultCell_vel{i} = vel_grid;
    disp(T)
end
toc
% save([curFolder '\Result_24h'],...
%   'ResultCell_H','ResultCell_vel','gauge_t_h','gauge_t_hU','Z_DEM','R','TimeSeries','-v7.3')
clear,clc
cd('F:\LondonFlood\cudaSWEsSolver\')%\output
fileNo = 2;
outpu_dt = 3600;
t = outpu_dt*(fileNo-1);
RunType = 'cudaOverlandflowSolver';%'cudaSWEsSolver';% SWEsSolver1st
curF = cd; %current folder
result_filepath = [curF '\output_nohr\'];
result_h = load([result_filepath 'h_' num2str(fileNo,'%03u') '_' num2str(t) '.dat']);
result_hU = load([result_filepath 'hU_' num2str(fileNo,'%03u') '_' num2str(t) '.dat']);
x = result_h(:,1);
y = result_h(:,2);
h = result_h(:,3);
u = result_hU(:,3)./(h + 0.000001);
v = result_hU(:,4)./(h + 0.000001);
vel = sqrt(u.*u + v.*v);

DEMfilepath = [curF '\input\mesh\dem.txt'];
[Z,R] = arcgridread(DEMfilepath);% R colume1: [0,dx,x11], colume2: [dy,0,y11]
nrows = size(Z,1);
ncols = size(Z,2);
dx = abs(R(2,1)); %size of the grid
dy = abs(R(1,2));


%%transfer XY coordinates to index according to matrix Z
[row,col] = map2pix(R,x,y); row = round(row); col = round(col);
row(row==0)=1; col(col==0)=1;
row(row>nrows)=nrows; col(col>ncols)=ncols;
index = round(sub2ind(size(Z),row,col));
% data for mapping
Z_DEM = Z;
Z_Result = nan(size(Z)); 
Z_Result(index) = h; ReName = 'Height';%Z_DEM(index)=h;
%Z_Result(index) = vel; ReName = 'Velocity';%Z_DEM(index)=vel;
Z_Result(Z_Result==0)=nan;

%%plot model results
figure;
addpath 'C:\Users\b4042552\freezeColors'
%mapshow(Z_DEM,R,'DisplayType','surface');demcmap(Z);freezeColors; 
hold on
mapshow(Z_Result,R,'DisplayType','surface');
colorbar; colormap parula
%caxis([0,1])
%caxis([min(Z_Result(:)) max(Z_Result(:))]); 
title([ReName ' (t = ' num2str(t) 's)--' RunType])
hold off

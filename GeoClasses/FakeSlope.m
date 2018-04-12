%% Generate DEM
clear,clc
dx = 10;
dy = -dx;
slopeX = 0.05;
slopeY = 0.02;
xlength = 800*cos(atan(slopeX))*2+20;
ylength = 1000*cos(atan(slopeY));
ncols = round(xlength/abs(dx));
nrows = round(ylength/abs(dy));
x11 = 0;
y11 = 1000;
R = makerefmat(x11, y11, [0 dx], [dy 0]);
H = zeros(nrows,ncols);
[Rows, Cols] = meshgrid(1:nrows,1:ncols);
[X,Y] =  pix2map(R, Rows(:), Cols(:));

for i=1:numel(X)
    x = X(i); y = Y(i);
    if x <= 800*cos(atan(slopeX))
        z = (800*cos(atan(slopeX))-x)*slopeX + 20;%y*slopeY;
    elseif x <= 800*cos(atan(slopeX))+20
        z = y*slopeY;
    else
        z = (x-800*cos(atan(slopeX))-20)*slopeX + 20;%y*slopeY;
    end
    H(Rows(i),Cols(i)) = z;
end
%3D surface, move it down
mapshow(H,R,'DisplayType','surface'); demcmap(H);axis normal;  axis equal; 
xlabel('x(m)');ylabel('y(m)');zlabel('z(m)')
%%
DEM_path = 'G:\cudaSWEsSolver\FakeSlope\input\mesh\';
fileID = fopen([DEM_path 'DEM.txt'],'w');
fprintf(fileID,'ncols    %d\n', ncols);
fprintf(fileID,'nrows    %d\n', nrows);
fprintf(fileID,'xllcorner    %f\n', 0.0);
fprintf(fileID,'yllcorner    %f\n', 0.0);
fprintf(fileID,'cellsize    %f\n', abs(dx));
fprintf(fileID,'NODATA_value    -9999\n');
fclose(fileID);
dlmwrite([DEM_path 'DEM.txt'],H,'-append','delimiter','\t')

%% *Generate cell ID
[Z,R] = arcgridread([DEM_path 'DEM.txt']);
IO_Bound = [810 1015; 820 1000];
[Valid_ID,Bound_ID] = gen_VB_ID(Z,R,IO_Bound);
%% Generate field files
writePath = 'G:\cudaSWEsSolver\FakeSlope\input\field\';
%*write z.dat; h.dat; hU.dat;
fileNames = {'z','h','hU','precipitation'};
I_h = Z*0;
I_h(Bound_ID==2)=0.001;
InitialValue = {Z,I_h,[Z*0, Z*0],Z*0};
% boundary type vector for [h; hU]* [bound1, bound2]
BT_vec = {[2 0 0; 2 2 0],[2 0 0; 2 1 0]};
for i = 1:4;
    WriteInitialValue(writePath,fileNames{i},Valid_ID,Bound_ID,InitialValue{i},BT_vec)
end

%*precipitation_mask.data
Pre_Mask = [Valid_ID(:) zeros(numel(Valid_ID),1)];
Pre_Mask(Pre_Mask(:,1)==-1,:)=[];
Pre_Mask = sortrows(Pre_Mask,1);

%*print file: precipitation_mask.dat
filename = [writePath 'precipitation_mask.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'$Element Number\n'); fprintf(fileID,'%d\n',length(Pre_Mask));
fprintf(fileID,'$Element_id  Value\n'); fprintf(fileID,'%-12d %d\n',Pre_Mask');
fclose(fileID);

%*source data
Ts = (0:600:12000)';
uq = Ts*0; 
vq = uq;
h_Tide = Ts*0; %water depth
pre = Ts*0; % precipitation m/s
pre(Ts<=5400)= 10.8/1000/3600;
filename = [writePath 'h_BC_0.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'%-10d %.8f \n',[Ts h_Tide]');
fclose(fileID);
filename = [writePath 'hU_BC_0.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'%-10d %.8f %.8f\n',[Ts uq vq]');
fclose(fileID);
filename = [writePath 'precipitation_source_0.dat']; %precipitation_source_0.dat
fileID = fopen(filename,'w');
fprintf(fileID,'%-10d %.12f \n',[Ts pre]');
fclose(fileID);
%% gauges_pos.dat
gaugeXY = [1 810 9.801;
    2 840 59.8];
filename = [writePath 'gauges_pos.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'%.1f %.1f\n',gaugeXY(1:end-1,2:3)');
fprintf(fileID,'%.1f %.1f',gaugeXY(end,2:3)');
fclose(fileID);
type(filename)
%% load gauge value
clear,clc
gaugeValue = dlmread('h_gauges.dat');
TimeHis = gaugeValue(:,1);
gauge_h = gaugeValue(:,2:end);
gaugeValue = dlmread('hU_gauges.dat');
gauge_hu = gaugeValue(:,2:2:end-1);
gauge_hv = gaugeValue(:,3:2:end);
% gauge_UV = (gauge_hu.^2 + gauge_hv.^2).^0.5;
% gauge_Q = gauge_UV(:,2).*gauge_h(:,2)*1000;
gauge_hu1 = sum(gauge_hu(:,3:end),2); %gauge in the vertical line
gauge_hu2 = sum(gauge_hu(:,1:2),2); % gauge at the two output cells
gauge_hv1 = sum(gauge_hv(:,3:end),2); %gauge in the vertical line
gauge_hv2 = sum(gauge_hv(:,1:2),2); % gauge at the two output cells
gaugeValue = dlmread('hU_NS_gauges.dat');
gauge_hUN = gaugeValue(:,2:2:end-1);
gauge_hUS = gaugeValue(:,3:2:end);
Fs = sum(gauge_hUS(:,1:2),2)*-1;
TotalHU = ( gauge_hv2)*(-10);
gaugeValue = dlmread('hU_EW_gauges.dat');
gauge_hUE = gaugeValue(:,2:2:end-1);
gauge_hUW = gaugeValue(:,3:2:end);
TotalHU_Hill = sum(gauge_hUE(:,3:end),2)*10;
%%
figure
plot(TimeHis,TotalHU,'LineWidth',1.2)
xlabel('Time (s)');ylabel('hU (m^2/s)')


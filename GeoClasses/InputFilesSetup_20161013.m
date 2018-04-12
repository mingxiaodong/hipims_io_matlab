clear,clc
addpath('C:\Users\b4042552\Dropbox\Matlab\GeoClasses');
% addpath('/Users/b4042552/Dropbox/Matlab/GeoClasses');
CaseName = 'UrbanDistrict';
CaseFolder = ['G:\cudaSWEsSolver\' CaseName '\'];
writePath = [CaseFolder 'input/field/'];
% cd(CaseFolder);
h_Or_eta = 'h';

%% ***read initial conditions from xls file
Initial_XLS = [CaseFolder 'InitialData.xlsx'];
%WaterDepth = xlsread(Initial_XLS,'Depth','A:B');
Discharge = xlsread(Initial_XLS,'Discharge','A:B'); %m^3/s
%Rain = xlsread(Initial_XLS,'Rainfall','A:B'); % precipitation m/s
gaugeXY = xlsread(Initial_XLS,'Gauge','A:C');
IO_Bound = xlsread(Initial_XLS,'BoundFrame','B:E');% the frame range of boundary cells
BoundType = xlsread(Initial_XLS,'BoundType');
%% ***read DEM file
[Z,R] = arcgridread([CaseFolder 'input/mesh/DEM.txt']);

%***Get ID of valid cells and bound cells 
[Valid_ID,Bound_ID] = gen_VB_ID(Z,R,IO_Bound);%
% figure;subplot(2,1,1); mapshow(Z,R,'DisplayType','surface'); subplot(2,1,2); mapshow(Bound_ID,R,'DisplayType','mesh');

%% *****write z.dat; h.dat(eta.dat); hU.dat; precipitation.dat; manning.dat
Manning = 0.035;
Z_Initial = Z; %initial Z 
initialwaterdepth = arcgridread([CaseFolder '/initial_h.txt']);
h_Or_eta_Initial =  initialwaterdepth; %initial h water depth or eta water surface elevation  
hU_Initial = [zeros(size(Z)) zeros(size(Z))]; % initial hU discharge per unit width
Pre_Initial = zeros(size(Z)); % initial precipitation
Manning_Initial = Manning+ zeros(size(Z)); %initial manning

I_cumul_depth = Z_Initial*0;
I_hydraulic_conductivity = Z_Initial*0;
I_capillary_head = Z_Initial*0 + 0.06;
I_water_content_diff = Z_Initial*0 + 0.4;
sewer_sink = arcgridread([CaseFolder '/streetmask.txt']);
sewer_sink(isnan(sewer_sink)) = 0;
I_sewer_sink = sewer_sink*0.000006667;
InitialValue = { Z_Initial, h_Or_eta_Initial, hU_Initial, Pre_Initial, Manning_Initial,...
    I_cumul_depth, I_hydraulic_conductivity, I_capillary_head, I_water_content_diff, I_sewer_sink }; %initial values
fileNames = {'z',h_Or_eta,'hU','precipitation','manning',...
    'cumulative_depth','hydraulic_conductivity','capillary_head','water_content_diff',...
    'sewer_sink'};
% sewer_sink
% boundary type vector for [h; hU]* [bound1, bound2, bound3]

%%invoke the writing function
for i = 1:length(InitialValue)
    WriteInitialValue(writePath,fileNames{i},Valid_ID,Bound_ID,InitialValue{i},BoundType)
    disp([fileNames{i} '.dat has been written'])
end

%% ***Water Depth or Water Surface elevation
% h_BC_0.dat (eta_BC_0.dat)
TD0 = [0 0; 3600 0];
filename = [writePath h_Or_eta '_BC_0.dat']; fileID = fopen(filename,'w');
fprintf(fileID,'%-10.1f %15.10f \n',TD0'); fclose(fileID);
% h_BC_1.dat (eta_BC_0.dat)
load('TideDepth.mat')
filename = [writePath h_Or_eta '_BC_1.dat']; fileID = fopen(filename,'w');
fprintf(fileID,'%-10.1f %15.10f \n',TideDepth'); fclose(fileID);

%% ***Discharge
% hU_BC_0.dat
filename = [writePath 'hU_BC_0.dat']; fileID = fopen(filename,'w');
Ts_Q = Discharge(:,1); hu = Discharge(:,2)*0; hv = hu; 
% BoundWidth = abs(R(2))*sum(Bound_ID(:)==2); %m 
% calculate unit width discharge in X and Y direction
% if BoundWidth > 0
%     Q_per_m = Discharge(:,2)/BoundWidth; % m^2/s 
%     [I,J] = ind2sub(size(Bound_ID),find(Bound_ID(:)==2)); %range(cols):delta_x, range(a):delta_y
%     theta = atan(range(J)/range(I)); %theta: the angle between bound line and the vertical line    
%     hu = Q_per_m*cos(theta); hv = Q_per_m*sin(theta); %m^2/s
% else hu = Ts_Q*0; hv = Ts_Q*0; % no inflow bound, width is zero
% end
fprintf(fileID,'%-10.1f %15.10f %15.10f\n',[Ts_Q hu hv]'); fclose(fileID);

%% ***Precipitation Mask 
% precipitation_mask.data
[RainMask,~] = arcgridread([CaseFolder 'rainmask.txt']);
VID_RMask = [Valid_ID(:) RainMask(:)]; % zeros(numel(Valid_ID),1) rainfall mask value is 0 
VID_RMask(isnan(Valid_ID(:))|Valid_ID(:)==-1,:)=[]; % delete nodata value
VID_RMask = sortrows(VID_RMask,1);
% print file
filename = [writePath 'precipitation_mask.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'$Element Number\n'); fprintf(fileID,'%d\n',length(VID_RMask));
fprintf(fileID,'$Element_id  Value\n'); fprintf(fileID,'%-12d %d\n',VID_RMask');
fclose(fileID);
%% *Rainfall Source
% precipitation_source_0.dat
load([CaseFolder 'RainRecords72.mat'])
RainSource = cell(1,41);
NumOfHour = 72;
T = 0:3600:3600*(NumOfHour-1);
for i =1:41
    T_Rain = [T; RainRecords(i,3:end)/3600/1000]';
    T_Rain(isnan(T_Rain(:,2)),:)=[];
    RainSource{i} = T_Rain;
end
%%
for s = 0:length(RainSource)-1    
    filename = [writePath 'precipitation_source_' num2str(s) '.dat']; %precipitation_source_0.dat
    fileID = fopen(filename,'w');
    T_Rain = RainSource{s+1};
    fprintf(fileID,'%-10.1f %18.14f\n',T_Rain');
    fclose(fileID);
end
%% *****gauges_pos.dat
%coordinate of gauge
% gaugeXY = xlsread(Initial_XLS,'Gauge','A:C');
filename = [writePath 'gauges_pos.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'%.3f %.3f\n',gaugeXY(1:end-1,2:3)');
fprintf(fileID,'%.3f %.3f',gaugeXY(end,2:3));
fclose(fileID);
type(filename)
clc,clear
addpath /Users/b4042552/Dropbox/Matlab/GeoClasses % add the path of FieldSetup fucntion
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
Location = 'F:\30M';
CaseFolder = [Location '\0\'];
cd(CaseFolder)
[Z,R] = arcgridread([CaseFolder 'input/mesh/DEM.txt']);
% Initial_XLS = [Location 'InputData.xlsx'];
% IO_BoundFrame = xlsread(Initial_XLS,'BoundFrame','B:E');% the frame range of boundary cells
% BoundCode = xlsread(Initial_XLS,'BoundType');
% GaugeCoor = xlsread(Initial_XLS,'Gauge','B:C');
% initial_h = arcgridread('initial_h.txt');
% initial_eta = initial_h+Z;
% StreetMask = arcgridread([CaseFolder 'streetmask.txt']);
% RainMask = arcgridread([CaseFolder 'rainmask.txt']);
% RainSource = load([CaseFolder 'RainSource_41p_72h.mat']); RainSource = RainSource.RainSource;
% load([CaseFolder 'TideDepth.mat'])
% TideDepth = xlsread(Initial_XLS,'Depth','A:B');
% Discharge = xlsread(Initial_XLS,'Discharge','A:B');
%% prepare function parameters
h_Eta = 'eta'; % write h file or eta file
h_BC_0 = TideDepth;
% h_BC_1 = [TideDepth(:,1),TideDepth(:,2)-7.72];
hU_BC_0 = Discharge;
h_BC_Source = {h_BC_0};
hU_BC_Source = {hU_BC_0};
initial_hE_hU_pre = {initial_h, 0, 0};
sewer_sink = 0;%0.000006667*StreetMask;
hydro_params_Value = {0.015,sewer_sink,0,0,0,0};

%% invoke function
FieldSetup(CaseFolder, Z, R,'h_Eta',h_Eta,...
    'initial_hE_hU_pre', initial_hE_hU_pre, ...
    'BoundCode',BoundCode,'IO_BoundFrame',IO_BoundFrame,'h_BC_Source',h_BC_Source,'hU_BC_Source',hU_BC_Source,...
    'GaugeCoor',GaugeCoor,'WriteAllFiles',true);
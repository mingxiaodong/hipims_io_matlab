%% create field files for TyneCoast
clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
CaseFolder = 'G:\cudaSWEsSolver\TyneCoast\';
cd(CaseFolder)
load('FluvialEvents.mat') %boundary conditon for two flood events in Dec 2015
load('TideSurge20131205.mat') %boundary conditon for tide surge in Dec 2013
load('JointExtreme.mat') %boundary conditon for Joint extreme flow and tide
DEMName = 'dsm.txt';
%*****prepare boundary condition
%**1.event 2015-12-03
% Flow = Bywell15minFlow031215071215.Value_m3s; Tide = NorthShields15minTidal031215071215.Value_m;
% times_setup = [64800 345600 3600 345600];
% times_setup = [0 199800 199800 199800];
% Tide = 0.55+NorthShields15minTidal031215071215.Value_m;
% Flow = Tide*0+100;

%**2.event 2015-12-25
% Flow = Bywell15minFlow251215271215.Value_m3s; Tide = NorthShields15minTidal251215271215.Value_m;
% times_setup = [64800 226800 3600 226800];

%**3.event Tide surge Dec 2013
% Flow = TideSurge20131205.Flow_Bywell; Tide = TideSurge20131205.Tide_NorthShield;
% times_setup = [0 117900 3600 117900];

%**4.event Joint extreme
% Flow = JointExtreme.Flow; Tide = JointExtreme.Tide;
% times_setup = [0 85500 3600 85500];

%**4.event Future Joint extreme
Flow = JointExtreme.Flow; Tide = JointExtreme.Tide+0.55;
times_setup = [0 85500 3600 85500];


TimeSeries = (0:15*60:15*60*(length(Flow)-1))';
figure; plotyy(TimeSeries/3600,Flow,TimeSeries/3600,Tide)
h_BC_Source = {[TimeSeries Tide]};
hU_BC_Source = {[TimeSeries Flow]};
t_range = (times_setup(2)-times_setup(1))/3600;
disp(['t_range:' num2str(t_range) 'h'])
%%
dlmwrite([CaseFolder '\input\times_setup.dat'],times_setup,' ')
h_Eta = 'eta'; % write h file or eta file
BoundCode = [2 0 0,	2 0 0, 3 0 0;...
             2 2 0, 3 0 0, 2 1 0];
IO_BoundFrame = [411474	564334	411489	564241;
                 436320	571809	443239	560280];
GaugesPosition = [418537	563627;
                  422313	563013;
                  426633	563848;
                  432127	565960;
                  437175	568468;
                  412148	564682;
                  436419    568449];
% EventDuration = Bywell15minFlow031215071215.Timestamp(end)-Bywell15minFlow031215071215.Timestamp(1);
% EventDuration = duration(EventDuration,'Format','s');
% XlsName = 'TyneCoast_InitialData.xlsx';
% initialwaterdepth = arcgridread([CaseFolder '/initial_h.asc']);
% landuse = arcgridread([CaseFolder '/landuse.asc']);
% StreetMask = arcgridread([CaseFolder '/streetmask.asc']);
% RainMask = arcgridread([CaseFolder '/rainmask.asc']);
[Z,R] = arcgridread(DEMName);
copyfile(DEMName,[CaseFolder '\input\mesh\DEM.txt'])
Initial_h = arcgridread('initial_h_futureSeaLevel.asc');
h_Or_eta_Initial = Z+Initial_h; h_Or_eta_Initial(Z<0)=0;
h_Or_eta_Initial(isnan(h_Or_eta_Initial)) = 0;
initial_hE_hU_pre = {h_Or_eta_Initial, 0, 0};
%% write all files
FieldSetup(CaseFolder, Z, R,'h_Eta',h_Eta,...
        'initial_hE_hU_pre',initial_hE_hU_pre,... %Initial values of h, hU and rainfall
        'BoundCode',BoundCode,'IO_BoundFrame',IO_BoundFrame,... % boundary code
        'h_BC_Source',h_BC_Source,'hU_BC_Source',hU_BC_Source,... % boundary source
        'GaugeCoor',GaugesPosition,'WriteAllFiles','true');
%% ONLY write boundary condition files
FieldSetup(CaseFolder, Z, R,'h_Eta',h_Eta,...
        'BoundCode',BoundCode,'IO_BoundFrame',IO_BoundFrame,... % boundary code
        'h_BC_Source',h_BC_Source,'hU_BC_Source',hU_BC_Source); % boundary source
%%
FieldSetup(CaseFolder, Z, R, 'GaugeCoor',GaugesPosition)
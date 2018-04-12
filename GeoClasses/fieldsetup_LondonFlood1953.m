%% prepare function parameters
clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
CaseFolder = 'G:\cudaSWEsSolver\LondonFlood1953';
cd(CaseFolder)
% load KingstonFlow20140205
% load GaugesPositionUpperLee

%%
h_Eta = 'eta'; % write h file or eta file
DEMName = 'London1953_RemoveDyke.asc';
[Z,R] = arcgridread(DEMName);
copyfile(DEMName,[CaseFolder '/input/mesh/DEM.txt'])
times_setup = [100000 213000 7200 43200];
dlmwrite([CaseFolder '/input/times_setup.dat'],times_setup,' ')
Depth = xlsread('LondonFlood1953_InitialData.xlsx','Depth','A:B');
IO_BoundFrame = [597717 199748 598313 164228];
BoundCode = [2 0 0, 3 0 0;...
             2 2 0, 2 1 0];
h_BC_Source = {Depth};
% hU_BC_Source = {[0 0; 60 0],KingstonFlow20140205};
h_initial = Z; h_initial(Z>0)=0; h_initial(isnan(Z))=0; h_initial = -h_initial;
initial_hE_hU_pre = {h_initial,0,0};
GaugesPosition = xlsread('LondonFlood1953_InitialData.xlsx','Gauge','B:C');
%%
FieldSetupV2(CaseFolder, Z, R,'GaugeCoor',GaugesPosition,'WriteAllFiles','true',...
        'initial_hE_hU_pre',initial_hE_hU_pre, 'h_Eta',h_Eta,...
        'BoundCode',BoundCode,'IO_BoundFrame',IO_BoundFrame,... % boundary code
        'h_BC_Source',h_BC_Source); % boundary source
       
%% convert to multi_GPU files
OrginalInputLocation = 'G:\cudaSWEsSolver\London38001';
NewMultiInputLocation = 'G:\cudaSWEsSolver\London38001\MultiGPU\';
NumSec = 8;
DomainDecomposite(OrginalInputLocation,NewMultiInputLocation,NumSec)
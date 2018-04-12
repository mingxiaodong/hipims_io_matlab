%% prepare function parameters
clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
CaseFolder = 'G:\cudaSWEsSolver\Handan';
cd(CaseFolder)
%%
% RainSource = 'precipitation_source_all.dat';
h_Eta = 'h'; % write h file or eta file
DEMName = 'dem.asc';
[Z,R] = arcgridread(DEMName);
copyfile(DEMName,[CaseFolder '/input/mesh/DEM.txt'])
% times_setup = [0 82800 3600 43200];
% dlmwrite([CaseFolder '/input/times_setup.dat'],times_setup, ' ')
RainMask = arcgridread([CaseFolder '/rainfall_mask.txt']);
BoundCode = [3 0 0;...
             3 0 0];
h_BC_Source = {[0 0; 60 0]};
%%
FieldSetup(CaseFolder, Z, R,'h_Eta',h_Eta,...
        'BoundCode',BoundCode,... % boundary code
        'h_BC_Source',h_BC_Source,... % boundary source
        'RainMask',RainMask,... %rainfall mask
        'WriteAllFiles','true');
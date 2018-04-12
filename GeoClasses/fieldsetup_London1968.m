%% prepare function parameters
clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
% addpath /Users/b4042552/Dropbox/Matlab/GeoClasses
CaseFolder = 'G:\cudaSWEsSolver\LondonUperLee';
cd(CaseFolder)
%%
h_Eta = 'h'; % write h file or eta file
DEMName = 'dsm.asc';
[Z,R] = arcgridread(DEMName);
copyfile(DEMName,[CaseFolder '/input/mesh/DEM.txt'])
times_setup = [0 345600 7200 43200];
dlmwrite([CaseFolder '/input/times_setup.dat'],times_setup,' ')
RainMask = arcgridread([CaseFolder '/rainfall_mask_radar.asc']);
IO_Bound_Frame = [517630.946  169204.047 517790.596  169241.799];
BoundCode = [3 0 0;...
             3 0 0];
h_BC_Source = {[0 0; 60 0]};
h_initial = Z;
h_initial(Z>0)=0; h_initial(isnan(Z))=0; h_initial = -h_initial;
%%
initial_hE_hU_pre = {h_initial,0,0};

FieldSetup(CaseFolder, Z, R,'h_Eta',h_Eta,'initial_hE_hU_pre',initial_hE_hU_pre,...
        'BoundCode',BoundCode,... % boundary code
        'h_BC_Source',h_BC_Source,... % boundary source%'RainMask',RainMask,... %rainfall mask
        'WriteAllFiles','true');
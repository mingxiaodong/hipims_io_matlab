%% prepare function parameters
clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
% addpath /Users/b4042552/Dropbox/Matlab/GeoClasses
CaseFolder = 'G:\cudaSWEsSolver\Eden20m';
load('C:/Users/b4042552/Google Drive/Data/RainSourceFuzhou3Events')
% CaseFolder = '/Users/b4042552/Google Drive/Data/data_test/20M/';
cd(CaseFolder)
RainSource = RainSource911;
h_Eta = 'h'; % write h file or eta file
DEMName = 'dsm.txt';
[Z,R] = arcgridread(DEMName);
copyfile(DEMName,[CaseFolder '/input/mesh/DEM.txt'])
times_setup = [0 43200 3600 43201];
dlmwrite([CaseFolder '/input/times_setup.dat'],times_setup)
% initialwaterdepth = arcgridread([CaseFolder '/initialwaterdepth.txt']);
initialwaterdepth = arcgridread([CaseFolder '/h_21600.asc']);
initialwaterdepth(isnan(initialwaterdepth))=0;
landuse = arcgridread([CaseFolder '/landuse.asc']);
StreetMask = arcgridread([CaseFolder '/streetmask.asc']);
RainMask = arcgridread([CaseFolder '/rainmask.asc']);

BoundCode = [3 0 0 3 0 1;...
             3 0 0 2 1 0];
IO_BoundFrame = [436633.09	2880713.06	436856.93	2880882.92];
if strcmp(h_Eta,'h')
    h_Or_eta_Initial = initialwaterdepth;
else
    h_Or_eta_Initial = initialwaterdepth + Z;
end
h_BC_Source = {[0 0; 60 0],[0 8; 60 8]};
landIndicate = ones(size(landuse)); % 0:impermeable 1:permeable
landIndicate(ismember(landuse,[0 1 2 3])) = 0; 
I_cumul_depth = landIndicate*0.00001;
I_hydraulic_conductivity = landIndicate*0.0000038;
I_capillary_head = landIndicate*0.0889; I_capillary_head(landIndicate==0)=0.06; clear landIndicate
I_water_content_diff = 0.4;

I_sewer_sink = StreetMask;
I_sewer_sink(isnan(StreetMask))=0.000002667;
I_sewer_sink(StreetMask==2)    =0.000006667;


initial_hE_hU_pre = {h_Or_eta_Initial, 0, 0}; %clear h_Or_eta_Initial
hydro_params_Value = {0.035,I_sewer_sink,I_cumul_depth,...
    I_hydraulic_conductivity,I_capillary_head,I_water_content_diff};
% clear I_sewer_sink I_cumul_depth I_hydraulic_conductivity I_capillary_head I_water_content_diff
%%
FieldSetup(CaseFolder, Z, R,'h_Eta',h_Eta,...
        'initial_hE_hU_pre',initial_hE_hU_pre,... %Initial values of h, hU and rainfall
        'hydro_params_Value',hydro_params_Value,... % hydro parameter files
        'BoundCode',BoundCode,'IO_BoundFrame',IO_BoundFrame,... % boundary code
        'h_BC_Source',h_BC_Source,... % boundary source
        'RainMask',RainMask,'RainSource',RainSource,... %rainfall mask
        'WriteAllFiles','true');
%%
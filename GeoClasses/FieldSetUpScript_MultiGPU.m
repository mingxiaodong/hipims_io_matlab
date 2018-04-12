%%
clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
load( 'C:\Users\b4042552\Google Drive\Data\RainSourceFuzhou3Events')
RainSource = RainSource911;
CaseFolder = 'C:/Users/b4042552/Google Drive/Data/30m';
cd(CaseFolder)
load SectionRowNumber
h_Eta = 'h'; % write h file or eta file
[Z,R] = arcgridread([CaseFolder '/dem.txt']);
initialwaterdepth = arcgridread([CaseFolder '/initial_h.asc']);
landuse=arcgridread([CaseFolder '/landuse.asc']);
StreetMask = arcgridread([CaseFolder '/streetmask.asc']);
RainMask = arcgridread([CaseFolder '/rainmask.asc']);
BoundCode = [3 0 0 3 0 1; 3 0 0 2 1 0];
IO_BoundFrame = [436633.09	2880713.06	436856.93	2880882.92];
if strcmp(h_Eta,'h')
    h_Or_eta_Initial = initialwaterdepth;
else
    h_Or_eta_Initial = initialwaterdepth + Z;
end
h_Or_eta_Initial(isnan(h_Or_eta_Initial))=0;
h_BC_Source = {[0 0; 60 0],[0 8; 60 8]};
landIndicate = ones(size(landuse)); % 0:impermeable 1:permeable
landIndicate(ismember(landuse,[0 1 2 3])) = 0; 
I_cumul_depth = landIndicate*0.00001;
I_hydraulic_conductivity = landIndicate*0.0000038;
I_capillary_head = landIndicate*0.0889; I_capillary_head(landIndicate==0)=0.06;
I_water_content_diff = 0.4;
I_sewer_sink = StreetMask*0.000004+0.000002667; I_sewer_sink(isnan(I_sewer_sink))=0;

%% prepare function parameters

NumSec = length(SectionRowNumber);
for i = 1:NumSec
    SubCaseFolder = [CaseFolder '/' num2str(i-1)];
    [Z_sec,R_sec] = arcgridread([SubCaseFolder '/input/mesh/DEM.txt' ]);
    IndSec = SectionRowNumber{i};
    initial_hE_hU_pre = {h_Or_eta_Initial(IndSec,:), 0, 0};
    hydro_params_Value = {0.035,...
        I_sewer_sink(IndSec,:),...
        I_cumul_depth(IndSec,:),...
        I_hydraulic_conductivity(IndSec,:),...
        I_capillary_head(IndSec,:),...
        I_water_content_diff};
    RainMask_sec = RainMask(IndSec,:);
%%invoke function
    if i==1
        BoundOption = 'top';
    elseif i==NumSec
        BoundOption = 'bottom';
    else
        BoundOption = 'all';
    end
    FieldSetup(SubCaseFolder, Z_sec, R_sec,'h_Eta',h_Eta,...
        'initial_hE_hU_pre',initial_hE_hU_pre,... %Initial values of h, hU and rainfall
        'hydro_params_Value',hydro_params_Value,... % hydro parameter files
        'BoundCode',BoundCode,'IO_BoundFrame',IO_BoundFrame,... % boundary code
        'h_BC_Source',h_BC_Source,... % boundary source
        'RainMask',RainMask_sec,'RainSource',RainSource,... %rainfall mask
        'ExportSharedBound',BoundOption,'WriteAllFiles','true');
end


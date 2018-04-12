%% prepare function parameters
clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
CaseFolder = 'G:\cudaSWEsSolver\CylinderWater\';
% load('C:\Users\b4042552\Google Drive\Data\RainSourceFuzhou3Events')
% RainSource = RainSource911;
h_Eta = 'h'; % write h file or eta file
[Z,R] = arcgridread([CaseFolder '\input\mesh\dem.txt']);
initialwaterdepth = arcgridread('initial_h.txt');
%landuse = arcgridread([CaseFolder '/landuse.asc']);
%StreetMask = arcgridread([CaseFolder '/SewerPara1.asc']);
%RainMask = arcgridread([CaseFolder '/rainmask.asc']);
cd(CaseFolder)
initial_hE_hU_pre = {initialwaterdepth,{0,0},0};
hydro_params_Value = {0,0,0,0,0,0,};
%%
FieldSetup(CaseFolder, Z, R,'initial_hE_hU_pre',initial_hE_hU_pre,...
    'hydro_params_Value',hydro_params_Value,'WriteAllFiles','true');
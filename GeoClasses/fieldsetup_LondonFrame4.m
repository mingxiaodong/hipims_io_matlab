%% prepare function parameters
clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
% addpath /Users/b4042552/Dropbox/Matlab/GeoClasses
CaseFolder = 'G:\cudaSWEsSolver\LondonFrame4';
load C:\Users\b4042552\Dropbox\Matlab\London\multiEventData.mat
cd(CaseFolder)
% DEMName = 'dsm.asc';
% LevelBound(DEMName, [588441 174599 594262 184075] ,-20)
DEMName = 'dsm_boundLevelled.asc';
[Z,R] = arcgridread(DEMName);
%%cut DEM
domainFrame = [514014 184553;...
    594270  171554;...
    513768  167813];
[z_New,r_New] = ClipDEM(Z,R,domainFrame);
times_setup = [0 345600 3600 43200];
dlmwrite([CaseFolder '/input/times_setup.dat'],times_setup,' ')
landuse = arcgridread([CaseFolder '/landuse.asc']);
gaugeCoor = dlmread('gauges_pos.dat');
load landuseCode0805
%% set manning based on landuse
% landuseCode = readtable('C:\Users\b4042552\Google Drive\MyResearch\London\LCM2015_code.txt');

manning_Z = zeros(size(Z))+0.025;
for i=1:height(landuseCode)
    lia = ismember(landuse,landuseCode.Value(i));
    manning_Z(lia) = landuseCode.manning(i);
end
[manning_Z,~] = ClipDEM(manning_Z,R,domainFrame);
%%boundary conditions
Z = z_New; R = r_New;
h_Eta = 'h'; % write h file or eta file
IO_Bound_Frame = [517614 169204 517797 169250;...
                  588441 174599 594262 184075 ];
BoundCode = [3 0 0, 2 0 0, 3 0 1;... %1. outline bound 2. river flow in
             3 0 0, 3 0 1, 2 1 0]; %3. tide bound
tide = [multiEventData(2).tide(:,1) multiEventData(2).tide(:,2)+20];
flow = multiEventData(2).flow;
h_BC_Source = {[0 0; 60 0],tide};
hU_BC_Source = {[0 0; 60 0],flow};
%%*initial conditions
h_initial = Z;h_initial(Z>0)=0; h_initial(isnan(Z))=0; h_initial = -h_initial;
%%
InputSetup(CaseFolder, Z, R,'h_Eta',h_Eta,...
        'IO_BoundFrame',IO_Bound_Frame,'BoundCode',BoundCode,... % boundary code
        'initial_hE',h_initial,...
        'h_BC_Source',h_BC_Source,'hU_BC_Source',hU_BC_Source,... % boundary source
        'manning',manning_Z,...
        'GaugeCoor',gaugeCoor,...
        'WriteAllFiles','true');
%%
InputSetup(CaseFolder, Z, R,'manning',manning_Z);
%%
[z_h,r_h] = arcgridread('h_86400.asc');
eta_all = dlmread('eta_gauges.dat');
z_h(z_h<=0) = nan;
%%
figure
h_map = mapshow(z_h,R,'DisplayType','Surface');
zdatam(h_map,Z-500)
mapshow(gaugeCoor(:,1),gaugeCoor(:,2),'DisplayType','Point')
textStr = num2cell((1:length(gaugeCoor))');
text(gaugeCoor(:,1),gaugeCoor(:,2),textStr)
axis image
%%
figure; plot(eta_all(:,1)/3600,eta_all(:,2:10))
legendStr = strsplit(num2str(1:9,'%02d\n'),'\n');
legend(legendStr)
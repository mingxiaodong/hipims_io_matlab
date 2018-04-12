% clear,clc
% addpath('/Users/b4042552/Dropbox/Matlab/GeoClasses');
% CaseFolder = 'G:\cudaSWEsSolver\Fuzhou30m\';
[Z,R] = arcgridread([CaseFolder '\input\mesh\DEM.txt']);
[Valid_ID,Bound_ID] = gen_VB_ID(Z,R);%
ind1 = find(Valid_ID>=0);
[a,b] = ind2sub(size(Z),ind1);
ab = [a,b];
ab = sortrows(ab,[-1 2]);
ind2 =  sub2ind(size(Z),ab(:,1),ab(:,2));
%%
filename = [CaseFolder '\input\field\hU.dat'];
NumOfCells = dlmread(filename,'',[1,0,1,0]);
validID_FieldValue = dlmread(filename,'',[3,0,NumOfCells+2,2]);
Z_value = nan(size(Z));
Z_value(ind2) = ((validID_FieldValue(:,2).^2+validID_FieldValue(:,3).^2).^0.5;
% Z_value(Z_value-Z==0)=nan;
figure;mapshow(Z_value,R,'DisplayType','surface')
%% find the initial water points
H1 = arcgridread([CaseFolder '\initial_h.asc']);
H2 = arcgridread([CaseFolder '\h_21600.asc']);
ind = H1>0&H2>0;
validID_initial_h = Valid_ID(ind);
validID_initial_h = sort(validID_initial_h);
validID_initial_h(isnan(validID_initial_h)) = [];
Geo_initial_h = Coor30mGeo(validID_initial_h+1,:);
%% rectify the value of field file: e.g. z.dat
clc
filename = 'hU.dat';
%*****read
fileID = fopen(filename);
if filename(2) == 'U'
    formatSpec1 = '%-12d %f %f\n';
else
    formatSpec1 = '%-12d %f\n';
end
NumOfValidCells = textscan(fileID,'%d','HeaderLines',1,'delimiter','');
NumOfValidCells = NumOfValidCells{1};
CellID_InitialValue = textscan(fileID,formatSpec1,'HeaderLines',1,'delimiter',''); 
CellID_OfInitialValue_read = CellID_InitialValue{1};
InitialValue_read = cell2mat(CellID_InitialValue(2:end));
formatSpec2 = '%-12d %d %d %d\n';
NumOfBoundCells = textscan(fileID,'%d','HeaderLines',1,'delimiter',''); NumOfBoundCells = NumOfBoundCells{1};
CellID_BoundVectors = textscan(fileID,formatSpec2,'HeaderLines',1,'delimiter','');
CellID_OfBoundVectors_read = CellID_BoundVectors{1};
BoundVectors_read = cell2mat(CellID_BoundVectors(2:end));
fclose(fileID);
%% *****write
fileID = fopen(filename,'w');
% print valid cell and their initial value
fprintf(fileID,'$Element Number\n'); 
fprintf(fileID,'%d\n',NumOfValidCells);
fprintf(fileID,'$Element_id  Value\n');
fprintf(fileID,formatSpec1,[CellID_OfInitialValue_read InitialValue_read]');
% print boundary cell and their vector
fprintf(fileID,'$Boundary Numbers\n');
fprintf(fileID,'%d\n',NumOfBoundCells);
fprintf(fileID,'$Element_id Boundary_type\n');
fprintf(fileID,formatSpec2,[CellID_OfBoundVectors_read BoundVectors_read]');
fclose(fileID);

%% read precipitation source data
RainCell = cell(26,1);
figure
hold on
for i = 1:26
    filename = [CaseFolder '\input\field\precipitation_source_' num2str(i-1) '.dat'];
    Rain = dlmread(filename);
    dlmwrite(filename,[Rain(:,1),Rain(:,2)*0])
    RainCell{i} = Rain;
    plot(Rain(:,1),Rain(:,2))
end
hold off
%%
clear,clc
% addpath('/Users/b4042552/Dropbox/Matlab/GeoClasses');
CaseFolder = 'G:\cudaSWEsSolver\Fuzhou10m\';
[Z,R] = arcgridread([CaseFolder '\input\mesh\DEM.txt']);
h_result = arcgridread([CaseFolder '\output\h_max_46800.asc']);
h_initial = arcgridread([CaseFolder '/initial_h.asc']);
ind = h_initial>0;
h_extra = h_result; h_extra(ind)=0;
% figure;mapshow(h_extra,R,'DisplayType','surface')
%% write txt file for bed initial value
X0_arc = R(3,1)+R(2,1)/2; 
Y0_arc = R(3,2)+R(1,2)*size(Z,1)+R(1,2)/2;
nrows = size(Z,1);
ncols = size(Z,2);
cellsize = abs(R(2));
NODATA_value = -9999;
Z_arc = Z_DEM;
Z_arc(isnan(Z_arc)) = NODATA_value;
filename = 'dsm.asc';
fileID = fopen(filename,'w');
fprintf(fileID,'ncols    %d\n', ncols);
fprintf(fileID,'nrows    %d\n', nrows);
fprintf(fileID,'xllcorner    %f\n', X0_arc);
fprintf(fileID,'yllcorner    %f\n', Y0_arc);
fprintf(fileID,'cellsize    %f\n', cellsize);
fprintf(fileID,'NODATA_value    -9999\n');
% fprintf(fileID,[repmat(' %.2f ',1,ncols) '\n'],Z_arc');
dlmwrite(filename,Z_arc,'-append','delimiter','\t')
fclose(fileID);

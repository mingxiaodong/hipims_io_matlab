clear,clc
NumSec = 5; %number of section
[Z,R] =  arcgridread('/Users/b4042552/Google Drive/Data/30m/input/mesh/dem.txt');
ValidCell_Ind = ~isnan(Z); % 1: valid cell 0: invalid cell
NumValidCell = sum(ValidCell_Ind(:)); % number of valid cells
NumOfValidCellRows = sum(ValidCell_Ind,2);
CumNumOfValidRows = cumsum(NumOfValidCellRows,'reverse'); % cumulative number of valid cell in each row
CumNumValidCell = (1:NumSec)/(NumSec)*NumValidCell; %cumulative number of valid cell in each section
ind = repmat(CumNumOfValidRows,[1,NumSec])>repmat(CumNumValidCell,size(CumNumOfValidRows));
SecRow = sum(ind); SecRow = [size(Z,1) SecRow];
SectionRowNumber = cell(NumSec,1);
for i = 1:NumSec
    SectionRowNumber{i} = SecRow(i+1)+1:SecRow(i);
end
% save SectionRowNumber SectionRowNumber
%% write subDEM
parentFolder = 'F:\30M\';
% figure
% hold on
X0_arc = R(3,1)+R(2,1)/2; 
Y0_arc = R(3,2)+R(1,2)*size(Z,1)+R(1,2)/2;
NODATA_value = -9999;
BoundaryRow = cell(2,NumSec);
for i = 1:NumSec    
    Z_sub = Z(SectionRowNumber{i},:);
    nrows = size(Z_sub,1);
    ncols = size(Z_sub,2);
    ind1 = ~isnan(Z_sub);
    disp(sum(ind1(:)))
    Z_sub(isnan(Z_sub)) = NODATA_value;
    X0_arc_sec = X0_arc;
    Y0_arc_sec = Y0_arc - R(1,2)*(size(Z,1)-SectionRowNumber{i}(end));
    
    folderName = [num2str(i-1) '\input\mesh\'];
%     mkdir(parentFolder,folderName) % create mesh folder
%     mkdir([parentFolder num2str(i-1) '\input\field']) % create field folder
     DEMName = [parentFolder folderName 'DEM.txt'];
%     fileID = fopen(DEMName,'w');
%     fprintf(fileID,'ncols    %d\n', ncols);
%     fprintf(fileID,'nrows    %d\n', nrows);
%     fprintf(fileID,'xllcorner    %.2f\n', X0_arc_sec);
%     fprintf(fileID,'yllcorner    %.2f\n', Y0_arc_sec);
%     fprintf(fileID,'cellsize    %.2f\n', abs(R(1,2)));
%     fprintf(fileID,'NODATA_value    %d\n', NODATA_value);
%     dlmwrite(DEMName,Z_sub,'-append','delimiter','\t')
%     fclose(fileID);
    [Z_sec, R_sec] = arcgridread(DEMName);
    %invoke function
    [Valid_Cell_ID,Bound_Cell_ID] = gen_VB_ID(Z_sec,R_sec);
    
    top_row = Valid_Cell_ID(1,:);
    top_row(isnan(top_row))=[];
    bottom_row = Valid_Cell_ID(end,:);
    bottom_row(isnan(bottom_row))=[];
    
    BoundaryRow{1,i} = top_row;
    BoundaryRow{2,i} = bottom_row;
%     mapshow(Z_sec,R_sec,'DisplayType','surface')
end
% hold off
% axis equal
%% generate shared_boundary.dat
filename = [parentFolder 'shared_boundary.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'No. of Domains\n');
fprintf(fileID,'%d \n',NumSec);
for i = 1:NumSec
    %write section number
    fprintf(fileID,'#%d\n',i-1);
    %write bottom shared_boundary cell ID
    if i==1
        fprintf(fileID,'%d\n',[]); fprintf(fileID,'%d\n',[]); % no bottom, write two empty row
    else
        fprintf(fileID,'%d ',BoundaryRow{2,i});fprintf(fileID,'%d\n',[]);%bottom_self
        fprintf(fileID,'%d ',BoundaryRow{1,i-1});fprintf(fileID,'%d\n',[]); %bottom_neighbour
    end
    %write top shared_boundary cell ID
    if i==NumSec
        fprintf(fileID,'%d\n',[]); % no top, write two empty row
    else        
        fprintf(fileID,'%d ',BoundaryRow{1,i});fprintf(fileID,'%d\n',[]); %top_self
        fprintf(fileID,'%d ',BoundaryRow{2,i+1});fprintf(fileID,'%d\n',[]); %top_neighbour       
    end    
end
fclose(fileID);
function DomainDecomposite(orginalInput,multiInput,numGPU,varargin)
% DOMAINDECOMPOSITE  convert iuput files of single GPU model to input files
% for multi-GPU modelling
%   DomainDecomposite(orginalInput,multiInput,numGPU)
%   input files converted from single model to multi-GPU model.
%   orginalInput is the location of the input folder of single-GPU model.
%   multiInput is the location of the input folder of multi-GPU model.
%   numGPU is the number of GPUs.
%   DomainDecomposite(orginalInput,multiInput,numGPU,field_file_names) only
%   convert the files specified by the field_file_names
%   DomainDecomposite(orginalInput,multiInput,numGPU,sectionPropotion)
%   divide the files into subsection inputs specified by a section
%   propotion. The propotion should be a vector with numGPU elements and  
%   the sum of the vector should be 1. 
%
%   Example
%   -------
%   % convert manning.dat & h.dat
%   DomainDecomposite(orginalInput,multiInput,numGPU,'manning','h')
%   DomainDecomposite(orginalInput,multiInput,4,[0.1 0.4 0.3 0.2],'manning','h')
%   Updated by Xiaodong Ming on 8 Mar 2018.
NumOverlappedRows = 1; %currently only support 1 interchange row

%% number of GPU and Propotion
if numel(numGPU)==1&&numGPU>1&&round(numGPU)==numGPU
    sectionPropotion = repmat(1/numGPU,[1 numGPU]);
    sectionPropotion(end) = 1-sum(sectionPropotion(1:end-1));
else % input argument is a vector, gives the propotion of division
    sectionPropotion = numGPU;
    numGPU = numel(sectionPropotion);
    disp(['The input files are dividing into ' num2str(numGPU) ' sections base on the propotion:']);
end
CumSectionPropotion = cumsum(sectionPropotion); % cumulative propotion
%% files to be divided
fieldFileList = dir([orginalInput '/input/field/*.dat']);
fieldFileNames = struct2cell(fieldFileList); 
fieldFileNames = fieldFileNames(1,:);
globalFileNames = {'h',...
    'eta',...
    'hU',...
    'precipitation',...
    'z',...
    'manning',...
    'sewer_sink',...
    'cumulative_depth',...
    'hydraulic_conductivity',...
    'capillary_head',...
    'water_content_diff',...
    'gauges_pos',...
    'precipitation_source_all',...
    'precipitation_mask'...
    };
globalFileWrittingFlag = zeros(size(globalFileNames));
if isempty(varargin)
    globalFileWrittingFlag = globalFileWrittingFlag+1;
    list = fieldFileList;
else
    inputFileNames = varargin;
    lia = ismember(globalFileNames,inputFileNames);
    globalFileWrittingFlag(lia)=1;
    list = struct('name',[]);
    for i=1:length(inputFileNames)
        filename = [inputFileNames{i} '.dat'];
        list(i).name = filename;
        if ~ismember(filename,fieldFileNames)
            error(['Can not find ' filename ' in ' orginalInput '/input/field/'])
        end
    end
end
%% get global files to be divided
if ~isempty(varargin)
    list = struct('name',[]);
    for i=1:length(varargin)
        list(i).name = [varargin{i} '.dat'];
    end
end
%% divide dem matrix
[Z,R] =  ArcgridreadM([orginalInput '/input/mesh/DEM.txt']);

ValidCell_Ind = ~isnan(Z); % 1: valid cell 0: invalid cell
NumValidCell = sum(ValidCell_Ind(:)); % number of valid cells
NumOfValidCellRows = sum(ValidCell_Ind,2); clear ValidCell_Ind
CumNumOfValidRows = cumsum(NumOfValidCellRows,'reverse'); % cumulative number of valid cell in each row
CumNumValidCell = CumSectionPropotion*NumValidCell; %cumulative number of valid cell in each section
ind = repmat(CumNumOfValidRows,[1,numGPU])>repmat(CumNumValidCell,size(CumNumOfValidRows));
SecRow = sum(ind); SecRow = [size(Z,1) SecRow];
SectionRowNumber = cell(numGPU,1);
for i = 1:numGPU
    if i == numGPU
        first = 1;
    else
        first = 1 - NumOverlappedRows;
    end
    if i == 1
        last = SecRow(i);
    else
        last = SecRow(i) + NumOverlappedRows;
    end
    SectionRowNumber{i} = SecRow(i+1)+first:last;
end
Shared_ValidID_BoundVec = cell(numGPU,1);
GlobalToLocal_Valid_ID_Convert = zeros(numGPU,2);
BoundaryRow = cell(4,numGPU);
NODATA_value = -9999;
%% divide DEM.txt & gauges_pos.dat
% Local_Valid_ID = cell(NumSec,1);
% Local_Bound_ID = cell(NumSec,1);
gauges_xy_global = dlmread([orginalInput '/input/field/gauges_pos.dat']);
gauges_x = gauges_xy_global(:,1);
gauges_y = gauges_xy_global(:,2);
for i = 1:numGPU
    SectionRowNumberOfThis = SectionRowNumber{i};
    GlobalToLocal_Valid_ID_Convert(i,1) = CumNumOfValidRows(SectionRowNumberOfThis(end))...
        - NumOfValidCellRows(SectionRowNumberOfThis(end));
    GlobalToLocal_Valid_ID_Convert(i,2) = CumNumOfValidRows(SectionRowNumberOfThis(1));
    if exist([multiInput '/' num2str(i-1) '/output'], 'dir')~=7
        mkdir([multiInput '/' num2str(i-1) '/output']) % create output folder
    end
    folderName = [num2str(i-1) '/input/mesh/'];
    DEMName = [multiInput '/' folderName 'DEM.txt'];
    if exist([multiInput '/',folderName], 'dir')~=7
        mkdir([multiInput '/'],folderName) % create mesh folder
    end
    if exist([multiInput  '/' num2str(i-1) '/input/field'], 'dir')~=7
        mkdir([multiInput  '/' num2str(i-1) '/input/field']) % create field folder
    end
    
    Z_sub = Z(SectionRowNumber{i},:);
    Z_sub(isnan(Z_sub)) = NODATA_value;
    R_sub = R;
    R_sub(3,2) =  R(3,2)-abs(R(1,2))*(min(SectionRowNumber{i})-1);
    if min(globalFileWrittingFlag)>0 %whether to write sub DEM
        Arcgridwrite(DEMName,Z_sub,R_sub)
        if i== numGPU
            fprintf('%s\n ', 'DEM.txt divided')
        end
        if exist([orginalInput '/input/times_setup.dat'],'file')
        copyfile([orginalInput '/input/times_setup.dat'],multiInput);
        end
    end
    xmin = R_sub(3,1)+0.5*R_sub(2,1); xmax = xmin+size(Z_sub,2)*R_sub(2,1);
    ymax = R_sub(3,2)+0.5*R_sub(2,2); ymin = ymax+size(Z_sub,1)*R_sub(1,2);
    ind = gauges_x>=xmin&gauges_x<=xmax&gauges_y>=ymin&gauges_y<=ymax;
    if min(globalFileWrittingFlag)>0||globalFileWrittingFlag(12)
        gauges_xy_local = gauges_xy_global(ind,:);
        gauge_index_local = find(ind==1);
        
        writePath = [multiInput '/' num2str(i-1) '/input/field/'];
        writename = 'gauges_pos.dat';
        fileID = fopen([writePath writename],'w');
        if isempty(gauges_xy_local)
            gauges_xy_local = [0 0];
        end
        fprintf(fileID,'%.3f %.3f\n',gauges_xy_local(1:end-1,:)');
        fprintf(fileID,'%.3f %.3f',gauges_xy_local(end,:));
        fclose(fileID);
        if isempty(gauge_index_local)
            gauge_index_local=-1;
        end
        dlmwrite([writePath 'gauges_index.dat'],gauge_index_local,'delimiter',' ');
        if i==numGPU
            fprintf('%s\n ','gauges_pos.dat')
            fprintf('%s\n ','gauges_index.dat')
        end
    end
    if i==1
        BoundOption = 'top';
    elseif i==numGPU
        BoundOption = 'bottom';
    else
        BoundOption = 'all';
    end
    [Valid_ID,Bound_ID] = gen_VB_ID(Z_sub,R_sub,'ExportSharedBound',BoundOption);
    %     Local_Valid_ID{NumSec-i+1} = Valid_ID;
    %     Local_Bound_ID{NumSec-i+1} = Bound_ID;
    SharedValidID = Valid_ID(Bound_ID==-1);
    Shared_ValidID_BoundVec{i} =  [SharedValidID repmat([4 0 0],length(SharedValidID),1)];
    top_row = Valid_ID(1,:);
    top_row(isnan(top_row))=[];
    second_top_row = Valid_ID(2,:);
    second_top_row(isnan(second_top_row))=[];
    bottom_row = Valid_ID(end,:);
    bottom_row(isnan(bottom_row))=[];
    second_bottom_row = Valid_ID(end-1,:);
    second_bottom_row(isnan(second_bottom_row))=[];
    
    BoundaryRow{1,i} = top_row;
    BoundaryRow{2,i} = bottom_row;
    BoundaryRow{3,i} = second_top_row;
    BoundaryRow{4,i} = second_bottom_row;
    
end
%% generate shared_boundary 'halo.dat'
if min(globalFileWrittingFlag)>0 %whether to write halo.dat
    filename = [multiInput  '/' 'halo.dat'];
    fileID = fopen(filename,'w');
    fprintf(fileID,'No. of Domains\n');
    fprintf(fileID,'%d \n',numGPU);
    for i = 1:numGPU
        %write section number
        fprintf(fileID,'#%d\n',i-1);
        %write bottom shared_boundary cell ID
        if i==1
            fprintf(fileID,'%d\n',[]); fprintf(fileID,'%d\n',[]); % no bottom, write two empty row
        else
            fprintf(fileID,'%d ',BoundaryRow{2,i});fprintf(fileID,'%d\n',[]);%bottom to receive
            fprintf(fileID,'%d ',BoundaryRow{4,i});fprintf(fileID,'%d\n',[]); %bottom to send
        end
        %write top shared_boundary cell ID
        if i==numGPU
            fprintf(fileID,'%d\n',[]); % no top, write two empty row
        else
            fprintf(fileID,'%d ',BoundaryRow{1,i});fprintf(fileID,'%d\n',[]); %top to receive
            fprintf(fileID,'%d ',BoundaryRow{3,i});fprintf(fileID,'%d\n',[]); %top to send
        end
    end
    fclose(fileID);
    fprintf('%s\n ','halo.dat')
end

% if min(globalFileWrittingFlag)>0
%     list = dir([orginalInput '/input/field/*.dat']);
% else
%     list = struct('name',[]);
%     ind = find(globalFileWrittingFlag==1);
%     for i=1:length(ind)
%         filename = [globalFileNames{ind(i)} '.dat'];
%         list(i).name = filename;
%     end
% end

for n = 1:length(list)
    if length(list(n).name)>3
        if strcmp(list(n).name(end-2:end),'dat')
            filename = list(n).name;
            fileID = fopen([orginalInput '/input/field/' filename]);
            FirstChar = fread(fileID,[1,1],'*char');
            fclose(fileID);
            if FirstChar=='$'
                % divide global file to local files
                %delimiter = '';
                if ~strcmp(filename,'precipitation_mask.dat')
                    has_boundary = 1;
                else
                    has_boundary = 0;
                end
                if filename(2) == 'U'
                    value_dim = 2;
                    %formatSpec1 = '%-12d %.8f %.8f\n';
                else
                    value_dim = 1;
                    %formatSpec1 = '%-12d %.8f\n';
                end
                %formatSpec2 = '%-12d %d %d %d\n';
                [GlobalID_InitialValue,GlobalID_BoundVectors] = ReadField([orginalInput '/input/field/' filename], value_dim, has_boundary);
                for i_sec = 1:numGPU
                    ID_Convert = GlobalToLocal_Valid_ID_Convert(i_sec,1);
                    if i_sec< numGPU
                        ind1 = GlobalID_InitialValue(:,1)>=GlobalToLocal_Valid_ID_Convert(i_sec,1)&GlobalID_InitialValue(:,1)<GlobalToLocal_Valid_ID_Convert(i_sec,2);
                    else
                        ind1 = GlobalID_InitialValue(:,1)>=GlobalToLocal_Valid_ID_Convert(i_sec,1);
                    end
                    LocalID_InitialValue = [GlobalID_InitialValue(ind1,1)-ID_Convert,GlobalID_InitialValue(ind1,2:end)];
                    if ~strcmp(filename,'precipitation_mask.dat')
                        if i_sec< numGPU
                            ind2 = GlobalID_BoundVectors(:,1)>=GlobalToLocal_Valid_ID_Convert(i_sec,1)...
                                &GlobalID_BoundVectors(:,1)< GlobalToLocal_Valid_ID_Convert(i_sec,2);
                        else
                            ind2 = GlobalID_BoundVectors(:,1)>=GlobalToLocal_Valid_ID_Convert(i_sec,1);
                        end
                        LocalID_BoundVectors = [GlobalID_BoundVectors(ind2,1)-ID_Convert,GlobalID_BoundVectors(ind2,2:end)];
                        LocalID_BoundVectors_AddShare = [LocalID_BoundVectors;Shared_ValidID_BoundVec{i_sec}];
                        LocalID_BoundVectors_AddShare = sortrows(LocalID_BoundVectors_AddShare,1);
                    else
                        LocalID_BoundVectors_AddShare = [];
                    end
                    % write file
                    writePath = [multiInput  '/' num2str(i_sec-1) '/input/field/'] ;
                    WriteField([writePath filename], LocalID_InitialValue, LocalID_BoundVectors_AddShare);
                    if i_sec==numGPU
                        fprintf('%s\n ',filename)
                    end
                end
            else
                %copy to local files
                for i_sec = 1:numGPU
                    writePath = [multiInput  '/' num2str(i_sec-1) '/input/field/'];
                    if ~strcmp(filename,'gauges_pos.dat')
                        copyfile([orginalInput '/input/field/' filename],writePath);
                        if i_sec==numGPU
                            fprintf('%s\n ',filename)
                        end
                    end
                end
            end
        end
    end
end
disp('Mission accomplished!');
end
%%
function [Valid_Cell_ID,Bound_Cell_ID] = gen_VB_ID(Z,R,varargin)
% generate the ID of Valid cells and Boundary cells with different bound type
% Valid_Cell_ID: 0,1,2,...,number_of_valid_cells-1
% Boundary_Cell_ID: outline_bound=1; shared_bound=-1; IO_bound:2,3,...,number_of_IO_bound-1
% Parameters: Z(DEM),R(Reference Matrix), IO_BoundFrame(Bound Frame Coordinates)
% Name-VAlue Pair Arguments
%       'IO_BoundFrame' -- [](default)|n*4 numeric matrix, each row represents one bound frame
%       'ExportSharedBound' -- 'none'(default),'all','top','bottom'| export certain row of the shared bound
% Created by Xiaodong Ming on 28/2/2016.
% Updated on 30/11/2016

NumOfVarargin = length(varargin);
DefaultParamNames = {'IO_BoundFrame','ExportSharedBound'};
DefaultParamValues = {[],'none'};
if NumOfVarargin == 1 % three input parameters, the 3rd is IO_BoundFrame
    DefaultParamValues{1} = varargin{1};
elseif NumOfVarargin == 2 || NumOfVarargin == 4
    for i = 1:NumOfVarargin/2
        if strcmp(varargin(i*2-1),DefaultParamNames{1})
            DefaultParamValues{1} = varargin(i*2);
        elseif strcmp(varargin(i*2-1),DefaultParamNames{2})
            DefaultParamValues{2} = varargin(i*2);
        else
            error('Wrong argument')
        end
    end
elseif NumOfVarargin == 3
    if ismatrix(varargin{1})
        DefaultParamValues{1} = varargin{1};
        if strcmp(varargin(2),DefaultParamNames{2})
            DefaultParamValues{2} = varargin(3);
        else
            error('Wrong arguments: argument 4')
        end
    else
        error('Wrong number of arguments: argument 3 should be matrix')
    end
elseif NumOfVarargin == 0
else
    error('Wrong number of arguments')
end
IO_BoundFrame = DefaultParamValues{1};

if isempty(Z)
    error('DEM matrix is empty!')
end
nrows = size(Z,1);
ncols = size(Z,2);
nan_value = -9999;
% ****1.assignment of the valid cell ID
% Valid_Cell_ID = Z*0-1; valid_cell_cnt = 0; % number of valid cell
% for i = nrows:-1:1
%     for j = 1:ncols;
%         if Z(i,j)>nan_value
%             Valid_Cell_ID(i,j) = valid_cell_cnt;
%             valid_cell_cnt = valid_cell_cnt+1;
%         end
%     end
% end
indicate_valid_cell = Z > nan_value;
sum_in_eachrow = sum(indicate_valid_cell,2);
cumsum_in_eachrow = cumsum(sum_in_eachrow,'reverse') - sum_in_eachrow;
auxilliary_matrix = [cumsum_in_eachrow indicate_valid_cell];
new_index = cumsum(auxilliary_matrix,2);
new_index(:,1) = [];
new_index(indicate_valid_cell == 0) = NaN;
Valid_Cell_ID = new_index - 1;
% valid_cell_cnt = cumsum_in_eachrow(1) + sum_in_eachrow(1);

% ****2. assignment of the boundary cell: outline, IO, shared
% 2.1 assignment of the outline boundary cell
% Bound_Cell_ID = Z*0; boundary_cell_cnt = 0; % number of boundary cell
% for i = nrows:-1:1
%     for j = 1:ncols;
%         if Valid_Cell_ID(i,j)>=0
%             if i == 1 || i == nrows || j == 1 || j == ncols % cells on the box margin
%                 Bound_Cell_ID(i,j) = 1; % outline bound cell is assigned as 1
%                 boundary_cell_cnt = boundary_cell_cnt + 1;
%             else
%                 L_cell = Valid_Cell_ID(i,j-1);
%                 R_cell = Valid_Cell_ID(i,j+1);
%                 U_cell = Valid_Cell_ID(i+1,j);
%                 D_cell = Valid_Cell_ID(i-1,j);
%                 if isnan(L_cell+R_cell+U_cell+D_cell)
%                     Bound_Cell_ID(i,j) = 1;
%                     boundary_cell_cnt = boundary_cell_cnt + 1;
%                 end
%             end
%         end
%     end
% end

shift_up = [indicate_valid_cell(2:end,:); zeros(1, ncols)];
shift_down = [zeros(1, ncols); indicate_valid_cell(1:end-1,:)];
shift_right = [zeros(nrows, 1) indicate_valid_cell(:,1:end-1)];
shift_left = [indicate_valid_cell(:,2:end) zeros(nrows, 1)];
auxilliary_matrix = shift_up + shift_down + shift_right + shift_left;
Bound_Cell_ID = Z*0;
Bound_Cell_ID(auxilliary_matrix < 4) = 1;
Bound_Cell_ID(indicate_valid_cell == 0) = NaN;

% 2.2 assignment of the I/O boundary cell
if ~isempty(IO_BoundFrame)
    f_N = 1:size(IO_BoundFrame,1); % number of IO bound
    b_N = f_N+1; % values assigned to the IO bounds
    for i = 1:length(f_N)
        % transfer map coordinate to pixel coordinate of the Z matrix (col and row)
        row_col = map2pix(R,[IO_BoundFrame(i,1:2);IO_BoundFrame(i,3:4)]);
        row_col = ceil(row_col);
        row_col(row_col<1)=1;
        row_col(row_col(:,1)>nrows,1) = nrows;
        row_col(row_col(:,2)>ncols,2) = ncols;
        
        frame_row = min(row_col(:,1)):max(row_col(:,1));
        frame_col = min(row_col(:,2)):max(row_col(:,2));
        
        [rectX,rectY] = meshgrid(frame_row,frame_col); % matrix of row and col number
        % convert row and column to index of matrix
        ind = sub2ind(size(Bound_Cell_ID),rectX(:),rectY(:));
        BM = Bound_Cell_ID(ind);
        Bound_Cell_ID(ind(BM==1)) = b_N(i);
    end
end
% 2.3 assignment of the shared boundary cell
top_row = Valid_Cell_ID(1,:);
extend_row = [nan top_row nan];
indcate_row = top_row+extend_row(1:end-2)+extend_row(3:end);
top_row_ind = ~isnan(indcate_row);

bottom_row = Valid_Cell_ID(end,:);
extend_row = [nan bottom_row nan];
indcate_row = bottom_row+extend_row(1:end-2)+extend_row(3:end);
bottom_row_ind = ~isnan(indcate_row);

if strcmp(DefaultParamValues{2},'all')
    Bound_Cell_ID(1,top_row_ind) = -1;
    Bound_Cell_ID(end,bottom_row_ind) = -1;
elseif strcmp(DefaultParamValues{2},'top')
    Bound_Cell_ID(1,top_row_ind) = -1;
elseif strcmp(DefaultParamValues{2},'bottom')
    Bound_Cell_ID(end,bottom_row_ind) = -1;
elseif strcmp(DefaultParamValues{2},'none')
else
    warning('no shared boundary ID because of the wrong ExportSharedBound argument')
end
end

%%
function [GlobalID_InitialValue,GlobalID_BoundVectors] = ReadField(filename,value_dim,has_boundary)
fid = fopen(filename,'rt');
s = fread(fid,'*char');
out = textscan(s,'%f','CommentStyle',{'$', newline},'MultipleDelimsAsOne',1);
skip = value_dim + 1;
C = out{1,1};
NumElements = C(1);
GlobalID_InitialValue = reshape(C(2:NumElements*skip+1),skip,NumElements)';
GlobalID_BoundVectors = [];
if has_boundary == 1
    NumBound = C(NumElements*skip+2);
    GlobalID_BoundVectors = reshape(C(NumElements*skip+3:end),4,NumBound)';
end
fclose(fid);
end
%%
function [Z,R] = ArcgridreadM(filename)
%ArcgridreadM Simplified version of arcgridread. asc file must have 6 head lines
% Created by Xiaodong Ming on 2017-11-27
% See also arcgridread
fileID = fopen(filename);
C = textscan(fileID,'%s %f',6);
ncols = C{2}(1);
nrows = C{2}(2);
xllcorner = C{2}(3);
yllcorner = C{2}(4);
cellsize = C{2}(5);
nanValue = C{2}(6);
formatSpec = [repmat('%f',[1,ncols]) '%[^\n\r]'];
Z = textscan(fileID,formatSpec,nrows);
fclose(fileID);
Z = [Z{1:end-1}];
Z(Z==nanValue) = nan;
x11 = xllcorner+0.5*cellsize;
y11 = yllcorner+(nrows-0.5)*cellsize;
R = makerefmat(x11,y11,cellsize,-cellsize);
end

%%
function [] = WriteField(filename, GlobalID_InitialValue, GlobalID_BoundVectors)
[rows, cols] = size(GlobalID_InitialValue);
str_head = sprintf('$Element Number$\n%d\n$Element_id  Value$\n',rows);
if cols == 2
    formatSpec = '%-12d %.8f\n';
else
    formatSpec = '%-12d %.8f %.8f\n';
end
str_value = sprintf(formatSpec,GlobalID_InitialValue');
if ~isempty(GlobalID_BoundVectors)
    [rows, ~] = size(GlobalID_BoundVectors);
    str_bound_head = sprintf('$Boundary Numbers$\n%d\n$Element_id Boundary_type$\n',rows);
    str_bound_value = sprintf('%-12d %d %d %d\n',GlobalID_BoundVectors');
else
    str_bound_head = [];
    str_bound_value = [];
end
str = [str_head str_value str_bound_head str_bound_value];
fid = fopen(filename,'w');
fwrite(fid,str,'char');
fclose(fid);
end

function Arcgridwrite(filename,Z,R)
%ARCGRIDWRITE Write gridded data set to Arc ASCII Grid Format
%
%   ARCGRIDWRITE(filename,Z,R) write a matlab format grid with its
%   reference into an Arc ASCII file Grid format.
%   Z is a 2D array containing the data values.  R is a
%   referencing matrix (see MAKEREFMAT).  -9999 is assigned to elements
%   of Z corresponding to null data values in the grid file.
%
%   See also ARCGRIDREAD.

% Created by Ming 2017-06-09.

% x11 and y11 specify the map location of the center of the first (1,1)
% pixel in the image or the first element of the data grid
x11 = R(3,1);
y11 = R(3,2);
gridSize = abs(R(2));
nrows = size(Z,1);
ncols = size(Z,2);

% XLLCORNER and YLLCORNER are the coordinates of the
% lower left corner of the lower left cell of Z.
xllcorner = x11 + gridSize/2;
yllcorner = y11 - gridSize*nrows - gridSize/2;
NODATA_value = -9999;
Z(isnan(Z)) = NODATA_value;

fileID = fopen(filename,'w');
fprintf(fileID,'ncols    %d\n', ncols);
fprintf(fileID,'nrows    %d\n', nrows);
fprintf(fileID,'xllcorner    %f\n', xllcorner);
fprintf(fileID,'yllcorner    %f\n', yllcorner);
fprintf(fileID,'cellsize    %f\n', gridSize);
fprintf(fileID,'NODATA_value    -9999\n');
dlmwrite(filename,Z,'-append','delimiter','\t')
fclose(fileID);
end
function [Z_value,R] = Field2Raster(CaseFolder,fieldFile)
%[Z_value,R] = Field2Raster(CaseFolder,fieldFile) convert field file to
%   raster file. CaseFolder is the location of the input folder. fieldFile 
%   is the name of the field file to be converted. If field file is hU, the
%   return Z_value is a cell consisted of two matrix. Othervise, the 
%   Z_value is a matrix and R is its reference mat.

% Created by Xiaodong Ming on 2017-8-5
% See also InputSetup, FieldSetup
if length(fieldFile)<3
    fieldFile = [fieldFile '.dat'];
elseif ~strcmp(fieldFile(end-3:end),'.dat')
    fieldFile = [fieldFile '.dat'];
end
fieldFileName = [CaseFolder '\input\field\' fieldFile];
demFileName = [CaseFolder '\input\mesh\DEM.txt'];
if exist(demFileName,'file')
    [Z,R] = arcgridread(demFileName);
else
    error(['Cannot find DEM file:' demFileName])
end
Z_value = nan(size(Z));
[Valid_ID,~] = gen_VB_ID(Z,R);%
ind1 = find(Valid_ID>=0);
[a,b] = ind2sub(size(Z),ind1);
ab = [a,b];
ab = sortrows(ab,[-1 2]);
ind2 =  sub2ind(size(Z),ab(:,1),ab(:,2));
NumOfCells = dlmread(fieldFileName,'',[1,0,1,0]);
if fieldFileName(end-4)=='U'
    colNum = 2;
    validID_FieldValue = dlmread(fieldFileName,'',[3,1,NumOfCells+2,colNum]);
    Z_value_x = Z_value;
    Z_value_y = Z_value;
    Z_value_x(ind2) = validID_FieldValue(:,1);
    Z_value_y(ind2) = validID_FieldValue(:,2);
    Z_value = {Z_value_x,Z_value_y};
else
    colNum = 1;
    validID_FieldValue = dlmread(fieldFileName,'',[3,1,NumOfCells+2,colNum]);
    Z_value(ind2) = validID_FieldValue;
end
% figure;mapshow(Z_value,R,'DisplayType','surface')
end

function [Valid_Cell_ID,Bound_Cell_ID] = gen_VB_ID(Z,R,varargin)
% gen_VB_ID generate the ID of Valid cells and Boundary cells with
% different bound type
% [Valid_Cell_ID,Bound_Cell_ID] = gen_VB_ID(Z,R,varargin)
% Valid_Cell_ID: a matrix with the same size of Z. Its values are
% 0,1,2,...,number_of_valid_cells-1
% Boundary_Cell_ID: a matrix with the same size of Z. Its values are
% outline_bound=1; shared_bound=-1; IO_bound:2,3,...,number_of_IO_bound-1.
% Parameters:
% Z(DEM), R(Reference Matrix)
% IO_BoundFrame(Points-pair Coordinates)[x11,y11,x12,y12; x21,y21,x22,y22;...]
% Name-Value Pair Arguments
%  'IO_BoundFrame' -- [](default)|n*4 numeric matrix, each row represents one bound frame
%  'ExportSharedBound' -- 'none'(default),'all','top','bottom'| export certain row of the shared bound
% Created by Xiaodong Ming on 28/2/2016.
% Updated on 30/11/2016
% Updated on 28/06/2017

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
Valid_Cell_ID = Z*0-1; valid_cell_cnt = 0; % number of valid cell
for i = nrows:-1:1
    for j = 1:ncols;
        if Z(i,j)>nan_value
            Valid_Cell_ID(i,j) = valid_cell_cnt;
            valid_cell_cnt = valid_cell_cnt+1;
        end
    end
end
% ****2. assignment of the boundary cell: outline, IO, shared
% 2.1 assignment of the outline boundary cell
Bound_Cell_ID = Z*0; boundary_cell_cnt = 0; % number of boundary cell
for i = nrows:-1:1
    for j = 1:ncols;
        if Valid_Cell_ID(i,j)>=0
            if i == 1 || i == nrows || j == 1 || j == ncols % cells on the box margin
                Bound_Cell_ID(i,j) = 1; % outline bound cell is assigned as 1
                boundary_cell_cnt = boundary_cell_cnt + 1;
            else
                L_cell = Valid_Cell_ID(i,j-1);
                R_cell = Valid_Cell_ID(i,j+1);
                U_cell = Valid_Cell_ID(i+1,j);
                D_cell = Valid_Cell_ID(i-1,j);
                if isnan(L_cell+R_cell+U_cell+D_cell)
                    Bound_Cell_ID(i,j) = 1;
                    boundary_cell_cnt = boundary_cell_cnt + 1;
                end
            end
        end
    end
end

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
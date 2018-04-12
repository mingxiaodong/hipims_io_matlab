function FieldSetup(CaseFolderAddress, Z, R, varargin)
%FIELDSETUP     set up field files for hydrodynamic model
%   FieldSetUp(CaseFolderAddress, Z, R) set up field files of a case, all other parameters are default value. 
%   CaseFolderAddress is the location of case folder. Z is the DEM value matrix. R is the DEM reference.
%   FieldSetUp(CaseFolderAddress, Z, R, Name, Value) 
%   Name-VAlue Pair Arguments
%       'h_Eta' -- identification of writing  whether water depth or elevation files
%           'h'(default)|'eta'
%       'IO_BoundFrame' -- coordinates of the input or output bound frame
%           [](default)|n*4 numeric matrix, each row represents one bound frame
%       'BoundCode' -- boundary code
%           [2 0 0; 2 1 0](default)|2*3n numeric matrix, n is the number of bounds
%       'h_BC_Source' -- h boundary condition source data
%           {[0 0]}(default)|cell matrix consists of certain number of 2-column numeric matrixs
%       'hU_BC_Source' -- hU boundary condition source data
%           {[0 0 0]}(default)|cell matrix consists of certain number of 
%           3-column numeric matrix for discharge per unit width(m2/s)
%           2-column numeric matrix for total discharge in one bound (m3/s)
%       'initial_hE_hU_pre' -- initial value of h/eta, hU, and precipitation
%           {0, {0 0}, 0}(default)| 3-element cell: scalar/matrix with the
%           same size of Z
%       'hydro_params_Value' -- hydrological parameters value: 1-manning,
%               2-sewer_sink, 3-cumul_depth, 4-hydraulic_conductivity,
%               5-capillary_head, 6-water_content_diff
%           {0.035,0,0,0,0,0}(default)| 4-element cell: scalar or 
%               matrix with the same size of Z
%       'RainMask' -- rainfall mask: the No. of rainfall source
%           zeros(size(Z))(default)|integral matrix with the same size of Z
%       'RainSource' -- rainfall source values
%           {[0 0]}(default)| cell consists of two-column numeric matrices
%       'GaugeCoor' -- coordinates of gauge points
%           [](default)| Two-column(X Y) numerical matrix
%       'WriteAllFiles' -- logical value to indicate whether to write all files
%           false(default)| true: write all 14 files
%   See also gen_VB_ID, WriteInitialValue.
%   Created by Xiaodong Ming on 30 Sep 2016.
%   Updated by Xiaodong Ming on 08 May 2017.
%   Updated by Xiaodong Ming on 28 Jun 2017.

%% location and write path
% SourceValueCell, MaskValueCell, InitialValueCell, h_eta_string
if nargin < 3
    error('FieldSetUp::WrongNumberOfInputs')
end
writePath = [CaseFolderAddress '/input/field/'];
Delimiter = ' ';
%% define default input param names and values
% Name-Value pairs: 11
AllInputParamNames = {'h_Eta',      'IO_BoundFrame',    'BoundCode', 'h_BC_Source',    'hU_BC_Source',...
         'initial_hE_hU_pre', ...
         'hydro_params_Value',...
         'RainMask',  'RainSource',       'GaugeCoor', 'WriteAllFiles',...
         'initial_hE', 'initial_hU',...
         'manning','sewer_sink','cumul_depth','hydraulic_conductivity','capillary_head','water_content_diff',};
OriginParamValues  = {'h'   ,       []            , [2 0 0; 2 1 0], {[0 0; 60 0]}, {[0 0 0; 60 0 0]},...
               {0, {0 0}, 0},          {0},        {0 0},...
               {0.035,0,0,0,0,0} ,0.035, 0,0,0,0,0,...
               zeros(size(Z)), {[0 0; 60 0]},            []    ,    false      };
n = length(varargin)/2; % number of Name-Value pairs from the input parameter
if n~=fix(n)
    error('Wrong number of varargin parameters')
end
%% Writting Flags
% field files 15
WrittingFlag = zeros(15,1); % flags about whether to write the 15 field files
%  1. h_BC_num.dat % could be more than one files
%  2. hU_BC_num.dat % could be more than one files
%  3. h.dat
%  4. hU.dat
%  5. precipitation.dat
%  6. manning.dat
%  7. sewer_sink.dat
%  8. cumulative_depth.dat
%  9. hydraulic_conductivity.dat
% 10. capillary_head.dat
% 11. water_content_diff.dat
% 12. precipitation_mask.dat
% 13. precipitation_source_num.dat % could be more than one file
% 14. gauges_pos.dat
% 15. z.dat
%% Parameters passing
if n > 0
    InputParamNames = varargin(1:2:n*2-1); % all the input parameter names
    InputParamValues = varargin(2:2:n*2);  % the values of corresponding input parameters
    for i = 1:n
        ParamName = InputParamNames{i};
        switch ParamName
            case AllInputParamNames{1} % h_Eta string: h or eta
                OriginParamValues{1} = InputParamValues{i};
%                 WrittingFlag(3) = 1;
            case AllInputParamNames{2} % IO_BoundFrame
                OriginParamValues{2} = InputParamValues{i};
%                WrittingFlag(1:2) = 1; % write bound source file
            case AllInputParamNames{3} % BoundCode
                OriginParamValues{3} = InputParamValues{i};
            case AllInputParamNames{4} % h_BC_Source
                OriginParamValues{4} = InputParamValues{i};
            case AllInputParamNames{5} % hU_BC_Source
                OriginParamValues{5} = InputParamValues{i};
            case AllInputParamNames{6} % initial_hE_hU_pre
                OriginParamValues{6} = InputParamValues{i};
                WrittingFlag([3,4,5]) = 1; %3. h.dat  %4. hU.dat %5. precipitation.dat
            case AllInputParamNames{7} % hydro_params_Value
                OriginParamValues{7} = InputParamValues{i};
                WrittingFlag([6,7,8,9,10,11]) = 1; % six files about hydro_params 
            case AllInputParamNames{8} % RainMask
                OriginParamValues{8} = InputParamValues{i};
                WrittingFlag(12) = 1;
            case AllInputParamNames{9}% RainSource
                OriginParamValues{9} = InputParamValues{i};
                WrittingFlag(13) = 1;
            case AllInputParamNames{10} % GaugeCoor
                OriginParamValues{10} = InputParamValues{i};
                WrittingFlag(14) = 1;
            case AllInputParamNames{11} % WriteAllFiles
                OriginParamValues{11} = InputParamValues{i};
            otherwise
                error(['Wrong parameters name ' ParamName])
        end
    end
end
%% Initialize internal variables
if OriginParamValues{11};
    WrittingFlag = ones(15,1);
end
%write water depth(h) or water level(eta) file
h_Eta = OriginParamValues{1};
% coordinates of IO Bound Frame
IO_BoundFrame = OriginParamValues{2}; % coordinates of IO Bound Frame
if ~isempty(IO_BoundFrame)
    if size(IO_BoundFrame,2)~=4
        error('Wrong IO_bound frame coordinates')
    end
end
% Bound_Attributes_Code
Bound_Code = OriginParamValues{3}; 
% boundary source data of water depth/level
h_BC_source_cell  = OriginParamValues{4};
% boundary source data of water velocity
hU_BC_source_cell = OriginParamValues{5};
% initial h, vel, and precipitation value
initial_hE_hU_pre = OriginParamValues{6};

% hydro params value (1-manning, 2-sewer_sink, 3-cumulative_depth,
%   4-hydraulic_conductivity, 5-capillary_head, 6-water_content_diff;
hydro_params_Value = OriginParamValues{7}; % 6-element cell
hydro_params_Names = {'manning','sewer_sink','cumulative_depth','hydraulic_conductivity','capillary_head','water_content_diff'};

% rainfall mask and rainfall source data
RainMask = OriginParamValues{8};
if size(RainMask)==size(Z)|numel(RainMask)==1;
    RainMask = zeros(size(Z))+RainMask;
else
    error('Wrong Rain mask value')
end
RainSourceValue = OriginParamValues{9};
%% generate valid ID and boundary ID for each grid cell 
% invoke sub function "gen_VB_ID"
[Valid_ID,Bound_ID] = gen_VB_ID(Z,R,IO_BoundFrame);%
% disp('Valiad grid cell ID and Boundary grid cell ID have been generated')
%% ****** write files *********
%% * write boundary condition files: h/eta_BC.dat and hU_BC.dat 
%***output file No. 1,2: h_BC_num.dat, hU_BC_num.dat
% number of bounds in total
NumOfBound = size(IO_BoundFrame,1)+1;
if size(Bound_Code,2)==3&&NumOfBound>1
    Bound_Code = repmat(Bound_Code,[1,NumOfBound]);
    Warning('No Bound_Code for your boundary cell: set as open boundary')
end
% the first values of Bound_Attributes_Code for each bound
h_BC_Code_1  = Bound_Code(1,1:3:3*NumOfBound);
hU_BC_Code_1 = Bound_Code(2,1:3:3*NumOfBound);
h_BC_Num  = 1;
hU_BC_Num = 1;
for i = 1:NumOfBound % bound ID
    % judge whether it is a IO bound
    if h_BC_Code_1(i)  == 3
        h_BC_source_num  = Bound_Code(1,3+3*(i-1));
        h_BC_index_name  = [h_Eta '_BC_' num2str(h_BC_source_num) '.dat'];
        if length(h_BC_source_cell) < h_BC_Num
            error('lack of source matrix for h_BC')
        end
        h_BC_source_value = h_BC_source_cell{h_BC_Num}; h_BC_Num = h_BC_Num+1;
        if size(h_BC_source_value,2)~=2
            warning(['Wrong boundary condition value of hEta for Bound ' num2str(i)])
        end
        % output file No 1: h_BC_num.dat
        if WrittingFlag(1)==1
            dlmwrite([writePath h_BC_index_name],h_BC_source_value,...
                'precision',10,'delimiter',Delimiter)
            disp([h_BC_index_name ' for Bound' num2str(i)])
        end
    end
    if hU_BC_Code_1(i) == 3
        hU_BC_source_num = Bound_Code(2,3+3*(i-1));
        hU_BC_index_name = ['hU_BC_' num2str(hU_BC_source_num) '.dat'];
        if length(hU_BC_source_cell) < hU_BC_Num
            error('lack of source matrix for hU_BC')
        end
        TempValue = hU_BC_source_cell{hU_BC_Num}; 
        if size(TempValue,2)==2
            BoundWidth = abs(R(2))*sum(Bound_ID(:)==i);
            if BoundWidth > 0
                Q_per_m = TempValue(:,2)/BoundWidth; % m^2/s
                [I,J] = ind2sub(size(Bound_ID),find(Bound_ID(:)==i)); %range(cols):delta_x, range(a):delta_y
                theta = atan(range(J)/range(I)); %theta: the angle between bound line and the vertical line
                hu = Q_per_m*cos(theta); hv = Q_per_m*sin(theta); %m^2/s
            else
                hu = TempValue(:,1)*0; hv = TempValue(:,1)*0; % no inflow bound, width is zero
            end
            hU_BC_source_value = [TempValue(:,1) hu hv];
        elseif size(TempValue,2)==3
            hU_BC_source_value = hU_BC_source_cell{hU_BC_Num};
        else
            warning(['Wrong boundary condition value hU for Bound ' num2str(i)])
        end
        hU_BC_Num = hU_BC_Num+1;
        % output file No 2: hU_BC_num.dat
        if WrittingFlag(2)==1
            dlmwrite([writePath hU_BC_index_name],hU_BC_source_value,...
                'precision',10,'delimiter',Delimiter)
            disp([hU_BC_index_name ' for Bound' num2str(i)])
        end
    end
end
if WrittingFlag(1)==1
    disp(['FileNo-01: h_BC_num has been created, total number of h source file: ' num2str(h_BC_Num-1)])
end
if WrittingFlag(2)==1
    disp(['FileNo-02: hU_BC_num has been created, total number of hU source file: ' num2str(hU_BC_Num-1)])
end
%% * write initial depth, velocity, and precipitation files
    % h/eta,dat, hU.dat, and precipitation.dat  
%hE_hU_pre_value = {hE_value,hU_value,pre_value};
hE_hU_pre_name = {h_Eta,'hU','precipitation'};
%***output file No. 3,4,5: h.dat, hU.dat, precipitation.dat  
for i=1:length(hE_hU_pre_name)
    write_value = initial_hE_hU_pre{i};
    if iscell(write_value)
        if (size(write_value{1})==size(Z)|numel(write_value{1})==1)...
                & (size(write_value{2})==size(Z)|numel(write_value{2})==1)
            write_value = [zeros(size(Z))+write_value{1}, zeros(size(Z))+write_value{2}]; % initial hU discharge per unit width
        else
            error(['Wrong initial ' hE_hU_pre_name{i} ' data: does not match with the size of DEM matrix'])
        end
    else
        if size(write_value)==size(Z)|numel(write_value)==1
            write_value =  zeros(size(Z)) + write_value; %initial h water depth or eta water surface elevation
        else
            error(['Wrong initial ' hE_hU_pre_name{i} ' data: does not match with the size of DEM matrix'])
        end
    end
    if WrittingFlag(2+i)==1
        if strcmp(hE_hU_pre_name{i},'hU')
            if size(write_value)==size(Z)
                write_value = repmat(write_value,[1 2]);
            end
        end
        WriteInitialValue(writePath,hE_hU_pre_name{i},Valid_ID,Bound_ID,...
            write_value,Bound_Code)
        disp(['FileNo-' num2str(2+i,'%02u') ': ' hE_hU_pre_name{i} '.dat has been created'])
    end
end

%% * write hydro params files
    % A manning.dat, B sewer_sink.dat, C cumulative_depth.dat,
    % D hydraulic_conductivity.dat, E capillary_head.dat,
    % F water_content_diff.dat
%***output file No. 6~11
for i = 1:length(hydro_params_Names)
    param_name = hydro_params_Names{i};
    param_value = hydro_params_Value{i};
    if numel(param_value)~=1 & size(param_value)~=size(Z)
        error(['Wrong ' param_name ' value'])
    end
    grid_param_Value = zeros(size(Z))+param_value;
    if WrittingFlag(5+i)==1
        WriteInitialValue(writePath,param_name,Valid_ID,Bound_ID,grid_param_Value,Bound_Code)
        disp(['FileNo-' num2str(5+i,'%02u') ': ' param_name '.dat has been created'])
    end
end

%% write Precipitation Mask file
% precipitation_mask.data
VID_RMask = [Valid_ID(:) RainMask(:)]; % zeros(numel(Valid_ID),1) rainfall mask value is 0 
VID_RMask(isnan(Valid_ID(:))|Valid_ID(:)==-1,:)=[]; % delete nodata value
VID_RMask = sortrows(VID_RMask,1);
% write file
%***output file No.12: precipitation_mask.dat 
if WrittingFlag(12)==1
    filename = 'precipitation_mask.dat';
    fileID = fopen([writePath filename],'w');
    fprintf(fileID,'$Element Number\n'); fprintf(fileID,'%d\n',length(VID_RMask));
    fprintf(fileID,'$Element_id  Value\n'); fprintf(fileID,'%-12d %d\n',VID_RMask');
    fclose(fileID);
    disp(['FileNo-12: ' filename ' has been created'])
end
%% write Precipitation Source files: 
%***output file No.13: precipitation_source_num.dat
if WrittingFlag(13)==1
    for i = 1:length(RainSourceValue)
        Time_Rain = RainSourceValue{i};
        Time_Rain(isnan(Time_Rain(:,2)),:)=[];
        filename = ['precipitation_source_' num2str(i-1) '.dat'];
        dlmwrite([writePath filename],Time_Rain,'precision',10,'delimiter',Delimiter);
    end
    disp(['FileNo-13: ' num2str(i) ' precipitation_source files have been created'])
end
RainMaskNum = unique(RainMask(:)); RainMaskNum(isnan(RainMaskNum))=[];
if numel(RainMaskNum)~=length(RainSourceValue)
    warning('The number of rain mask value is not the same as the number of rainfall source files.')
end
%% *****gauges_pos.dat
%coordinate of gauge
%***output file No.14:gauges_pos.dat
GaugeCoor = OriginParamValues{10};
if ~isempty(GaugeCoor)
    if size(GaugeCoor)~=2
        error('Wrong dimension of GaugeCoor')
    end
end
filename = 'gauges_pos.dat';
if WrittingFlag(14)==1
    fileID = fopen([writePath filename],'w');
    if isempty(GaugeCoor)
        GaugeCoor = [0 0];
    end
    fprintf(fileID,'%.3f %.3f\n',GaugeCoor(1:end-1,:)');
    fprintf(fileID,'%.3f %.3f',GaugeCoor(end,:));
    fclose(fileID);
%     dlmwrite([writePath filename],GaugeCoor,'precision',10,'delimiter',Delimiter);
    disp(['FileNo-14: ' filename ' has been created'])
end
%% **** z.dat
if WrittingFlag(15)==1
   WriteInitialValue(writePath,'z',Valid_ID,Bound_ID,Z,Bound_Code)
   disp('FileNo-15: z.dat has been created')
end
%% *****summary information
% if min(WrittingFlag)==1
%     disp('All field files have been created')
%     figure;
%     fig_hd = mapshow(Z,R,'DisplayType','surface');
%     for NOB = 1:NumOfBound-1
%         rect_CoorX = IO_BoundFrame(NOB,[1 3]);
%         rect_CoorY = IO_BoundFrame(NOB,[2 4]);
%         rec_hd = rectangle('Position',[min(rect_CoorX),min(rect_CoorY),range(rect_CoorX),range(rect_CoorY)]);
%         rec_hd.LineWidth = 3; rec_hd.EdgeColor = [1 0 0];% 1 0 0: red
%     end
%     zdatam(fig_hd,Z-max(Z(:)))
% end
end

function WriteInitialValue(writePath,fileName,Cell_ID,Bound_ID,InitialValue,BoundType)
% default bound type : 1: normal bound; 2: input bound; 3: output bound
% Bound_ID: -1: shared bound; 1: outline bound; 2,3,4,...: IO bound
% Cell_ID: nan, 0 ~ number_of_valid_cell-1
% fileNames: 'h','eta','hU' are three special files, others including
% 'z','precipitation','manning'...
% BoundType_Vec is cell containing 2*3 matrix  [h;hU]*{bound1 2 3}
adjustValue = 0.001; % avoid all the value to be zero in 'h', 'eta', 'hU';
Outline_Bound_Value = 2; %[2 0 0]
Shared_Bound_Value = 4; %[4 0 0]
BoundType_Vec = cell(1,size(BoundType,2)/3);
for i = 1:length(BoundType_Vec)
    BoundType_Vec{i} = BoundType(:,(1:3)+(i-1)*3);
end

%% ****id_IV: valid cells' ID and their initial value
if strcmp(fileName,'h')||strcmp(fileName,'eta') % avoid the depth in the flow-input-bound cells to be nil
    InitialValue(Bound_ID>=2) = InitialValue(Bound_ID>=2)+adjustValue;
end 

ValidCellValue = InitialValue(:); %Valid cells' value
if strcmp(fileName,'hU') % file name is hU, hU consist of [hu hv]
    id_IV = [Cell_ID(:), ValidCellValue(1:end/2), ValidCellValue(end/2+1:end)];
    formatSpec1 = '%-12d %.8f %.8f\n';    % format for valid cells' ID and their initial values 
else  %if filename is z, h, precepitation, or manning   
    id_IV = [Cell_ID(:), ValidCellValue];
    formatSpec1 = '%-12d %.8f\n'; % format for valid cells' ID and their initial values
end
id_IV(isnan(id_IV(:,1))|id_IV(:,1)<0,:) = []; % remove the invalid cells and their value
id_IV(isnan(id_IV)) = 0; % replace the NaN initial value to zero
id_IV = sortrows(id_IV,1);

%% ****id_BTV: boundary cells' ID & their boundary type vec
formatSpec2 = '%-12d %d %d %d\n'; % format for bound cells' ID and their boundary vector
Cell_Bound_ID = [Cell_ID(:),Bound_ID(:)]; 
Cell_Bound_ID(isnan(Cell_Bound_ID(:,1)),:) = []; % remove the invalid cells
Cell_Bound_ID(Cell_Bound_ID(:,2)==0,:)=[]; % remove the non-boundary cells
Cell_Bound_ID = sortrows(Cell_Bound_ID,1); % only boundary cells left

bound_cnt = length(Cell_Bound_ID);
BondType_Vector_Print = zeros(bound_cnt,3); % default Value [0 0 0]
% default outline bound type vector [2 0 0]
BondType_Vector_Print(:,1) = Outline_Bound_Value;
% default shared bound type vector [4 0 0]
BondType_Vector_Print(Cell_Bound_ID(:,2)==-1,1) = Shared_Bound_Value;

if strcmp(fileName,'h')|| strcmp(fileName,'eta')||strcmp(fileName,'hU')
    if strcmp(fileName,'h')|| strcmp(fileName,'eta'); i = 1;
    elseif strcmp(fileName,'hU'); i = 2;
    else error('Please check the fileName of h, eta and hU')
    end
    % boundary type vector for [h; hU]* [bound1, bound2, bound3]
    for m = 1:length(BoundType_Vec) %m: number of bound types
        ind = (Cell_Bound_ID(:,2) == m); % find the the corresponding bound ID 
        BondType_Vector_Print(ind,:) = repmat(BoundType_Vec{m}(i,:),sum(ind),1);
    end    
end
BondType_Vector_Print = [Cell_Bound_ID(:,1), BondType_Vector_Print];

%% print files
% if the value in 'h.dat' or 'hU.dat' is pure nil,add a small positive value to the first element
if strcmp(fileName(1),'h') && sum(id_IV(:,2))==0 
    id_IV(1,2) = adjustValue;
end
fileID = fopen([writePath fileName '.dat'],'w');
% print valid cell and their initial value
fprintf(fileID,'$Element Number\n'); fprintf(fileID,'%d\n',length(id_IV));
fprintf(fileID,'$Element_id  Value\n'); fprintf(fileID,formatSpec1,id_IV');
% print boundary cell and their vector
fprintf(fileID,'$Boundary Numbers\n'); fprintf(fileID,'%d\n',bound_cnt);
fprintf(fileID,'$Element_id Boundary_type\n'); fprintf(fileID,formatSpec2,BondType_Vector_Print');
fclose(fileID);
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
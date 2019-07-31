function InputSetup(caseFolder, Z, R, varargin)
%INPUTSETUP     set up input files for HIPIMS (one mesh file and 15 types
%   of field files)
%   InputSetUp(caseFolder, Z, R) set up input files of a case, all other
%   parameters are default value. caseFolder is the location of input
%   folder. Z is the DEM value matrix. R is the DEM spatial reference.
%   InputSetUp(caseFolder, Z, R, Name, Value) 
%   Name-VAlue Pair Arguments
%       'h_Eta' -- identification of writing  whether water depth or elevation files
%           'h'(default)|'eta'
%************************Boundary Condition***********************
%       'IO_BoundFrame' -- coordinates of the IO bound frame
%           [](default)|n*4 numeric matrix, each row gives the extent of
%           one bound frame
%       'BoundType' -- Type of Boundary, corresponding to BoundCode. If
%           both BoundType and BoundCode are listed in the input parameters,
%           then BoundType is ignored.
%           'open'(default)| string: 'open','rigid','hgiven','Qgiven','hQgiven'
%       'h_BC_Source' -- h boundary condition source data
%           {[0 0]}(default)|cell matrix consists of certain number of
%           2-column numeric matrixs
%       'hU_BC_Source' -- hU boundary condition source data
%           {[0 0 0]}(default)|cell matrix consists of certain number of 
%           3-column numeric matrix for discharge per unit width(m2/s)
%           2-column numeric matrix for total discharge in one bound (m3/s)
%       'BoundCode' -- boundary code. If BoundType is given by input
%           paramters, BoundCode is unnecessary.
%           [2 0 0; 2 1 0](default)|2*3n numeric matrix, n is the number of bounds
%************************Initial Condition*********************************
%       'initial_hE' -- initial value of h/eta
%           0(default)| scalar/matrix with the same size of Z
%       'initial_hU' -- initial value of hU
%           {0 0}(default)| two-element cell consists of scalar/matrix with
%           the same size of Z
%       'initial_hE_hU_pre' -- initial value of h/eta, hU, and precipitation
%           {0, {0 0}, 0}(default)| 3-element cell: scalar/matrix with the
%           same size of Z. When it is included in the input paramters, then
%           initial_hE and initial_hU will be ignored if they are also
%           listed in the input parameters.
%************************Rainfall Source***********************************
%       'RainMask' -- rainfall mask: the No. of rainfall source
%           0(default)|integral matrix with the same size of Z
%       'RainSource' -- rainfall source values
%           {[0 0]}(default)| cell consists of two-column numeric matrices
%************************Gauges Points*************************************
%       'GaugeCoor' -- coordinates of gauge points
%           [](default)| Two-column(X Y) numerical matrix
%************************Hydrolical Parameters*****************************
%       'manning' -- manning coefficient
%           0.035(default)| scalar/matrix with the same size of Z
%       'sewer_sink' -- sewer_sink coefficient
%           0(default)| scalar/matrix with the same size of Z
%       'cumul_depth' -- cumul_depth coefficient
%           0(default)| scalar/matrix with the same size of Z
%       'hydraulic_conductivity' -- hydraulic_conductivity coefficient
%           0(default)| scalar/matrix with the same size of Z
%       'capillary_head' -- capillary_head coefficient
%           0(default)| scalar/matrix with the same size of Z
%       'water_content_diff' -- water_content_diff coefficient
%           0(default)| scalar/matrix with the same size of Z
%************************Writing Flags*****************************
%       'WriteAllFiles' -- to indicate whether to write all files
%           false(default)| logical value: 'true' write all input files
%***************************not recommended parameters*********************
%       'initial_hE_hU_pre' -- initial value of h/eta, hU, and precipitation
%           {0, {0 0}, 0}(default)| 3-element cell: scalar/matrix with the
%           same size of Z. When it is included in the input paramters, then
%           initial_hE and initial_hU will be ignored if they are also
%           listed in the input parameters.
%       'hydro_params_Value' -- hydrological parameters value: 1-manning,
%               2-sewer_sink, 3-cumul_depth, 4-hydraulic_conductivity,
%               5-capillary_head, 6-water_content_diff. When
%               hydro_params_Value is included in the input parameters, all
%               the six parameters will be ignored if they are also listed
%               in the input parameters.
%           {0.035,0,0,0,0,0}(default)| 4-element cell: scalar or 
%               matrix with the same size of Z
%   See also FieldSetUp, gen_VB_ID, WriteInitialValue.
%   Created by Xiaodong Ming on 03 Aug 2017. Based on FieldSetup
%   Updated by Xiaodong Ming on 28 Aug 2018. 

%% location and write path
% SourceValueCell, MaskValueCell, InitialValueCell, h_eta_string
if nargin < 3
    error('InputSetup: Wrong Number Of Inputs')
end
if ~exist([caseFolder '/input'],'dir')
    mkdir([caseFolder '/input']);
    mkdir([caseFolder '/input'],'field');
    mkdir([caseFolder '/input'],'mesh')
end
if ~exist([caseFolder '/input/field'],'dir')
    mkdir([caseFolder '/input/field']);
end
if ~exist([caseFolder '/input/mesh'],'dir')
    mkdir([caseFolder '/input/mesh']);
end
writePath = [caseFolder '/input/field/'];
Delimiter = ' ';
%% define default input param names and values
% Name-Value pairs: 11
allParamNames = {'h_Eta', ... %1
         'IO_BoundFrame', 'BoundCode',...%2,3
         'h_BC_Source', 'hU_BC_Source',...%4,5
         'initial_hE_hU_pre', ...%6
         'hydro_params_Value',...%7
         'RainMask', 'RainSource',...%8,9
         'GaugeCoor',...%10
         'WriteAllFiles',...%11
         'initial_hE', 'initial_hU',...%12,13
         'manning',...%14
         'sewer_sink',...%15
         'cumul_depth','hydraulic_conductivity','capillary_head','water_content_diff',...%16,17,18,19
         'BoundType'}; %20
originParamValues  = {'h',...%1:h_Eta
          [], [2 0 0; 2 1 0],...%2,3:'IO_BoundFrame','BoundCode'
          {[0 0; 3600 0]}, {[0 0 0; 3600 0 0]},...%4,5:'h_BC_Source','hU_BC_Source'
          [],...%6:initial_hE_hU_pre
          [],...%7:hydro_params_Value
          0, {[0 0; 3600 0]},...%8,9:'RainMask', 'RainSource'
          [],...%10: GaugeCoor
          false,...%11: WriteAllFiles
          0,{0 0},...%12,13:'initial_hE','initial_hU'
          0.035,...%14: manning
          0,...%15: sewer_sink
          0,0,0,0,...%16,17,18,19:'cumulative_depth','hydraulic_conductivity','capillary_head','water_content_diff'
          'open'}; % 20 bound type

%% Parameters passing
n = length(varargin)/2; % number of Name-Value pairs from the input parameter
if n~=fix(n)
    error('Wrong number of varargin parameters')
end
paramInputFlag = zeros(size(allParamNames)); % to check which param is input
if n > 0
    inputParamNames = varargin(1:2:n*2-1); % all the input parameter names
    inputParamValues = varargin(2:2:n*2);  % the values of corresponding input parameters
    for i = 1:n
        param_ind = strcmp(inputParamNames{i},allParamNames);
        originParamValues{param_ind} = inputParamValues{i};
        paramInputFlag(param_ind) = 1;
    end
end
if isempty(originParamValues{6})
    originParamValues{6} = {originParamValues{12},originParamValues{13},0};
end
if isempty(originParamValues{7})
    originParamValues{7} = originParamValues(14:19);
end
%% writting file flags
% field files 15
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
% flags about whether to write the 15 field files
if isnumeric(originParamValues{11}) % writting flag input argument
    writtingFlag = originParamValues{11};
elseif originParamValues{11}==true
    writtingFlag = ones(15,1);
else
    writtingFlag = zeros(15,1);
    writtingFlag(1) = paramInputFlag(4);
    writtingFlag(2) = paramInputFlag(5);
    writtingFlag(3:4) = max(paramInputFlag(6),paramInputFlag(12:13));
    writtingFlag(6:11) = max(paramInputFlag(7),paramInputFlag(14:19));
    writtingFlag(12) = paramInputFlag(8);
    writtingFlag(13) = paramInputFlag(9);
    writtingFlag(14) = paramInputFlag(10);
end
%% Initialize internal variables
%write water depth(h) or water level(eta) file
h_Eta = originParamValues{1};
% coordinates of IO Bound Frame
IO_BoundFrame = originParamValues{2}; % coordinates of IO Bound Frame
if ~isempty(IO_BoundFrame)
    if size(IO_BoundFrame,2)~=4
        error('Wrong IO_bound frame coordinates')
    end
    numOfBound = size(IO_BoundFrame,1)+1;
else
    numOfBound = 1;
end
% Create Bound Code
bound_Code_Cell = cell(1,numOfBound);
givenMark_h = 0;
givenMark_U = 0;
if paramInputFlag(20) == 1 % BoundType is input from outside
    boundTypeStr = originParamValues{20};
    if ~iscellstr(boundTypeStr)
        if ischar(boundTypeStr)
            boundTypeStr = {boundTypeStr};
        else
            error('Wrong boundType')
        end
    end
    if numOfBound ~= length(boundTypeStr)
        error('The number of boundType is not in accordance with the rows of IO_Bound Frame')
    end
    for i=1:numOfBound
        boundTypeStr_one = boundTypeStr{i};
        switch boundTypeStr_one
            case 'open'
                bound_Code_Cell{i} = [2 0 0; 2 1 0];
            case 'rigid'
                bound_Code_Cell{i} = [2 0 0; 2 2 0];
            case 'hgiven'
                bound_Code_Cell{i} = [3 0 givenMark_h; 2 1 0];
                givenMark_h = givenMark_h+1;
            case 'Qgiven'
                bound_Code_Cell{i} = [2 0 0; 3 0 givenMark_U];
                givenMark_U = givenMark_U+1;
            case 'hQgiven'
                bound_Code_Cell{i} = [3 0 givenMark_h; 3 0 givenMark_U];
                givenMark_h = givenMark_h+1;
                givenMark_U = givenMark_U+1;
            otherwise
                error(['cannot recognize string ' boundTypeStr_one ' from input BoundType'])
        end
    end
end

if paramInputFlag(3)==1 % boundCode imported while boundType is not
    bound_Code = originParamValues{3};
    if paramInputFlag(20)==1
        warning('input param BoundType is noneffective as BoundCode is also imported')
    end
else
    bound_Code = cell2mat(bound_Code_Cell);
end
% boundary source data of water depth/level
h_BC_source_cell  = originParamValues{4};
if ~iscell(h_BC_source_cell)
    h_BC_source_cell = {h_BC_source_cell};
end
% boundary source data of water velocity
hU_BC_source_cell = originParamValues{5};
if ~iscell(hU_BC_source_cell)
    hU_BC_source_cell = {hU_BC_source_cell};
end
% initial h, vel, and precipitation value
initial_hE_hU_pre = originParamValues{6};

% hydro params value (1-manning, 2-sewer_sink, 3-cumulative_depth,
%   4-hydraulic_conductivity, 5-capillary_head, 6-water_content_diff;
hydro_params_Value = originParamValues{7}; % 6-element cell
hydro_params_Names = {'manning','sewer_sink','cumulative_depth','hydraulic_conductivity','capillary_head','water_content_diff'};

% rainfall mask and rainfall source data
rainMask = originParamValues{8};
if or(size(rainMask)==size(Z),numel(rainMask)==1)
    rainMask = zeros(size(Z))+rainMask;
else
    error('Wrong Rain mask value')
end
rainSourceValue = originParamValues{9};
%% generate valid ID and boundary ID for each grid cell 
% invoke sub function "gen_VB_ID"
[Valid_ID,Bound_ID] = gen_VB_ID(Z,R,IO_BoundFrame);%
% disp('Valiad grid cell ID and Boundary grid cell ID have been generated')
%% ****** write files *********
%% * write mesh file: DEM.txt 
%***output file No. 0: mesh/DEM.txt
demName = [caseFolder '/input/mesh/DEM.txt'];
if originParamValues{11} %write_all_files is true 
    arcgridwrite(demName,Z,R)
    disp('FileNo-00: mesh file DEM.txt created')
end
%% * write boundary condition files: h/eta_BC.dat and hU_BC.dat 
%***output file No. 1,2: h_BC_num.dat, hU_BC_num.dat
% number of bounds in total

if  and(numOfBound==1,isempty(bound_Code))
    bound_Code = originParamValues{3};%No Bound_Code for the outline boundary cell: set as open boundary 
else
    if size(bound_Code,2)==3&&numOfBound>1
    bound_Code = repmat(bound_Code,[1,numOfBound]);
    warning('No Bound_Code for some of your boundaries: set the same as the outline boundary')
    end
end
% the first values of Bound_Attributes_Code for each bound
h_BC_Code_1  = bound_Code(1,1:3:3*numOfBound);
hU_BC_Code_1 = bound_Code(2,1:3:3*numOfBound);
h_BC_Num  = 1;
hU_BC_Num = 1;
for i = 1:numOfBound % bound ID
    % judge whether it is a IO bound
    if h_BC_Code_1(i)  == 3
        h_BC_source_num  = bound_Code(1,3+3*(i-1));
        h_BC_index_name  = [h_Eta '_BC_' num2str(h_BC_source_num) '.dat'];
        if length(h_BC_source_cell) < h_BC_Num
            error('lack of source matrix for h_BC')
        end
        h_BC_source_value = h_BC_source_cell{h_BC_Num}; h_BC_Num = h_BC_Num+1;
        if size(h_BC_source_value,2)~=2
            warning(['Wrong boundary condition value of hEta for Bound ' num2str(i)])
        end
        % output file No 1: h_BC_num.dat
        if writtingFlag(1)==1
            dlmwrite([writePath h_BC_index_name],h_BC_source_value,...
                'precision',10,'delimiter',Delimiter)
            disp([h_BC_index_name ' for Bound' num2str(i)])
        end
    end
    if hU_BC_Code_1(i) == 3
        hU_BC_source_num = bound_Code(2,3+3*(i-1));
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
        if writtingFlag(2)==1
            dlmwrite([writePath hU_BC_index_name],hU_BC_source_value,...
                'precision',10,'delimiter',Delimiter)
            disp([hU_BC_index_name ' for Bound' num2str(i)])
        end
    end
end
if writtingFlag(1)==1
    disp(['FileNo-01: h_BC_num has been created, total number of h source file: ' num2str(h_BC_Num-1)])
end
if writtingFlag(2)==1
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
        if and(size(write_value{1})==size(Z)|numel(write_value{1})==1,...
                size(write_value{2})==size(Z)|numel(write_value{2})==1)
            write_value = [zeros(size(Z))+write_value{1}, zeros(size(Z))+write_value{2}]; % initial hU discharge per unit width
        else
            error(['Wrong initial ' hE_hU_pre_name{i} ' data: does not match with the size of DEM matrix'])
        end
    else
        if or(size(write_value)==size(Z),numel(write_value)==1)
            write_value =  zeros(size(Z)) + write_value; %initial h water depth or eta water surface elevation
        else
            error(['Wrong initial ' hE_hU_pre_name{i} ' data: does not match with the size of DEM matrix'])
        end
    end
    if writtingFlag(2+i)==1
        if strcmp(hE_hU_pre_name{i},'hU')
            if size(write_value)==size(Z)
                write_value = repmat(write_value,[1 2]);
            end
        end
        WriteInitialValue(writePath,hE_hU_pre_name{i},Valid_ID,Bound_ID,...
            write_value,bound_Code)
        disp(['FileNo-' num2str(2+i,'%02u') ': ' hE_hU_pre_name{i} '.dat created'])
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
    if and(numel(param_value)~=1,size(param_value)~=size(Z))
        error(['Wrong ' param_name ' value'])
    end
    grid_param_Value = zeros(size(Z))+param_value;
    if writtingFlag(5+i)==1
        WriteInitialValue(writePath,param_name,Valid_ID,Bound_ID,grid_param_Value,bound_Code)
        disp(['FileNo-' num2str(5+i,'%02u') ': ' param_name '.dat created'])
    end
end

%% write Precipitation Mask file
% precipitation_mask.data
VID_RMask = [Valid_ID(:) rainMask(:)]; % zeros(numel(Valid_ID),1) rainfall mask value is 0 
VID_RMask(isnan(Valid_ID(:))|Valid_ID(:)==-1,:)=[]; % delete nodata value
VID_RMask = sortrows(VID_RMask,1);
% write file
%***output file No.12: precipitation_mask.dat 
if writtingFlag(12)==1
    filename = 'precipitation_mask.dat';
    fileID = fopen([writePath filename],'w');
    fprintf(fileID,'$Element Number\n'); fprintf(fileID,'%d\n',length(VID_RMask));
    fprintf(fileID,'$Element_id  Value\n'); fprintf(fileID,'%-12d %d\n',VID_RMask');
    fclose(fileID);
    disp(['FileNo-12: ' filename ' created'])
end
%% write Precipitation Source files: 
%***output file No.13: precipitation_source_num.dat
numOfRainSource = 1;
if writtingFlag(13)==1
    if iscell(rainSourceValue)||size(rainSourceValue,2)==2
        numOfRainSource = length(rainSourceValue);
        for i = 1:numOfRainSource
            Time_Rain = rainSourceValue{i};
            Time_Rain(isnan(Time_Rain(:,2)),:)=[];
            filename = ['precipitation_source_' num2str(i-1) '.dat'];
            fileID = fopen([writePath filename],'w');
            fprintf(fileID,'%09.2f %.15e\n',Time_Rain');
            fclose(fileID);
        end
        disp(['FileNo-13: ' num2str(i) ' precipitation_source files created'])
    else
        numOfRainSource = size(rainSourceValue,2)-1;
        formatSpec = repmat('%.15e ',[1 numOfRainSource] );
        formatSpec = ['%09.2f ' formatSpec(1:end-1) '\n'];
        fileID = fopen([writePath 'precipitation_source_all.dat'],'w');
        fprintf(fileID,'%u\n',numOfRainSource);
        fprintf(fileID,formatSpec,rainSourceValue');
        fclose(fileID);
        disp('FileNo-13: precipitation_source_all created')
    end
end
rainMaskNum = unique(rainMask(:)); rainMaskNum(isnan(rainMaskNum))=[];
if numel(rainMaskNum) > numOfRainSource
    warning('The number of rainfall sources is smaller than the number of rain mask values.')
end
%% *****gauges_pos.dat
%coordinate of gauge
%***output file No.14:gauges_pos.dat
GaugeCoor = originParamValues{10};
if ~isempty(GaugeCoor)
    if size(GaugeCoor)~=2
        error('Wrong dimension of GaugeCoor')
    end
end
filename = 'gauges_pos.dat';
if writtingFlag(14)==1
    fileID = fopen([writePath filename],'w');
    if isempty(GaugeCoor)
        GaugeCoor = [0 0];
    end
    fprintf(fileID,'%.3f %.3f\n',GaugeCoor(1:end-1,:)');
    fprintf(fileID,'%.3f %.3f',GaugeCoor(end,:));
    fclose(fileID);
%     dlmwrite([writePath filename],GaugeCoor,'precision',10,'delimiter',Delimiter);
    disp(['FileNo-14: ' filename ' created'])
end
%% **** z.dat
if writtingFlag(15)==1
   WriteInitialValue(writePath,'z',Valid_ID,Bound_ID,Z,bound_Code)
   disp('FileNo-15: z.dat created')
end
%% *****summary information
if min(writtingFlag)==1
    disp('Mission accomplished: all input files have been created')
%     figure;
%     fig_hd = mapshow(Z,R,'DisplayType','surface');
%     for NOB = 1:NumOfBound-1
%         rect_CoorX = IO_BoundFrame(NOB,[1 3]);
%         rect_CoorY = IO_BoundFrame(NOB,[2 4]);
%         rec_hd = rectangle('Position',[min(rect_CoorX),min(rect_CoorY),range(rect_CoorX),range(rect_CoorY)]);
%         rec_hd.LineWidth = 3; rec_hd.EdgeColor = [1 0 0];% 1 0 0: red
%     end
%     zdatam(fig_hd,Z-max(Z(:)))
end
end

function arcgridwrite(filename,Z,R)
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
% x00 and y00 specify the centre of the upper-left fake pixel of (x11,y11)
x00 = R(3,1); 
y00 = R(3,2);
gridSize = abs(R(2));
nrows = size(Z,1);
ncols = size(Z,2);

% XLLCORNER and YLLCORNER are the coordinates of the
% lower left corner of the lower left cell of Z.
xllcorner = x00 + gridSize/2;
yllcorner = y00 - gridSize*nrows - gridSize/2;
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

% avoid the depth in the flow-input-bound cells to be nil
if or(strcmp(fileName,'h'),strcmp(fileName,'eta')) 
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

if strcmp(fileName,'h')||strcmp(fileName,'eta')||strcmp(fileName,'hU')
    if or(strcmp(fileName,'h'), strcmp(fileName,'eta'))
        i = 1;
    elseif strcmp(fileName,'hU')
        i = 2;
    else
        error('Please check the fileName of h, eta and hU')
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
if strcmp(fileName,'h')&&sum(id_IV(:,2))==0 
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
    for j = 1:ncols
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
    for j = 1:ncols
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
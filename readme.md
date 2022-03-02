Description of InputSetup function. There is another object-oriented function [HiPIMS_Setup](class/HiPIMS_Setup.m) available.
# Function
***INPUTSETUP*** is a MATLAB function to set up the input folders and files for HiPIMS
# List of input files
## Mesh file: 
DEM.txt	spatial reference of the DEM raster of the model domain
## Field files:
There are typically 15 field files if one boundary condition was given for depth and velocity, respectively. More field files would be created if there were more boundary conditions. The fundamental files are introduced in Table 1.

*Table 1. List of field files*

|**File name**|**File type**|**Notes**|
| :- | :- | :- |
|**z.dat**|Basic terrain|ID, type, and elevation values of each cell|
|**h.dat/eta.dat**|Initial conditions|initial water depth/elevation (m) of each cell|
|**hU.dat**||initial water velocities (m/s) of each cell|
|**precipitation.dat**||initial rainfall rate (m/s) of each cell|
|**manning.dat**|Parameters|friction coefficient of each cell|
|**hydraulic\_conductivity.dat** ||hydraulic conductivity of each cell|
|**cumulative\_depth.dat**||cumulative depth of each cell|
|**capillary\_head.dat**||capillary head of each cell|
|**water\_content\_diff.dat**||water content difference of each cell|
|**sewer\_sink.dat**||sewer sink rate of each cell|
|**h\_BC\_0.dat**|Boundary conditions|water depth on the first boundary|
|**hU\_BC\_0.dat**||water velocities on the first boundary|
|**precipitation\_mask.dat**|Rainfall sources|rainfall source ID of each cell|
|**precipitation\_source\_all.dat**||time series of all rainfall sources|
|**gauges\_pos.dat**|Monitors|coordinates of monitoring points in the model|

# Calling format
```
InputSetUp(caseFolder,Z,R)
```
set up input files for a case. *caseFolder* is the location of the folder storing input and output files. *Z* is a matrix storing elevation values of DEM. *R* is a 3\*2 matrix of DEM spatial-reference information, including the coordinate of the original points of the raster and the size of the grid. *Z* and *R* could be created separately (*makerefmat*) or read from existing Arc ascii files via *arcgridread*. All other parameters of HiPIMS are default values if the input parameters are as listed above. 

```
InputSetUp(caseFolder, Z, R, Name, Value)
```
caseFolder* is the location of the folder storing input and output files. Z is a matrix of elevation value. R is a spatial-reference matrix of DEM. Name-Value Pair Arguments are listed as Table 2.

*Table 2 Name-Value pair arguments*

|**Parameter Type**|**Name(Case sensitive)**|**Values(Default \| Alternative)**|**Format of Values**|**Note**|
| :- | :- | :- | :- | :- |
|**Decision Flags**|h\_Eta|'h'\|'eta'|string|To decide whether the initial conditions are given as water depth (h) or water elevation (eta)|
||WriteAllFiles|false\| true|logical|Whether to generate all the input files|
|Gauge coordinates|GaugeCoor|[ ]|a numeric array (2-column) |Coordinates of the gauge points inside domain|
|**Boundary conditions**|IO\_BoundFrame|[ ]|4\*n numeric array|Extent of the input-output boundaries. *n* is the number of IO boundaries|
||BoundType|'open'\|'rigid', 'hgiven', 'Qgiven', 'hQgiven'|a string or a cell of multiple strings |Type of boundaries. ‘*hgiven’* means water depth /elevation in the bound is pre-defined; '*Qgiven'* means the discharge/water velocity in the bound is pre-defined; *'hQgiven'* means both depth and discharge in the bound is predefined.|
||h\_BC\_Source|{[0 0]}|a cell of 2-column numeric arrays|Data of pre-defined water depth/elevation. The number of 2- column arrays should be the same with the number of boundaries that h/eta has been given. The first column of the array is time(s) and the second column is the water depth/elevation (m).|
||hU\_BC\_Source|{[0 0 0]}|a cell of 2/3-column numeric arrays|Data of pre-defined water discharge/ water velocity. The number of numeric arrays should be the same with the number of boundaries that discharge/velocity has been given. The first column of the array is time(s) and if the array is 2-col, the second column is the discharge (m3/s) or if the array is 3-col, then the second and third column are water velocity (m/s) in x and y direction respectively.|
||BoundCode (not recommended)|[2 0 0; 2 1 0]|2\*3n numeric array|Not recommended unless the alternative bound types cannot fulfil your requirements. It conveys more specific information of BoundType with numeric arrays. |
|**Initial conditions**|initial\_hE|0|scalar or numeric array with the same size of Z|Initial water depth/elevation. If it is a scalar, then all the grids in the domain have the same initial h/eta value.|
||initial\_hU|{0 0}|scalar (0) or cell of two numeric arrays with the same size of Z|Initial water velocity. Two components of the cell represent initial velocity in x and y direction respectively. If it is a scalar, then all the grids in the domain have the same initial water velocity value in both x and y direction.|
||initial\_hE\_hU\_pre (not recommended) |{0, {0 0}, 0}|cell of three numeric arrays|It is a combination of all initial conditions, including initial h, hU and precipitation. The last one (precipitation) is always 0 at current version.|
|**Rainfall**|RainMask|0|scalar or numeric array with the same size of Z|It is the serial number of rainfall source starting from 0. Grids with the same serial number will have the same rainfall from the same source.|
||RainSource|[0 0; 3600 0]|a numeric array or a cell of 2-column numeric arrays|To give rainfall value for different region of the domain.  If it is a numeric array, the first column is time(s) and the second and right forward columns are the rainfall rate(m/s), and the output rainfall source file will be a single file named ‘precipitation\_source\_all.dat’. The number of the single array column should be in accordance with the number of rainfall source in RainMask. If it is a cell of 2-column numeric arrays, each array conveys the time (s, 1st column) and rainfall rate (m/s, 2nd column) of one single rainfall source. Multiple files of rainfall source will be generated and named as ‘precipitation\_source\_n.dat’. The number of 2-column numeric arrays should be in accordance with the number of rainfall source.|
|**Hydro Parameter Values**|manning|0.035|scalar or numeric array with the same size of Z|It is manning coefficient. If it is a scalar, then all the grids in the domain have the same manning value.|
||sewer\_sink|0|scalar or numeric array with the same size of Z|It is sewer sink rate (m/s). If it is a scalar, then all the grids in the domain have the same sewer sink value.|
||cumul\_depth|0|scalar or numeric array with the same size of Z|It is one of the infiltration parameters. If it is a scalar, then all the grids in the domain have the same cumulative depth value.|
||hydraulic\_conductivity|0|scalar or numeric array with the same size of Z|It is one of the infiltration parameters. If it is a scalar, then all the grids in the domain have the same hydraulic conductivity value.|
||capillary\_head|0|scalar or numeric array with the same size of Z|It is one of the infiltration parameters. If it is a scalar, then all the grids in the domain have the same capillary head value.|
||water\_content\_diff|0|scalar or numeric array with the same size of Z|It is one of the infiltration parameters. If it is a scalar, then all the grids in the domain have the same water content diff value.|
||hydro\_params\_Value (not recommended)|{0.035, 0, 0, 0, 0, 0}|scalar or numeric array with the same size of Z|It is a combination of all the six hydro parameter parameters.|
# Example
```
%% Example to create input files based on a peaks DEM
%% creat a DEM with Z and R
Z = peaks(500); % elevation values of DEM
gridL = 1; % length of each square grid
x11 = 0; % coordinates of the center of the upper left point
y11 = (size(Z,1)-1)*gridL;
R = makerefmat(x11,y11,gridL,-gridL); %spatial reference of DEM
mapshow(Z,R,'DisplayType','surface'); 
colorbar; 
box on;
title('DEM');
xlabel('meter towards east');
ylabel('meter towards north');
```

![Domain map](https://github.com/mingxiaodong/hipims_io_matlab/blob/master/doc/domain_map01.jpeg)

```
%% define boundary condition
% outline boundary
outlineBoundType = 'open';
% coordinates of the end row/col of Z
x_end = x11+(size(Z,2)-1)*gridL; 
y_end = y11+(size(Z,1)-1)*(-gridL);
% input-output boundary 1
IO_Bound1_Frame = [x11-2*gridL 2*y11/5, x11+2*gridL 3*y11/5];
IO_Bound1_Type = 'Qgiven';
dischage = [0 30; 3600 300; 7200 30; 10800 30];
% input-output boundary 2
IO_Bound2_Frame = [x_end-2*gridL 2*y11/5, x_end+2*gridL 3*y11/5];
IO_Bound2_Type = 'hgiven';
depth = [0 1; 3600 3; 7200 1; 10800 1];
IO_BoundFrame = [IO_Bound1_Frame; IO_Bound2_Frame];
boundType = {outlineBoundType,IO_Bound1_Type,IO_Bound2_Type};
h_source = depth;
Q_source = dischage;
% show the IO bound frames
mapshow(Z,R,'DisplayType','surface');
box on;
axis off
rectangle('Position',[x11-2*gridL 2*y11/5 gridL*4 y11/5],'EdgeColor','r')
rectangle('Position',[x_end-2*gridL 2*y11/5 gridL*4 y11/5],'EdgeColor','r')
```

![boundary map](https://github.com/mingxiaodong/hipims_io_matlab/blob/master/doc/domain_map02.jpeg)


```
%% define rainfall condition
% rainfall mask: two rainfall source, north(0) and south(1)
rainMask = zeros(size(Z)); rainMask(round(size(Z,1)/2):end,:) = 1;
rainSource = [0   ,  0, 100/3600/4;...
              3600,  0, 200/3600/4;...
              7200,  0, 100/3600/4;...
              7201,  0, 100/3600/4];
%% generate input files
caseFolder = cd;
InputSetup(caseFolder, Z, R,...
        'IO_BoundFrame',IO_BoundFrame,'BoundType',boundType,... 
        'h_BC_Source',h_source,...
        'hU_BC_Source',Q_source,...
        'RainMask',rainMask,'RainSource',rainSource,...
        'WriteAllFiles','true');
```


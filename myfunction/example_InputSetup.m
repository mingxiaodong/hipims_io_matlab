%% Example to create input files based on a peaks DEM
%% creat a DEM with Z and R
Z = peaks(100); % elevation values of DEM
gridL = 1; % length of each square grid
x11 = 0; % coordinates of the center of the upper left point
y11 = (size(Z,1)-1)*gridL;
R = makerefmat(x11,y11,gridL,-gridL); %spatial reference of DEM
mapshow(Z,R,'DisplayType','surface'); colorbar; box on;
title('DEM'); xlabel('meter towards east'); ylabel('meter towards north');
%% define gauges position data
gaugesCoor = [(0:10:100)' (0:10:100)'];
%% define boundary condition
% outline boundary
outlineBoundType = 'hQgiven';
% coordinates of the end row/col of Z
x_end = x11+(size(Z,2)-1)*gridL; 
y_end = y11+(size(Z,1)-1)*(-gridL);
% input-output boundary 1
IO_Bound1_Frame = [x11-2*gridL 2*y11/5, x11+2*gridL 3*y11/5];
IO_Bound1_Type = 'Qgiven'; %dischage is pre-defined
dischage1 = [0 30; 3600 300; 7200 30; 10800 30];
dischage2 = [0 30; 3600 300; 7200 30; 10800 30];
% input-output boundary 2
IO_Bound2_Frame = [x_end-2*gridL 2*y11/5, x_end+2*gridL 3*y11/5];
IO_Bound2_Type = 'hgiven'; %water depth is pre-defined
depth1 = [0 1; 3600 3; 7200 1; 10800 1];
depth2 = [0 1; 3600 3; 7200 1; 10800 1];
% show the IO bound frames
mapshow(Z,R,'DisplayType','texturemap');box on; axis off
rectangle('Position',[x11-2*gridL 2*y11/5 gridL*4 y11/5],'EdgeColor','r')
rectangle('Position',[x_end-2*gridL 2*y11/5 gridL*4 y11/5],'EdgeColor','r')

IO_BoundFrame = [IO_Bound1_Frame; IO_Bound2_Frame];
boundType = {outlineBoundType,IO_Bound1_Type,IO_Bound2_Type};
h_source = {depth1,depth2};
Q_source = {dischage1,dischage2};
%% define rainfall
% rainfall mask: two rainfall source, north(0) and south(1)
rainMask = zeros(size(Z)); rainMask(round(size(Z,1)/2):end,:) = 1;
rainSource = [0   ,  0, 100/3600/4;...
              3600,  0, 200/3600/4;...
              7200,  0, 100/3600/4;...
              7201,  0, 100/3600/4]; %unit m/s
%% generate input files
caseFolder = cd;
InputSetup(caseFolder, Z, R,...
        'IO_BoundFrame',IO_BoundFrame,'BoundType',boundType,... % boundary code
        'h_BC_Source',h_source,... % boundary source
        'hU_BC_Source',Q_source,... % boundary source
        'RainMask',rainMask,'RainSource',rainSource,... % rainfall
        'GaugeCoor',gaugesCoor,...
        'WriteAllFiles',true);

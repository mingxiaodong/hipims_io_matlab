clear,clc
% Example 1: fake DEM 
Z = peaks(500);
dx = 2; xllcorner = 0; yllcorner = 0;
x11 = xllcorner+dx/2; y11 = yllcorner+size(Z,1)*dx-dx/2;
R = makerefmat(x11,y11,dx,-dx);
obj = HiPIMS_Setup(Z,R);
boundStruct = obj.BoundStruct;
boundStruct(1).Type = 'rigid';
boundStruct(2).Type = 'open';
boundStruct(2).Position = [-10,500; 10,700];
boundStruct(3).Type = 'open';
boundStruct(3).Position = [500,-10; 700,10];
boundStruct(3).T_huv = [0,10; 700,100];
obj = setBoundary(obj,boundStruct);
obj.GaugeCoor = [100,500;300,300];
writeInputFile(obj);
figure(1)
GeneralMap(obj)
figure(2)
BoundID = PlotBoundID(obj);
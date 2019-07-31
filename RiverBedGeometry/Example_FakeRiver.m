% example:
clear,clc
%*** create river bank lines
xCoor = linspace(1,100,100)';
yCoor = cos(xCoor/(2*pi))*2+50;
bankLine1 = [xCoor yCoor];
bankLine2 = [xCoor yCoor+10];
zCoor = linspace(20,5,numel(xCoor))'; % not necessary, only for plotting
% figure % plot 3D line for the two banks
% hold on
% plot3(xCoor,yCoor,zCoor)
% plot3(xCoor,yCoor+10,zCoor)
% hold off
% axis equal
%*** create cross-section lines
y = linspace(61.97,51.97,10)';
z = linspace(-1,1,numel(y))';
z = 5*z.^2-5;
x = y*0+1;
% cross-section lines
CL1 = [x, y, z+20];
CL2 = [x+38, y, z+14.24];
CL3 = [x+78, y, z+8.18];
CL4 = [x+99, y-4, z+5];
crossSectionLine = {CL1,CL2,CL3,CL4};
figure(1)
hold on
% cross-section lines
plot3(CL1(:,1),CL1(:,2),CL1(:,3),'b','LineWidth',1.5)
plot3(xCoor,yCoor,zCoor,'g-','LineWidth',2)
plot3(CL2(:,1),CL2(:,2),CL2(:,3),'b','LineWidth',1.5)
plot3(CL3(:,1),CL3(:,2),CL3(:,3),'b','LineWidth',1.5)
plot3(CL4(:,1),CL4(:,2),CL4(:,3),'b','LineWidth',1.5)
plot3(xCoor,yCoor+10,zCoor,'g-','LineWidth',2)
box on
grid on
ax = gca;
ax.BoxStyle = 'full';
axis image
legend({'Cross-section line','Bank line'})
view([60 8])

%% Generate river bed elevation grid
cellsize = 0.3; % meter
errorDistance = 0.1;% meter
[Z,R] = GenerateBedElevationRaster(bankLine1,bankLine2,cellsize,crossSectionLine,errorDistance);
%%
figure(2)
hold on
% cross-section lines
plot3(CL1(:,1),CL1(:,2),CL1(:,3),'b','LineWidth',1.5)
plot3(xCoor,yCoor,zCoor,'g-','LineWidth',2)
plot3(CL2(:,1),CL2(:,2),CL2(:,3),'b','LineWidth',1.5)
plot3(CL3(:,1),CL3(:,2),CL3(:,3),'b','LineWidth',1.5)
plot3(CL4(:,1),CL4(:,2),CL4(:,3),'b','LineWidth',1.5)
plot3(xCoor,yCoor+10,zCoor,'g-','LineWidth',2)
axis image
view([60 8])
hold off
hold on
mapshow(Z,R,'DisplayType','mesh')
hold off
box on
grid on
ax = gca;
ax.BoxStyle = 'full';
legend({'Cross-section line','Bank line'})

%%
Arcgridwrite('fakeriver.asc',Z,R)
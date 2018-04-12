%% generate a one dimension DEM and initial condition
clear,clc
domainL = 10000;
cellsize = 10;
domainW = cellsize*3;
x11 = -domainL*0.5;
y11 = domainW/2-0.5*cellsize;
dx = cellsize;
dy = -cellsize;
r_dem =makerefmat(x11, y11, dx, dy);

h0 = 10;
a = 3000;
x_vec = -domainL/2:cellsize:domainL/2;
b_vec = h0*((x_vec/a).^2);
t = 0;
B = 5;
g = 9.8;
s = sqrt(8*g*h0/a^2/2);
eta_vec = h0 - a^2*B^2*s^2*cos(2*s*t)/(8*g^2*h0) - B^2/(4*g)- B*s*x_vec/g;
h_vec = eta_vec-b_vec; h_vec(h_vec<0) = 0;
z_dem = repmat(b_vec,[3,1]);
h_grid = repmat(h_vec,[3,1]);
gauges_pos = [x_vec',x_vec'*0+y11-cellsize];
% plot(x_vec,b_vec,x_vec,h_vec)
% legend('bed','\eta value')
%% write file
fileName = 'G:\cudaSWEsSolver\WRRTest1\dem.txt';
arcgridwrite(fileName,z_dem,r_dem)
fileName = 'G:\cudaSWEsSolver\WRRTest1\h_initial.asc';
arcgridwrite(fileName,h_grid,r_dem)
fileName = 'G:\cudaSWEsSolver\WRRTest1\gauges.txt';
dlmwrite(fileName,gauges_pos)
%%
figure;
hold on
mapshow(z_dem,r_dem,'DisplayType','surface');
scatter3(gauges_pos(:,1),gauges_pos(:,2),gauges_pos(:,2)*0+20,'r')
axis square
hold off
%% modify the elevation in boundary cells of the DEM file
clear,clc
%cd('G:\cudaSWEsSolver\ThamesValley\input\mesh')
cd('G:\Data')
filename = 'emodnet-mean.asc';
delimiter = ' ';
endRow = 6;
formatSpec = '%s%f';
fileID = fopen(filename,'r');
fileRef = textscan(fileID, formatSpec, endRow, 'Delimiter', delimiter,'MultipleDelimsAsOne', true);
fclose(fileID);
RefValues = fileRef{2};
[Z,R] = arcgridread(filename);
currentPath = cd;
folder_path = 'C:\Users\b4042552\Dropbox\Matlab\GeoClasses';%;%
%folder_path = '/Users/Sheldon1/Dropbox/Matlab/GeoClasses';
cd(folder_path)
CaseName = 'ThamesValley';
writePath = ['G:\cudaSWEsSolver\' CaseName '\input\field\'];
Initial_XLS = [writePath(1:end-12) CaseName '_InitialData.xlsx'];
IO_Bound = xlsread(Initial_XLS,'BoundFrame','B:E');
[Valid_ID,Bound_ID] = gen_VB_ID(Z,R,IO_Bound);%  
cd(currentPath)
Z_N = Z;
Z_N(isnan(Z))= RefValues(end);
Z_N(Z>200)=RefValues(end);
Z_N(Z<0)=RefValues(end);
fileID = fopen('DEM.txt','w');
for i=1:6
    forspec = '%-12s  %d\n';
    if i==3||i==4; forspec = '%-12s  %.4f\n'; end
    fprintf(fileID,'%-12s  %d\n', fileRef{1}{i}, RefValues(i));
end
dlmwrite('DEM.txt',Z_N,'-append','delimiter',delimiter)
fclose(fileID);

%% generate smooth bowl DEM and initial water depth
clear,clc
dx = 50;
dy = -dx;
xlength = 10000;
ylength = 10000;
ncols = round(xlength/abs(dx));
nrows = round(ylength/abs(dy));
x11 = -5000;
y11 = -5000;
R = makerefmat(x11, y11+ylength, [0 dx], [dy 0]);
H = zeros(nrows,ncols);
[Rows, Cols] = meshgrid(1:nrows,1:ncols);
[X,Y] =  pix2map(R, Rows(:), Cols(:));
a = 3000; h0 = 10; B = 5;
Z = h0*(X.^2+Y.^2)/a^2;

% a = 0.5*xlength; b = 0.5*ylength; c = 600;
% Z = -((1-(X-500).^2/a^2-(Y-500).^2/b^2)).^0.5*c;
% ind = (1-(X-500).^2/a^2-(Y-500).^2/b^2)<0;
% Z(ind) = 0;
%g=9.8; s = sqrt(8*g*h0)/(a*2);
%h = h0- 0.5/g*B^2-B/g*s.*X -0.5/g*B^2-B/g*s.*X;
H(:) = Z;
h = h0-H;
h(h<0)=0;
u = H*0;
v = H*0-B; v(h==0) = 0;
hv = v.*h; hu = u.*h;
figure; mapshow(h+H,R,'DisplayType','mesh'); axis square
I_Z = H; I_h = h; I_hU = [hu hv];
%hold on; quiver(u,v,'r'); hold off
%%
fileID = fopen('DEM.txt','w');
fprintf(fileID,'ncols    %d\n', ncols);
fprintf(fileID,'nrows    %d\n', nrows);
fprintf(fileID,'xllcorner    %f\n', x11);
fprintf(fileID,'yllcorner    %f\n', y11-ylength);
fprintf(fileID,'cellsize    %f\n', abs(dx));
fprintf(fileID,'NODATA_value    -9999\n');
dlmwrite('DEM.txt',H,'-append','delimiter','\t')
fclose(fileID);

clc
dx = 2;
dy = -dx;
slopeX = 0.05;
slopeY = 0.02;
xlength = 800*cos(atan(slopeX))*2+20;
ylength = 1000*cos(atan(slopeY));
ncols = round(xlength/abs(dx));
nrows = round(ylength/abs(dy));
x11 = 0.001;
y11 = 0.001;
R = makerefmat(x11, y11+ylength, [0 dx], [dy 0]);

H = zeros(nrows,ncols);
[Rows, Cols] = meshgrid(1:nrows,1:ncols);
[X,Y] =  pix2map(R, Rows(:), Cols(:));

for i=1:numel(X)
    x = X(i); y = Y(i);
    if x <= 800*cos(atan(slopeX))
        z = (800*cos(atan(slopeX))-x)*slopeX + y*slopeY;
    elseif x <= 800*cos(atan(slopeX))+20
        z = y*slopeY;
    else
        z = (x-800*cos(atan(slopeX))-20)*slopeX + y*slopeY;
    end
    H(Rows(i),Cols(i)) = z;
end
%3D surface, move it down
figure; mapshow(H,R,'DisplayType','surface');
demcmap(H);axis normal;  axis equal; 
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
view(3); grid on
%% Generate a DEM Quiescent flow with wet–dry interface over uneven bed
% create Referece and Z value
clear,clc
dx = 50;
dy = -dx;
xlength = 8000;
ylength = 8000;
ncols = round(xlength/abs(dx));
nrows = round(ylength/abs(dy));
x11 = 0;
y11 = 0;
R = makerefmat(x11, y11+ylength, [0 dx], [dy 0]);
Z = zeros(nrows,ncols);
[Rows, Cols] = meshgrid(1:nrows,1:ncols);
[X,Y] =  pix2map(R, Rows(:), Cols(:));
% calculate Z value
z_bB1 = 2000 - 0.00032*((X-3000).^2 + (Y-5000).^2);
z_bB2 = 900 - 0.000144*((X-5000).^2 + (Y-3000).^2);
z_b = max([z_bB1*0,z_bB1,z_bB2],[],2);
Z(:) = z_b;
Ele0 = 1000; % initial elevation of water
h0 = Ele0-Z; 
h0(h0<0)=0; % initial water depth
ele0 = h0 + Z;
% figure; mapshow(Z,R,'DisplayType','surface'); axis square
%% write txt file for bed Z value
ncols = 7426;
nrows = 7726;
x11 = 423132.83288574;
y11 = 2880051.6378784;
dx = 2;
Z(isnan(Z)) = -9999;
filename = 'FuzhouDEM2m.txt';
fileID = fopen(filename,'w');
fprintf(fileID,'ncols    %d\n', ncols);
fprintf(fileID,'nrows    %d\n', nrows);
fprintf(fileID,'xllcorner    %f\n', x11);
fprintf(fileID,'yllcorner    %f\n', y11);
fprintf(fileID,'cellsize    %f\n', abs(dx));
fprintf(fileID,'NODATA_value    -9999\n');
dlmwrite(filename,Z,'-append','delimiter','\t')
fclose(fileID);
%% write txt file for bed initial value
filename = 'ele.txt';
fileID = fopen(filename,'w');
fprintf(fileID,'ncols    %d\n', ncols);
fprintf(fileID,'nrows    %d\n', nrows);
fprintf(fileID,'xllcorner    %f\n', x11);
fprintf(fileID,'yllcorner    %f\n', y11);
fprintf(fileID,'cellsize    %f\n', abs(dx));
fprintf(fileID,'NODATA_value    -9999\n');
dlmwrite(filename,ele0,'-append','delimiter','\t')
fclose(fileID);
%%
clear,clc
[Z,R] = arcgridread('DEM.txt');
%3D surface, move it down
figure; mapshow(Z,R,'DisplayType','surface');
demcmap(Z);axis normal;  axis equal; 
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
view(3); grid on;

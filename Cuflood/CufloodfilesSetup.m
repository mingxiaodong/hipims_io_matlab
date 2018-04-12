%create and write the txt files for running Cuflood model
% renewed on 11-Mar-2016
clear,clc
file_path = 'K:\CufloodCase\';
cd(file_path)
%% 1 DEM.txt was output from ArcGIS
% set the resolution and scope of running area
% read DEM data
filename = [file_path 'DEM.txt'];
fileID = fopen(filename,'r');
[Z,R] = arcgridread(filename);% R colume1: [0,dx,x11], colume2: [dy,0,y11]
% plot DEM map
mapshow(Z,R,'DisplayType','Surface')

%% ****initial water depth
% [M,I] = min(Z(:));
hValue = arcgridread('h0.txt');
filename = [file_path 'h.txt'];
Arcgridwrite(filename,hValue,R)

%% ****velocity in X direction
filename = [file_path 'u.txt'];
fileID = fopen(filename,'w');
for i=1:6
    fprintf(fileID,'%s\n',txthead{i});
end% write the head
hValue = zeros(nrows, ncols);
hValue(I)=0.001;
%value(index) = ux; % ux
dlmwrite(filename,hValue,'delimiter',' ','-append');
fclose(fileID);

%****velocity in Y direction
filename = [file_path 'v.txt'];
fileID = fopen(filename,'w');
for i=1:6
    fprintf(fileID,'%s\n',txthead{i});
end% write the head 
hValue = zeros(nrows, ncols);
%value(index) = uy; %uy
hValue(I)=0.001;
dlmwrite(filename,hValue,'delimiter',' ','-append');
fclose(fileID);

%% 3.mask
% mask.txt stores integer representing the grid types
% 0: non-boundary grid; 
% 1: rigid boundary grid; 
% 2: invalid grid; 
% <0: open boundary grid, the absolute value represents the serial number 
%of the boundary 
%**********find the surrounding boundary cell
Z1 = Z;
for i=1:nrows
    if i==1||i==nrows
        Z1(i,:) = Z(i,:)+nan(1,ncols);
    else
        Z1(i,:) = Z(i-1,:)+Z(i+1,:);
    end
end
Z2 = Z;
for j=1:ncols
    if j==1||j==ncols
        Z2(:,j) = Z(:,j)+nan(nrows,1);
    else
        Z2(:,j) = Z(:,j-1)+Z(:,j+1);
    end
end
Z3 = Z1+Z2;
clear Z1 Z2
Z3(~isnan(Z3))=0; % non-boundary cell
Z3(isnan(Z)) = nan_value; % the cells out of the research area
Z3(isnan(Z3)) = 1; %rigid boundary cell
Z3(Z3==nan_value) = 2;
%Z3(Z3==1) = -1;
%% **********find the open boundary
% -1: input boundary; -2: output boundary
x = [516276; 516352; 592700; 592990];%[516201; 516521; 553525; 554591]; %
y = [171983; 171906; 183950; 174900];%[167328; 167267; 179114; 178362]; %
[row,col] = map2pix(R,x,y);
row = ceil(row);
col = ceil(col);
[rect1X,rect1Y] = meshgrid( min(row(1:2)):max(row(1:2)), min(col(1:2)):max(col(1:2)));
[rect2X,rect2Y] = meshgrid( min(row(3:4)):max(row(3:4)), min(col(3:4)):max(col(3:4)));
ind1 = sub2ind(size(Z3),rect1X(:),rect1Y(:));
ind2 = sub2ind(size(Z3),rect2X(:),rect2Y(:));
Z4 = Z3;
Z3_s1 = Z3(ind1);
Z4(ind1(Z3_s1==1))=-1; %input boundary
Z3_s2 = Z3(ind2);
Z4(ind2(Z3_s2==1))=-2; %output boundary
%figure;mapshow(Z4,R,'DisplayType','surface');
%%
filename = [file_path 'mask.txt'];
fileID = fopen(filename,'w');
%**********write the head
for i=1:6
    fprintf(fileID,'%s\n',txthead{i});
end
dlmwrite(filename,Z3,'delimiter',' ','-append');
fclose(fileID);
%% 4. inflow   N.dat
%TideSimu = xlsread('FlowTide.xlsx','Tide','C2:C146');
%FlowSimu = xlsread('FlowTide.xlsx','Flow','C2:C146');
%simulated from tidemodel or somewhere else
Discharge = xlsread(['G:\cudaSWEsSolver\ThamesValley20m\' 'ThamesValley20m_InitialData.xlsx'],'Discharge','A:B'); 
Depth = xlsread(['G:\cudaSWEsSolver\ThamesValley20m\' 'ThamesValley20m_InitialData.xlsx'],'Depth','A:B');
h = Depth(:,2);%TideSimu(1:length(Ts)); %water depth
theta = 0;
vq = Discharge(:,2)*cos(theta);%zeros(size(Ts)); %velocity(u) or discharge(qx) in X direction
uq = Discharge(:,2)*sin(theta);%FlowSimu(1:length(Ts)); %velocity(v) or discharge(qy) in Y direction
Ts_h = Depth(:,1); Ts_f = Discharge(:,1);
%% *decide which inflow file to generate 

% 1/0: fixed/or not [h,u,v,qx,qy]
% the direction of the bound: North(1), East(2), South(3), West(4)
%inflow 1, the input discharge is not fixed as 0
%inflow 2, the tide height is not fixed as 0
hValue = [Ts_f zeros([length(Ts_f) 3])]; N = 1; d = 3; fixed = [0 0 0 1 1]; hValue(:,3) = uq; hValue(:,4) = vq;%inflow 1
%value = [Ts_h zeros([length(Ts_h) 3])]; N = 2; d = 2; fixed = [1 0 0 0 0]; value(:,2) = h; %inflow 2
%create and write the file
filename = [file_path 'inflow   ' num2str(N) '.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'%u\n',d);
fprintf(fileID,'%u %u %u %u %u\n',fixed);
fprintf(fileID,'% u\n',length(Ts_f));
fprintf(fileID,'%-7u %7.4f %7.4f %7.4f\n',hValue');
fclose(fileID);
clc;type(filename)
%% 5.source term: psc.txt
filename = [file_path 'psc.txt'];
fileID = fopen(filename,'w');
% write the head
for i=1:6
    fprintf(fileID,'%s\n',txthead{i});
end
hValue = zeros(nrows, ncols);
%value(~isnan(Z)) = 1;
dlmwrite(filename,hValue,'delimiter',' ','-append');
fclose(fileID);
%% psc   N.dat
N = 1;
Rain = xlsread(['G:\cudaSWEsSolver\LondonFlood\' 'InitialData.xlsx'],'Rainfall','A:B'); % precipitation m/s
Rain = Rain*0;
filename = [file_path 'psc   ' num2str(N) '.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'%u\n',length(Rain));
fprintf(fileID,'%-7u %15.13f\n',Rain');
fclose(fileID);
%% 6.manning factor
filename = [file_path 'Manning.txt'];
fileID = fopen(filename,'w');
% write the head
for i=1:6
    fprintf(fileID,'%s\n',txthead{i});
end
manning_factor = 0.015;
hValue = repmat(manning_factor, [nrows, ncols]);
dlmwrite(filename,hValue,'delimiter',' ','-append');
fclose(fileID);
%% 7. set_up.ste
    %*****set up the parameters******
order = 2; % 1 or 2, the order of the computational accuracy
tall = 3600*100;%Ts(end); %default value is Ts(end), means the total physical time(s)
initial_dt = 0.001; % initial time step
CFL = 0.5; %Courant number
output_inter = 3600; %output interval(s)
% type of bound in four direction1: open bound; 2: rigrid bound
bound_E = 2; bound_W = 2; bound_N = 2; bound_S = 2;
bound_body = 2; %the boundary type of invalid grid
nobr = 2; %number of input bounds
nops = 0; %number of point source batches
output_history = 'true'; %whether or not output the observations in specific points
all_cells_are_monitored = 'false';
console_in_one_line = 'true';
    %******write and save the file
filename = [file_path 'set_up.ste'];
fileID = fopen(filename,'w');
fprintf(fileID,'%s\n','order'); fprintf(fileID,'%u\n',order);           %1
fprintf(fileID,'%s\n','number of rows'); fprintf(fileID,'%u\n',nrows);  %2
fprintf(fileID,'%s\n','number of columns'); fprintf(fileID,'%u\n',ncols);%3
fprintf(fileID,'%s\n','tall'); fprintf(fileID,'%.1f\n',tall);             %4
fprintf(fileID,'%s\n','initial dt'); fprintf(fileID,'%.3f\n',initial_dt);%5
fprintf(fileID,'%s\n','CFL'); fprintf(fileID,'%f\n',CFL);             %6
fprintf(fileID,'%s\n','x length'); fprintf(fileID,'%.1f\n',ncols*dx);     %7
fprintf(fileID,'%s\n','y length'); fprintf(fileID,'%.1f\n',nrows*dy);     %8
fprintf(fileID,'%s\n','x cell size'); fprintf(fileID,'%u\n',dx);        %9
fprintf(fileID,'%s\n','y cell size'); fprintf(fileID,'%u\n',dy);        %10
fprintf(fileID,'%s\n','output interval'); fprintf(fileID,'%.1f\n',output_inter);%11
fprintf(fileID,'%s\n','bound_E'); fprintf(fileID,'%u\n',bound_E);       %12
fprintf(fileID,'%s\n','bound_W'); fprintf(fileID,'%u\n',bound_W);       %13
fprintf(fileID,'%s\n','bound_N'); fprintf(fileID,'%u\n',bound_N);       %14
fprintf(fileID,'%s\n','bound_S'); fprintf(fileID,'%u\n',bound_S);       %15
fprintf(fileID,'%s\n','bound_body'); fprintf(fileID,'%u\n',bound_body); %16
fprintf(fileID,'%s\n','nobr #number of input bounds'); fprintf(fileID,'%u\n',nobr);             %17
fprintf(fileID,'%s\n','nops #number of point sources'); fprintf(fileID,'%u\n',nops);             %18
fprintf(fileID,'%s\n','output_history'); fprintf(fileID,'%s\n',output_history); %19
fprintf(fileID,'%s\n','all_cells_are_monitored'); fprintf(fileID,'%s\n',all_cells_are_monitored); %20
fprintf(fileID,'%s\n','console_in_one_line'); fprintf(fileID,'%s\n',console_in_one_line); %21
fclose(fileID);
clc
type(filename)
%% 9. thist_gauge.dat
%coordinate of gauge
gaugeXY = xlsread(['G:\cudaSWEsSolver\ThamesValley20m\' 'ThamesValley20m_InitialData.xlsx'],'Gauge','B:C');
nog = size(gaugeXY,1); %number of gauges
%create and write the file
filename = [file_path 'thist_gauge.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'%u\n',nog);
fprintf(fileID,'%.1f %.1f\n',gaugeXY(1:end-1,:)');
fprintf(fileID,'%.1f %.1f',gaugeXY(end,:)');
fclose(fileID);
type(filename)
%% 2.initial water depth, water velocity in X and Y directions
    %%2.1 load a result file and make it as the initial conditions
    filename = [file_path 'ResultAsInitial_72h.txt'];
    result = load(filename);
    x = result(:,1); y = result(:,2); % XY coordinates
    h = result(:,4); %water depth
    ux = result(:,5); %velocity in X direction
    uy = result(:,6); %velocity in Y direction
    u = sqrt(ux.^2 + uy.^2); %water velocity
    % transfer XY coordinates to the index of matrix Z
    [row,col] = map2pix(R,x,y); row = round(row); col = round(col);
    index = sub2ind(size(Z),row,col); 
    Z_DEM = Z; Z_DEM(index) = h;
    Z_Result = nan(size(Z)); Z_Result(index) = h;    
    %plot DEM
    figure;mapshow(Z_DEM,R,'DisplayType','surface'); 
    xlabel('x (easting in meters)'); ylabel('y (northing in meters)');
    demcmap(Z); freezeColors; hold on
    % plot the data from a result file    
    mapshow(Z_Result,R,'DisplayType','surface'); title('Initial Condition');
    colorbar; caxis([0 max(Z_Result(:))]); colormap parula; hold off
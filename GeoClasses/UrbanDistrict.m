clear,clc
curFolder = 'F:\cudaSWEsSolver\UrbanDistrict\'; %current folder
cd(curFolder)
output_path = [curFolder '\output\'];
FlowType = 'Low';%'Medium', 'High'
LabData = xlsread([curFolder 'Original_staggered.xls'],['OS_' FlowType]);
gaugeValue = load([curFolder, '\output\', 'h_gauges.dat']);
TimeHis = gaugeValue(:,1);
simu_h = gaugeValue(:,2:end); %y_unit = 'm';
simu_h = simu_h*100; y_unit='cm';
lab_h = LabData(:,3:end);
%% plot gauge data
figure
for i = 1:10
    %subplot(3,3,i-1)        
    plot(TimeHis,simu_h(:,i),LabData(:,1),lab_h(:,i),'or','LineWidth',1.2,'MarkerSize',3)   
    ylim([-0.1,11]); axis square;
    legend1 = legend('Simulated','Experimental');
    set(legend1,'FontSize',10,'Position',[0.57 0.78 0.24 0.075]); %[left bottom width height]
    xlabel('Time(s)','FontSize',10,'position',[55 0.5])
    ylabel(['Depth(' y_unit ')'],'FontSize',10,'position',[6 10.3],'rotation',0)
    filename = ['Gauge No.' num2str(i)];
    title(filename,'position',[30 10.3])
    pause(2)
    %saveas(gcf,[file_path 'Gauge No_' num2str(i)],'epsc')
end
%% plot map
%%***load data
T = 1*10;
filename_h = [output_path 'h_' num2str(T) '.dat'];
filename_hU = [output_path 'hU_' num2str(T) '.dat'];
result_h = dlmread(filename_h);
result_hU = dlmread(filename_hU);
[Z_DEM,R] = arcgridread([curFolder '\input\mesh\dem.txt']);% R colume1: [0,dx,x11], colume2: [dy,0,y11]
nrows = size(Z_DEM,1); ncols = size(Z_DEM,2);
x = result_h(:,1); y = result_h(:,2);
h = result_h(:,3);
u = result_hU(:,3)./(h + 0.000001); v = result_hU(:,4)./(h + 0.000001);
vel = sqrt(u.*u + v.*v);

%transfer XY coordinates to index according to matrix Z
[row,col] = map2pix(R,x,y); 
row = ceil(row); col = ceil(col);
row(row==0) = 1; col(col==0) = 1;
row(row>nrows)= nrows; col(col>ncols) = ncols;
index = round(sub2ind(size(Z_DEM),row,col));
% Data for plotting
h_grid = nan(size(Z_DEM)); 
u_grid = nan(size(Z_DEM)); v_grid = nan(size(Z_DEM)); vel_grid = nan(size(Z_DEM));
h_grid(index) = h;
vel_grid(index) = vel;
u_grid(index) = u; 
v_grid(index) = v;

%% ***plot result map
addpath 'C:\Users\b4042552\freezeColors'
Mapping_Result = h_grid; MapTitle = 'Depth';
%Mapping_Result = vel_grid; MapTitle = 'Velocity';
figure;
%mapshow(Z_DEM,R,'DisplayType','surface');demcmap(Z_DEM); freezeColors;
xlabel('m towards East'); ylabel('m towards North');
hold on
mapshow(Mapping_Result,R,'DisplayType','surface')
%quiver(x,y,u,v,3,'r')
axis equal; colormap(parula); colorbar; 
axis([0 7 2 9])
caxis([min(Mapping_Result(:)) max(Mapping_Result(:))]); 
title([MapTitle ': (T = ' num2str(T) 's)']);
hold off
%%
X = min(x):0.050:max(x); Y = min(y):0.050:max(y);
[XX,YY] = meshgrid(X,Y);
Z_u = griddata(x,y,u,XX(:),YY(:));
Z_v = griddata(x,y,v,XX(:),YY(:));
hold on; quiver(XX(:),YY(:),Z_u,Z_v,2,'r'); hold off
clear,clc
curF = 'K:\CufloodCase\';
cd(curF)
%addpath '/Applications/MATLAB_R2014b.app/toolbox/freezeColors'
%addpath 'C:\Users\b4042552\freezeColors'

%% ****to find the newest output and get its time 
[Z,R] = arcgridread([curF '\dem.txt']);% R colume1: [0,dx,x11], colume2: [dy,0,y11]
cd([curF 'output'])
list = dir;% list of the files in current folder
%delete the non-result files
non_re = [0 0]; n = 1;
for i = 1:length(list)
    if list(i).name(1)~='r'
       non_re(n)=i;
       n=n+1;
    end
end
list(non_re)=[];
    %sort the files based on their created time
[~,index] = sortrows({list.datenum}.'); list = list(index); clear index non_re n i
curentT = list(end).name(8:end-5); % the name of the newest output
T = str2double(curentT); % the physical time of the newest file
disp(['Simulation time: ' num2str(T/3600) ' h'])
% running time 
t_run = datetime(datestr(list(end).datenum))-datetime(datestr(list(1).datenum));
disp('Elapse time:');disp(t_run)

%%*******read DEM data

%% ****load data from output file
T = 10000; %time of your file
filename = ['result ' num2str(T) 's.txt'];%[pathstr,name,ext] = fileparts(filename)
result = load(filename); %six column: X,Y,zb,h,ux,uy
x = result(:,1); %X coordinate
y = result(:,2); %Y coordinate
zb = result(:,3); %elevation of bed
h = result(:,4); %water depth
ux = result(:,5); %velocity in X direction
uy = result(:,6); %velocity in Y direction
s = zb + h; %elevation of water surface
u = sqrt(ux.^2 + uy.^2); %water velocity

%% transfer XY coordinates to index according to matrix Z
[row,col] = map2pix(R,x,y); row = round(row); col = round(col);
index = round(sub2ind(size(Z),row,col));
% data for mapping
Z_Result = nan(size(Z)); 
Z_Result(index) = h;
[Z,R] = arcgridread([curF 'DEM.txt']); % read DEM data
Z_DEM = Z; %Z_DEM(index)= h; Z_DEM(Z_DEM<=0) = min(Z(:));

%% ****mapping
figure
h0 = mapshow(Z_DEM,R,'DisplayType','surface'); zdatam(h0,Z_DEM-10000)
xlabel('x (m)'); ylabel('y (m)'); axis equal;
demcmap(Z_DEM);
%freezeColors
hold on;
mapshow(Z_Result,R,'DisplayType','surface'); axis equal
colorbar; caxis([min(Z_Result(:)) max(Z_Result(:))]); colormap parula
t_hour = floor(T/3600); t_min = (T/3600-floor(T/3600))*60;
h_min = [num2str(t_hour),'h ',num2str(t_min)];
title(['t = ',num2str(T),'s (', h_min, 'min)']);
Gauges = dlmread([curF 'thist_gauge.dat']);
h2 = plot(Gauges(2:end,1),Gauges(2:end,2),'r*');zdatam(h2,max(Mapping_Result(:)))
hold off

%% plot 'time history.dat'
gaugeXY = dlmread([curF 'thist_gauge.dat']);
gaugeXY(1,:) = [];
[row,col] = map2pix(R,gaugeXY(:,1),gaugeXY(:,2)); row = round(row); col = round(col);
index = round(sub2ind(size(Z),row,col));
gauge_Zb = Z(index);
%%
gaugeValue = dlmread([curF 'output\time_history.dat']);
TimeHis = gaugeValue(:,1)/3600;
gauge_eta = gaugeValue(:,2:end);%-gauge_Zb;%repmat(gauge_Zb,size(gaugeValue,2)-1);
gauge_h = gauge_eta - [repmat(gauge_Zb,1,size(gauge_eta,1))]';
%%
x_unit = 'Time (h)';
y_unit = 'Water depth (m)';
figure
PlotData = gauge_h; PlotTime = TimeHis;
legendName = cell(1,size(PlotData,2));
ShowPoints = 1:4;
for i=1:size(PlotData,2)
    legendName{i} = ['Point ' num2str(i)];
end
set(groot,'defaultAxesColorOrder',[1 0 0;0 1 0;0 0 1;0 1 1; 1 0 1],'defaultAxesLineStyleOrder','-|--|-.|:')

plot(PlotTime,PlotData(:,ShowPoints),'LineWidth',1.2)
set(groot,'defaultAxesLineStyleOrder','remove');set(groot,'defaultAxesColorOrder','remove')
legend(legendName{ShowPoints},'Location','best');xlabel(x_unit);ylabel(y_unit)
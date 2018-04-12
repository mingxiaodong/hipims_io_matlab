%% read files
clear,clc
% [z_dem, r_dem] = arcgridread('G:\cudaSWEsSolver\London38001\MultiGPU\DEM.txt');
catchmentPoly = shaperead('G:\Data\London\38001.shp');
riverPoly = shaperead('G:\Data\London\WaterSystemLondon.shp');
gaugeCoors = dlmread('G:\cudaSWEsSolver\London38001\input\field\gauges_pos.dat');
%% draw map
figure
hold on
fig_h1 = mapshow(z_dem, r_dem, 'DisplayType','Surface'); 
zdatam(fig_h1,z_dem-100)
demcmap(z_dem)
% linkaxes([ax1 ax2])
% mapshow(z_dem, r_dem, 'DisplayType','contour','LineColor','black','ShowText','on')
mapshow(catchmentPoly,'Facecolor','none') % show catchment boundary
axis image 
axis manual %off
riverSpec = makesymbolspec('Line',{'grid_code',5, 'LineWidth',2}, ...
                                  {'grid_code',4, 'LineWidth',1.5},...
                                  {'grid_code',3, 'LineWidth',1},...
                                  {'Default','Color',[0 0 1],'LineWidth',0.2});
mapshow(riverPoly,'SymbolSpec',riverSpec) % show rivers
mapshow(gaugeCoors(:,1),gaugeCoors(:,2),'DisplayType','point','Color','red') % show gauges
hold off

%************* ajust the Xticks and Yticks
ax_h = gca;
gridSize = 10000; % the grid to show in the map
xlim = ax_h.XLim;
ylim = ax_h.YLim;
ax_h.XTick  = ceil(xlim(1)/gridSize)*gridSize:gridSize:xlim(2);
ax_h.YTick  = ceil(ylim(1)/gridSize)*gridSize:gridSize:ylim(2);
ax_h.XTickLabel = ax_h.XTick/1000;
ax_h.YTickLabel = ax_h.YTick/1000;
ax_h.YTickLabelRotation = 90;
xlabel('km towards East');
ylabel('km towards North');
grid on
box on
%**********************************************

%% plot output of London
clear,clc
% read data
hU_gauges = dlmread('G:\cudaSWEsSolver\London38001\output\hU_gauges.dat');
load('G:\Data\London\LeeFeildesWeir15minFlow.mat')
%% define variables
hUy = hU_gauges(:,3:2:end);
hUy_1 = hUy(:,1:10);
hUy_2 = hUy(:,11:end);
hUy_sum1 = sum(hUy_1,2);
hUy_sum2 = sum(hUy_2,2);
t0 = datetime(2014,2,5,0,0,0);% starting date and time of the event simulated
eventDur = 24*4; %
ind = find(LeeFeildesWeir15minFlow.Timestamp==t0);

simu_Value = hUy_sum2*(-20);
simu_DateTime = hU_gauges(:,1)/3600/24+t0;

obv_DateTime = LeeFeildesWeir15minFlow.Timestamp(ind:ind+eventDur*4);
obv_Value = LeeFeildesWeir15minFlow.Discharge(ind:ind+eventDur*4);
%% plot lines for simulated and observed values
figure; plot(obv_DateTime,obv_Value,simu_DateTime,simu_Value)
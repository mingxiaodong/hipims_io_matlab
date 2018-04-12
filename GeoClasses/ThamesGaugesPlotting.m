clear,clc
addpath C:\Users\b4042552\Dropbox\Matlab\myfunction
cd('G:\cudaSWEsSolver\LondonFrame4\output')
load C:\Users\b4042552\Dropbox\Matlab\London\multiEventData.mat
load 'C:\Users\b4042552\Google Drive\MyResearch\London\matlab\ThamesGaugesTide2014_0331_1202.mat'
gaugeCoor = dlmread('G:\cudaSWEsSolver\LondonFrame4\gauges_pos.dat');
t0 = multiEventData(2).t0;
timeObv_all = ThamesGaugesTide20140331.record_datetime;
dataObv_all = table2array(ThamesGaugesTide20140331(:,2:end));
%%load output data
eta_gauged = dlmread('eta_gauges.dat');
timeStep = eta_gauged(:,1);
eta_simu = eta_gauged(:,2:9);
%%plot observation and simulated data
indObv = timeObv_all>=t0&timeObv_all<=t0+max(timeStep)/24/3600;
timeObv = timeObv_all(indObv); %timeObv = 3600*24*datenum(timeObv-t0);
dataObv = dataObv_all(indObv,:); 
timeSimu = timeStep; timeSimu = timeSimu/3600/24+t0;
gaugeNames = ThamesGaugesTide20140331.Properties.VariableNames;
gaugeNames = gaugeNames(2:9);
figure
for i = 1:8
    subplot(4,2,i)
    plot_h = plot(timeObv,dataObv(:,i),timeSimu,eta_simu(:,i));
    ylim([-5 5])
    title(gaugeNames{i})
end
legend_hl = legend('Observed','Simulated');
set(legend_hl,'Position',[0.485 0.485 0.06 0.05]);

%%
fileList = dir('*.asc');
[Z,R] = arcgridread(fileList(2).name);
Z(Z<=0)=nan;
figure
dem_h = mapshow(Z,R,'DisplayType','Surface');
zdatam(dem_h,Z-600)
mapshow(gaugeCoor(:,1),gaugeCoor(:,2),'DisplayType','Point')
textStr = cellstr(num2str((1:length(gaugeCoor))'));
text(gaugeCoor(:,1),gaugeCoor(:,2),textStr) 
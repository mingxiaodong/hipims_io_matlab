%% plot water level or water depth at gauges
%% import data
clear,clc
CaseName = 'Haltwhistle5m';
% cd(['G:\cudaSWEsSolver\' CaseName '\output']);
gaugeValue = dlmread('h_gauges.dat');
% cd('G:\Cuflood10m\output'); gaugeValue = dlmread('time_history.dat');
gaugeValue(end,:)=[];
TimeHis = gaugeValue(:,1); 
gauge_h = gaugeValue(:,2:end);
x_unit = 'h'; y_unit = 'm';
[num,txt,raw]= xlsread('G:\cudaSWEsSolver\ThamesValley20m\ObservedData20140331.xlsx');
gaugeName = txt(1,3:end);
%% compare the observasion and simulation at each gauge
z = [6.5, 13.5, 3.8, 3, 5.5, 5.8, 4, 1.5]*0; % gauge depth adjustment
T_range = [0 inf]; %time range(h) for plotting
ObservedTime = 0:900:max(TimeHis); %second
SimulatedData = [TimeHis/3600, gauge_h];
ObservedData = [ObservedTime'/3600, num(1:length(ObservedTime),:)];
%delete the data out of the time scale
SimulatedData(SimulatedData(:,1)<T_range(1)|SimulatedData(:,1)>T_range(2),:)=[];
ObservedData(ObservedData(:,1)<min(SimulatedData(:,1))|ObservedData(:,1)>max(SimulatedData(:,1)),:)=[];

figure
for i = 1:8
    subplot(4,2,i)
    %plot(SimulTime,gauge_h(:,i)-z(i),ObservedData.record_time(1:length(ObTime)),ObservedData{1:length(ObTime),i+1},'LineWidth',1.2)
    plot(SimulatedData(:,1),SimulatedData(:,i+1)-z(i),ObservedData(:,1),ObservedData(:,i+1),'LineWidth',1.2)
    title(gaugeName{i})
    xlabel(['Time (' x_unit ')']); 
    ylabel(['Water level (' y_unit ')'])
end
legend1 = legend('Simulated','Observed');
set(legend1,'Position',[0.49 0.5 0.05 0.05]);
%% plot time step
timestep = dlmread('timestep_log.txt');
T_range = [3600*0.5 inf];
ind = timestep(:,1)>T_range(1)&timestep(:,1)<T_range(2);
plotData = timestep(ind,:);
figure; plot(plotData(:,1),plotData(:,2),'-')
xlabel('t'); ylabel('dt');
disp(sum(plotData(:,2)<0.1)/length(plotData))
%% show all the simulation of all gauges in one plot
figure
ShowPoints = 1:10;%size(gauge_h,2);
legendName = cell(1,size(gauge_h,2));
for i=1:size(gauge_h,2)
    legendName{i} = ['Gauge ' num2str(i)];
end
set(groot,'defaultAxesColorOrder',[1 0 0; 0 0.8 0; 0 0 1; 0.8 0.8 0],...
    'defaultAxesLineStyleOrder','-|--|-.|:')
plot(TimeHis/3600,gauge_h(:,ShowPoints),'LineWidth',1.2)
set(groot,'defaultAxesLineStyleOrder','remove');set(groot,'defaultAxesColorOrder','remove')
legend(legendName{ShowPoints},'Location','best');
xlabel(['Time (' x_unit ')']); ylabel(['Water depth (' y_unit ')'])
%% show the simulation at each gauge with multiple plots
figure
ShowPoints = 1:10;
for i = 1:10
    subplot(5,2,i)
    plot(TimeHis/3600,gauge_h(:,ShowPoints(i)),'LineWidth',1.2)
    title(legendName{ShowPoints(i)})
    xlabel(['Time (' x_unit ')']); 
    ylabel(['Depth (' y_unit ')'])
end
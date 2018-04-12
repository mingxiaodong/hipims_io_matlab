clc
cd('C:\Users\b4042552\Dropbox\Matlab\OutputPlotting')
gaugeValue1 = dlmread('h_gauges_n.dat');
gaugeValue2 = dlmread('time_history_n.dat');
TimeHis1 = gaugeValue1(:,1)/3600; 
TimeHis2 = gaugeValue2(:,1)/3600;
gauge_h1 = gaugeValue1(:,2:end);
gauge_h2 = gaugeValue2(:,2:end);
x_unit = 'Time (h)';
y_unit = 'Water depth (m)';
%% compare the observations in the same gauge
figure
ShowPoints = 11:20;
for i = 1:length(ShowPoints)
    subplot(2,5,i)
    plot(TimeHis1,gauge_h1(:,ShowPoints(i)),TimeHis2,gauge_h2(:,ShowPoints(i)),'LineWidth',1.5)
    title(['gauge' num2str(ShowPoints(i))])
    legend('OverlandFlow','Cuflood','Location','best')
    xlabel(x_unit);ylabel(y_unit)
end
%% compare the observations between gauges
figure
PlotData = gauge_h1; PlotTime = TimeHis1;
%PlotData = gauge_h2; PlotTime = TimeHis2;
legendName = cell(1,size(PlotData,2));
ShowPoints = 11:20;
for i=1:size(PlotData,2)
    legendName{i} = ['Point ' num2str(i)];
end
set(groot,'defaultAxesColorOrder',[1 0 0;0 1 0;0 0 1;0 1 1; 1 0 1],'defaultAxesLineStyleOrder','-|--|-.|:')

plot(PlotTime,PlotData(:,ShowPoints),'LineWidth',1.2)
set(groot,'defaultAxesLineStyleOrder','remove');set(groot,'defaultAxesColorOrder','remove')
legend(legendName{ShowPoints},'Location','best');xlabel(x_unit);ylabel(y_unit)
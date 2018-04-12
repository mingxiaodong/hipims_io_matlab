% plot cuda model results based on saved data  
clear,clc
curFolder = 'G:\cudaSWEsSolver\Fuzhou';
cd(curFolder)
addpath 'C:/Users\b4042552/Documents/MATLAB/freezeColors'
load('Result_H.mat')
% load('Result_24h.mat', 'TimeSeries', 'R', 'Z_DEM', 'gauge_t_h', 'gauge_t_hU')
%% plot gauge value
% load('Result_24h.mat', 'gauge_t_h', 'gauge_t_hU')
gaugeValue = gauge_t_h;
TimeHis = gaugeValue(:,1)/3600; x_unit = 'h';
gauge_h = gaugeValue(:,2:end); y_unit = 'm'; y_name = 'Water depth';
ShowPoints = 1:size(gauge_h,2);
legendName = cell(1,size(gauge_h,2));
for i=1:size(gauge_h,2)
    legendName{i} = ['Gauge ' num2str(i)];
end; clear i;

figure
line_h = plot(TimeHis,gauge_h(:,ShowPoints),'LineWidth',1.2); 
legend(legendName{ShowPoints},'Location','best'); clear legendName ShowPoints
xlabel(['Time (' x_unit ')']); ylabel([y_name ' (' y_unit ')']); 
clear x_unit y_unit y_name
%% plot water depth map
t0 = datetime('2016-9-14 8:00','InputFormat','yyyy-MM-dd HH:mm','Format','yyyy-MM-dd HH:mm');
figure;
tic
for i = 55:56%:length(TimeSeries);
    T = TimeSeries(i);
    TimeString = ['Time = ' datestr(t0+T/3600/24,'mm/dd HH:MM am')];
    Mapping_Grid = ResultCell_H{i,1}; MapTitle = 'Water Depth'; Mapping_Grid(Mapping_Grid<=0.01)=nan;    
    h0 = mapshow(Z_DEM,R,'DisplayType','surface'); demcmap(Z_DEM); freezeColors; zdatam(h0,Z_DEM-600)
    hold on
    h1 = mapshow(Mapping_Grid,R,'DisplayType','surface');
    colormap(parula); colorbar;
    caxis([0 2]);
    text_hd = text(R(3,1),R(3,2),TimeString);
    xlabel('Meter towards East'); ylabel('Meter towards North'); title(MapTitle)
    axis([4.14e5 4.45e5 2.868e5 2.9e5])
    axis equal
    hold off
    PictureName = [curFolder '\pictures72\' num2str(i-1,'%02u')];
    print(PictureName,'-djpeg','-r300')
end
toc
%% make animation
cd([curFolder '\pictures72'])
AnimationName = [curFolder '\animation72.gif'];
list = dir('*.jpg');
N = size(list);
for i=1:N
    list(i).datenum = str2double(list(i).name(1:end-4));
end
[~,sortID] = sortrows([list.datenum].'); list = list(sortID); clear sortID
for i = 1:N
    img = imread(list(i).name);
    imshow(img);
    frame=getframe(gcf);
    im=frame2im(frame);
    [A,map]=rgb2ind(img,256);
    if i == 1
		imwrite(A,map,AnimationName,'gif','LoopCount',Inf,'DelayTime',0.5);
	else
		imwrite(A,map,AnimationName,'gif','WriteMode','append','DelayTime',0.5);
    end
end
%% show the largest value
[M,I] = max(vel);
% I = vel(:)>30;
hold on; plot(x(I),y(I),'*r'); hold off; %zdatam(h1,Mapping_Result-50)
%caxis([0 50]);
%%  load result data
clear,clc
OutputFolder = '/Users/b4042552/Google Drive/Data?10m/';
list = dir(OutputFolder);
ind = false(1,length(list));
for i = 1:length(list)
    if strcmp(list(i).name(end),'c')
        ind(i)=true;
    end
end
list1 = list(ind);
%% plot water depth map
% define the start time and date
t0 = datetime('2016-9-14 19:00','InputFormat','yyyy-MM-dd HH:mm','Format','yyyy-MM-dd HH:mm');
curFolder = cd;
PictureStoreFolder = [curFolder ''];
figure;
tic
for i = 1:length(TimeSeries);
    T = TimeSeries(i);
    TimeString = ['Time = ' datestr(t0+T/3600/24,'mm/dd HH:MM am')];
    % Mapping_Grid = ResultCell_H{i,1};
    MapTitle = 'Water Depth'; Mapping_Grid(Mapping_Grid<=0.01)=nan;    
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
    PictureName = num2str(i-1,'%02u');
    print([PictureStoreFolder PictureName],'-djpeg','-r300')
end
toc
%% make animation
cd(PictureStoreFolder)
AnimationName = [PictureStoreFolder '?animation.gif'];
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
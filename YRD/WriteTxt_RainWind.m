%% WriteTxtfiles for rain and wind in YRD
clc,clear
load('CountyEvents255010.mat')
load('CountyInfor.mat', 'CountyInfor')
RW_Events = table;
n = 1;
for i=1:140
    Rain = [CountyEvent{i}.MaxRain]';
    Wind = [CountyEvent{i}.MaxWind]';
    Name = repmat({CountyInfor(i).ename},size(Rain));
    Code = repmat(CountyInfor(i).code,size(Rain));
    NoC = repmat(i,size(Rain));
    RW_T = table(NoC,Code,Name,Rain,Wind);    
    RW_Events(n:n+length(Rain)-1,:) = RW_T;
    n = height(RW_Events)+1;
end
%%
writetable(RW_Events,['RW_Events' num2str(RainMin1) num2str(RainMin2) num2str(WindMin) '.txt'])
%% WriteTxtfiles for rain and wind in YRD
clc,clear
load('StationEvent.mat')
load('StationData.mat')
stationCode = [StationData.StationCode]';
stationEvents = table;
n = 1;
for i=1:numel(stationCode)
    Rain = [StationEvent{i}.MaxRain]';
    Wind = [StationEvent{i}.MaxWind]';
    Code = repmat(StationData(i).StationCode,size(Rain));
    NoC = repmat(i,size(Rain));
    RW_T = table(NoC,Code,Rain,Wind);    
    stationEvents(n:n+length(Rain)-1,:) = RW_T;
    n = height(stationEvents)+1;
end
%%
writetable(stationEvents,'StationEvent.txt')

%% Identify the rainfall with strong wind events in Yangtze River Delta (YRD), based on the daily weather station data
%18-8-2015
clc,clear

% load meteo data (from 1971 to 2011 with all missing values fixed)
load('FixedData_YRD')

% pick the meteorological data of each station
StationCode = unique(FixedData(:,1));
StationData = cell(length(StationCode),1);
for i = 1:length(StationCode)
    StationData{i} = FixedData(FixedData(:,1)==StationCode(i),:);
end


%% to store multi-hazard events based on station data
StationEvent = cell(length(StationCode),1);

% set threshold of rainfall and strong wind events
RainMin1 = 10; % mm, (medium rain): all days with rainfall larger than this will be considered within events
RainMin2 = 25; % mm, (large rain): at least one day rainfall should be larger than this
WindMin = 8;   % m/s

for i = 1:length(StationEvent)    
    UnitData = StationData{i};
    %*** select the large rain(mm>RainMin1) days
    RainDays = zeros(size(UnitData,1),1);
    RainDays((UnitData(:,6)>=RainMin1))=1;
    
    %*** transfer the one non-rain day between rain days as rain day
    % eg: 101 to 111, 1001 still as 1001
    RainDays2 = RainDays; %transfer rain days
    for j = 2:length(RainDays2)-1
        RainDays2(j) = RainDays(j)+RainDays(j-1)*RainDays(j+1);
    end
    
    %*** find the continuous rain process
    fid = find(RainDays2>0);
    RP_Start = [fid(1);0]; % record rain process start
    RP_End = [0;0];       % record rain process end
    j=1; % rainy day counter
    n=1; % rainy process counter
    for j = 1:length(fid)-1
        if fid(j)-fid(j+1)~=-1
            RP_End(n) = fid(j);
            RP_Start(n+1) = fid(j+1);
            n=n+1;
        end
    end
    RP_End(n) = fid(end);
    
    %*** identify multi-hazard events
    %%Identify multi-hazard events
    Event = struct('DateStart',0,'DateEnd',0,'WholeRec',zeros(2,10),...
    'TotalRain',0,'MaxRain',0,'MaxWind',0,'Duration',0);
    m=1; %multi-hazard event counter
    for n=1:length(RP_Start)
        WholeRec = UnitData(RP_Start(n):RP_End(n),1:7);
        MaxRain = max(WholeRec(:,6));
        MaxWind = max(WholeRec(:,7));
        % judge whether it is a multi-hazard event according to the max rain of
        % rainy days and max wind
        if MaxRain>=RainMin2 && MaxWind>=WindMin
            Event(m).DateStart = UnitData(RP_Start(n),5);
            Event(m).DateEnd = UnitData(RP_End(n),5);
            Event(m).Duration = Event(m).DateEnd-Event(m).DateStart+1;
            Event(m).WholeRec = WholeRec;
            Event(m).TotalRain = sum(WholeRec(:,6));
            Event(m).MaxRain = MaxRain;
            Event(m).MaxWind = MaxWind;
            m = m+1;
        end
    end    
    StationEvent{i} = Event;
end
%save StationEvent StationEvent


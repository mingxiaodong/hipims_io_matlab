%% 20-8-2015 ifentify multi-hazard events
%based on Daily rain and wind of each county from 01-Jan-1972 to 31-Dec-2011
clear,clc
%load for vars (DailyCountyRain, DailyCountyWind, DateEnd, DateStart)
load('DailyCountyHaz')

%%set threshold of rainfall and strong wind events
RainMin1 = 15; % mm, (large rain): all days with rainfall larger than this will be considered within events
RainMin2 = 30; % mm, (huge rain): at least one day rainfall should be larger than this
WindMin = 10;  % m/s
DateRange = (DateStart:DateEnd)';
[Year,Month,Day] = datevec(DateRange);
DateRange = [DateRange, Year, Month, Day];

CountyEvent = cell(140,1);
for i = 1:140    
    DailyRain = DailyCountyRain(i,:);
    DailyWind = DailyCountyWind(i,:);
    %*** select the large rain(mm>RainMin1) days
    RainDays = zeros(1,length(DailyRain));
    RainDays((DailyRain>=RainMin1))=1;
    
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
        RainRec = DailyRain(RP_Start(n):RP_End(n))';
        WindRec = DailyWind(RP_Start(n):RP_End(n))';
        MaxRain = max(RainRec);
        MaxWind = max(WindRec);
        WholeRec = [DateRange(RP_Start(n):RP_End(n),:), RainRec, WindRec];
        TotalRain = sum(RainRec);
        % judge whether it is a multi-hazard event according to the max rain of
        % rainy days and max wind
        if MaxRain>=RainMin2 && MaxWind>=WindMin
            Event(m).DateStart = DateRange(RP_Start(n),1);
            Event(m).DateEnd = DateRange(RP_End(n),1);
            Event(m).Duration = Event(m).DateEnd-Event(m).DateStart+1;
            Event(m).WholeRec = WholeRec;
            Event(m).TotalRain = TotalRain;
            Event(m).MaxRain = MaxRain;
            Event(m).MaxWind = MaxWind;
            m = m+1;
        end
    end    
    CountyEvent{i} = Event;
end
save(['CountyEvents' num2str(RainMin1) num2str(RainMin2) num2str(WindMin) '.mat'], 'CountyEvent', 'RainMin1', 'RainMin2', 'WindMin')
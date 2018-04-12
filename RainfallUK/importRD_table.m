%% similar with function importRDf_table
clear,clc
load('StationRec_RD.mat','StationRec'); %load information of daily recorded met stations in London
fileaddress = 'E:\Data\CEDA\RD\yearly_files\';
%6711: Cross Ness S WKS; 721:Kew; 695: Hampstead
GaugeID = 5590;%721;
GaugeName = 'BARNET';%'CrossNess';

%%information of the selected station
m = find([StationRec.src_id]==GaugeID); %met station: m th station in 'StationRec' 
%StationID = StationRec(m).src_id; % srd_id of the met station
RecYears = StationRec(m).RecYear; % recorded years in this Met station

%%***2015-12-7 use function importTable_RD to scan yearly orginal data files and extract records of the selected station
RD = table;
n = 1; %sequence of records in variables newly defined above
for i = 1:length(RecYears) %import data year by year
    tic
    y = RecYears(i); % years when the station has data recorded   
    % invoke subfunction
    RD1 = importTable_RD(y,fileaddress,GaugeID);             
    RD(n:n+size(RD1,1)-1,:) = RD1;
    n = n+size(RD1,1);
    disp(['Progress so far: ' num2str(y)])
    toc
end
disp(['Mission Accomplished! Met Station: "' StationRec(m).Name '"'])
RD_T = RD;
save(['RainD_' GaugeName '_FullT'], 'RD_T')

%% use function importRDtxt to scan yearly orginal data files and extract records of the selected station
for i = 102:length(RecYears) %import data year by year
    tic
    y = RecYears(i); % years when the station has data recorded   
    % invoke subfunction
    RD1 = importRDtxt(y,fileaddress,GaugeID);             
    RD(n:n+size(RD1,1)-1,:) = RD1;
    n = n+size(RD1,1);
    disp(['Progress so far: ' num2str(y)])
    toc
end
disp(['Mission Accomplished! Met Station: "' StationRec(m).Name '"'])
%Station_RD_T = RD;
%% extract daily rainfall data from RD data in table
clear,clc
GaugeName = 'Hampstead';
load(['RainD_' GaugeName '_FullT.mat'])
DateVec = (min(RD_T.OB_DATE):max(RD_T.OB_DATE))';
RainValue = zeros(length(DateVec),1);
%deal with the multiple recorded daily value
MultiRec = struct('Date',0,'Ind',0,'Recs',0,'Version',0);
n = 1;
for i = 1:length(DateVec)
    ind = find(RD_T.OB_DATE==DateVec(i));    
    if length(ind)>1
        MultiRec(n).Date = datestr(DateVec(i));
        MultiRec(n).Ind = ind;
        MultiRec(n).Recs = RD_T.PRCP_AMT(ind);
        MultiRec(n).Version = RD_T.VERSION_NUM(ind);
        RainValue(i) = MultiRec(n).Recs(end);
        n = n+1;
    elseif length(ind)==1
        RainValue(i) = RD_T.PRCP_AMT(ind);        
    else
        RainValue(i) = nan;
    end
end
DailyRain = table(DateVec,RainValue);
DailyRain.DateVec.Format = 'yyyy/MM/dd';
NanRec = DailyRain(isnan(DailyRain.RainValue),:);
save(['RainD_' GaugeName], 'DailyRain', 'MultiRec', 'NanRec')
%% calculate monthly rain based on the daily rain data
clear,clc
load('DailyRain_CrossNess.mat')
MonthSum = zeros(12,1);
MonthD_Mean = MonthSum;
MonthD_Max = MonthSum;
MonthD_HeavyR = MonthSum;
RecNum = MonthSum;
for i = 1:12
    MonthRain = RainValue(RainValue(:,2)==i,4);
    RecNum(i) = length(MonthRain) - sum(isnan(MonthRain));
    MonthSum(i) = sum(MonthRain,'omitnan');
    MonthD_Mean(i) = mean(MonthRain,'omitnan');
    MonthD_Max(i) = max(MonthRain);
    MonthD_HeavyR(i) = sum(MonthRain>25);
    
end
MonthData = table(MonthSum,MonthD_Mean,MonthD_Max,MonthD_HeavyR,RecNum);
MonthAver = MonthData.MonthD_Mean.*(eomday(1900, 1:12))';
bar(MonthAver)
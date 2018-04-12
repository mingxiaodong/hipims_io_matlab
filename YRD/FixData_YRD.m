%% Deal with the raw Txt data
clear,clc
load('RawTxtData_YRD.mat')
%   20-20????	0.1mm
%   ????	0.1m/s
%   32744      ????  ?? replaced by 0
%   32700      ????  ?? replaced by 0
%   32766      ????  ?? replaced by NaN
%   31XXX      ???    ?(???????replaced by XXX
%   30XXX      ???    ??? replaced by XXX
%   32XXX      ???    ??? replaced by XXX
%   1000+xxx   ??      ????xxx????????"1000"

%* deal with 32766 for rain and wind
RW = RawTxtData(:,5:6); %Rain and Wind Record column
RW(RW==32744) = 0;
RW(RW==32700) = 0;
RW(RW==32766) = NaN;

%*  deal with Rain
Rain = RW(:,1); %Rain column
FixV = zeros(size(Rain));
    %** deal with 32XXX
FixV1=FixV;
FixV1(Rain>=32000)=32000;
Rain = Rain-FixV1;
    %** deal with 31XXX
FixV1=FixV;
FixV1(Rain>=31000)=31000;
Rain = Rain-FixV1;
    %** deal with 30XXX
FixV1=FixV;
FixV1(Rain>=30000)=30000;
Rain = Rain-FixV1;
disp(['Rain range: [',num2str(min(Rain)),',',num2str(max(Rain)),']'])

%*  deal with wind
Wind = RW(:,2); %Wind column
FixV = zeros(size(Wind));
FixV(Wind>1000)=1000;
Wind = Wind-FixV;
disp(['Wind range: [',num2str(min(Wind)),',',num2str(max(Wind)),']'])

%* divided by 10 to match the unit mm, m/s
FixedData_NaN = [RawTxtData(:,1:4),Rain/10,Wind/10]; %only NaN values that mean '??' are not dealed m
save FixedData_NaN FixedData_NaN
%% produce weather station data struct, only NaN value is not been disposed
clear,clc
load FixedData_NaN 
load StationInfor
StationNum = unique(FixedData_NaN(:,1));
StationData = struct('Name',0,'Province',0,'StationCode',0,'GeoCoor',[0,0,0],'StartDate',[0,0,0],'EndDate',[0,0,0],'WholeRec',0,...
    'RecordDays',0);

WholeStation = cell2mat(StationInfor(:,3:6));
n = 1;% No Match station counter
ToDel = [0;0];
for i=1:length(StationNum)
    
    fid = find(WholeStation(:,1)==StationNum(i));       
    if isempty(fid)
        ToDel(n)=i;
        n=n+1;
    end
    
    StationData(i).Name = StationInfor{fid,1};
    StationData(i).Province = StationInfor{fid,2};
    StationData(i).GeoCoor = WholeStation(fid,2:4);
    StationData(i).StationCode = StationNum(i);
    StationData(i).WholeRec = FixedData_NaN(FixedData_NaN(:,1)==StationNum(i),:);
    SD = datevec(min(datenum(StationData(i).WholeRec(:,2:4))));
    ED = datevec(max(datenum(StationData(i).WholeRec(:,2:4))));
    StationData(i).StartDate = SD(1:3);
    StationData(i).EndDate = ED(1:3);
    StationData(i).RecordDays = size(StationData(i).WholeRec,1);
    
end
StationData(ToDel) = [];
save StationData StationData

%% replenish missing value 3276.6 by interpolating
clear,clc
% input raw data (column 6-7: rain-wind), missing value exists
load('sourcedata_hazard.mat', 'rawdata')
rawdata = rawdata(rawdata(:,2)>1970,:);

% find missing value
fid1 = find(rawdata(:,6)==3276.6);
fid2 = find(rawdata(:,7)==3276.6);

% interpolate the missing value of rain
for i= 1:numel(fid1)
    % Extract all stations' data on the missing date 
    MissingDate = rawdata(fid1(i),5);
    DateId = find(rawdata(:,5) == MissingDate);
    DayData = rawdata(DateId,:);
    DayData(DayData(:,6)==3276.6,:) = []; % delete the row of missing value
    % interpolate
    rawdata(fid1(i),6) = griddata(DayData(:,8),DayData(:,9),DayData(:,6),...
        rawdata(fid1(i),8),rawdata(fid1(i),9));
end

% interpolate the missing value of wind
for i= 1:numel(fid2)
    % Extract all stations' data on the missing date 
    MissingDate = rawdata(fid2(i),5);
    DateId = find(rawdata(:,5) == MissingDate);
    DayData = rawdata(DateId,:);
    DayData(DayData(:,7)==3276.6,:) = []; % delete the row of missing value
    % interpolate
    rawdata(fid2(i),7) = griddata(DayData(:,8),DayData(:,9),DayData(:,7),...
        rawdata(fid2(i),8),rawdata(fid2(i),9));
end
FixedData = rawdata;
save FixedData_YRD FixedData

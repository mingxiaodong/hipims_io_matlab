%% draw rain and wind data based on county code and date

CountyCode = 310108;
% from 1-1-1971 to 31-12-2011
StartDay = datenum(2009,4,30);
Duration = 90;
load CountyID
load( 'DailyCountyVar')
[Y,R,D] = datevec((StartDay:StartDay+Duration-1)');
n = find(CountyID==CountyCode);
Rain = DailyCountyRain(StartDay-DateStart+1:StartDay-DateStart+Duration)';
Wind = DailyCountyWind(StartDay-DateStart+1:StartDay-DateStart+Duration)';
WholeRec = [Y,R,D,Rain,Wind];

%% Draw rain and wind data based on weather station and date
StationCode = 000;
load('FixedData_NaN')
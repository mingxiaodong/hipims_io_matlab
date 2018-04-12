%% 19-8-2015 interpolate station data to county data based on the average value of km points inside each county
clear,clc
load YRD_GeoGrid % input Geo Coords of the km grid
% input YRD meteo data
load FixedData_NaN
% add staion coordinate to meteo data
load StationData
StationCode = [StationData.StationCode]; %all the recorded station in or around YRD
FixedData_NaN(:,10)= datenum(FixedData_NaN(:,2:4));

for i=1:length(StationCode)
    CodeMch = find(FixedData_NaN(:,1)==StationCode(i));
    FixedData_NaN(CodeMch,7) = StationData(i).GeoCoor(1);
    FixedData_NaN(CodeMch,8) = StationData(i).GeoCoor(2); 
    FixedData_NaN(CodeMch,9) = StationData(i).GeoCoor(3); 
end

%%screen out the points in YRD
Non0 = CkmG~=0;
X = LonG(Non0);
Y = LatG(Non0);
N = CkmG(Non0); % code of county that the grid belongs to

% store interpolated daily rain and wind in each county
DailyCountyRain = zeros(140,1);
DailyCountyWind = zeros(140,1);

%% Interpolation
tic
DateStart = datenum(1972,8,30);
DateEnd = datenum(1972,8,31);
for i = 0:DateEnd-DateStart
    DailyPoints = FixedData_NaN(FixedData_NaN(:,10)==i+DateStart,[7 8 5 6]); % rawdata in the selected day
    XYVrain = DailyPoints(~isnan(DailyPoints(:,3)),[1 2 3]); %delete the NaN value
    XYVwind = DailyPoints(~isnan(DailyPoints(:,4)),[1 2 4]);    
    % interpolation from station to grid points
    F1 = scatteredInterpolant(XYVrain(:,1), XYVrain(:,2), XYVrain(:,3),'natural','nearest'); %inner and external interpolation
    Rain = F1(X, Y);
    F2 = scatteredInterpolant(XYVwind(:,1), XYVwind(:,2), XYVwind(:,3),'natural','nearest');
    Wind = F2(X, Y);
    % calculate the mean value of each county based on the belonged grids value
    for j = 1:140
        GridRain = Rain(N==j);
        GridRain(isnan(GridRain)) = [];
        DailyCountyRain(j,i+1) = mean(GridRain);
        GridWind1 = Wind(N==j);
        GridWind = GridWind1;
        GridWind(isnan(GridWind)) = [];        
        DailyCountyWind(j,i+1) = mean(GridWind);
    end
    % show the process
    if i~=0
        % delete the last row in Command Window
        fprintf(repmat('\b',1,length(datestr(DateStart+i-1))+1))
    end
    disp(datestr(i+DateStart)) %show the current day
end
toc

%% mapping rain and wind on one day
RainG = NaN(size(CkmG));
WindG = RainG;
RainG(Non0) = Rain;
WindG(Non0) = Wind;
subplot(1,2,1)
geoshow(LatG,LonG,RainG,'DisplayType','mesh')
title(['Rainfall (mm)  ' datestr(i+DateStart)])
colorbar('location','East')
subplot(1,2,2)
geoshow(LatG,LonG,WindG,'DisplayType','mesh')
title(['Windspeed (m/s)  ' datestr(i+DateStart)])
colorbar('location','East')

%% show efficient weather station in YRD

fid = ~isnan(DailyPoints(:,4)); % efficient rain or wind station recorded
NaNPoints = NaN(size(CkmG));
NaNPoints(Non0) = Wind;
geoshow(DailyPoints(fid,2), DailyPoints(fid,1),'DisplayType','point')
hold on
geoshow(LatG,LonG,NaNPoints,'DisplayType','surface')
hold off
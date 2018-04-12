%% interpolate station data to county data based on the average value of km points inside each county
clear,clc
load YRD_GeoGrid % input Geo Coords of the km grid
load('sourcedata_hazard.mat') % input raw meteo data

% screen out the points in YRD
Non0 = CkmG~=0;
X = LonG(Non0);
Y = LatG(Non0);
N = CkmG(Non0); % code of county that the grid belongs to

% store interpolated daily rain and wind in each county
DailyCountyRain = zeros(140,1);
DailyCountyWind = DailyCountyRain;

%%Interpolation
tic
DateStart = datenum(1972,8,30);
DateEnd = datenum(1972,8,31);
for i = 0:DateEnd-DateStart
    DailyPoints = rawdata(rawdata(:,5)==i+DateStart,[8 9 6 7]); % rawdata in the selected day
    XYVrain = DailyPoints(DailyPoints(:,3)~=3276.6,[1 2 3]); %delete the value 3276.6
    XYVwind = DailyPoints(DailyPoints(:,4)~=3276.6,[1 2 4]);    
    % interpolation from station to grid points
    F1 = scatteredInterpolant(XYVrain(:,1), XYVrain(:,2), XYVrain(:,3),'natural');
    Rain = F1(X, Y);
    F2 = scatteredInterpolant(XYVwind(:,1), XYVwind(:,2), XYVwind(:,3),'natural');
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

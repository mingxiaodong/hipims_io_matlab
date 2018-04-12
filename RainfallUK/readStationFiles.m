%% load MIDIAS Station information
clear,clc
stationInfor = readtable('C:\Users\b4042552\Google Drive\MyResearch\London\MIDIAS_Station_LC3.txt');
srcID = stationInfor.src_id;
% check the existence of recorded station files
list = dir('F:\Data\CEDA\RD\station_files_London\*.txt');
listSrcIDs = zeros(length(list),1);
for i= 1:length(listSrcIDs)
    listSrcIDs(i) = str2double(list(i).name(1:end-4));
end
recordedFlag = ismember(srcID,listSrcIDs);
stationInfor.recordedFlag = recordedFlag;
%% read station files to data cell
filelocation = 'F:\Data\CEDA\RD\station_files_London\';
addpath C:\Users\b4042552\Dropbox\Matlab\RainfallUK
rainfall_RD_LC3 = cell(height(stationInfor),1);
for i = 1:height(stationInfor)
    if stationInfor.recordedFlag(i)
        staionID = stationInfor.src_id(i);
        RD_T = ReadStationFile_RD(filelocation,staionID);
        rainfall_RD_LC3{i} = RD_T;
    end
end
% save rainfall_RD_LC3 stationInfor rainfall_RD_LC3
%% write table
T = table(stationInfor.src_id,stationInfor.recordedFlag);
writetable(T,'stationSrcID_RecordedFlag.txt')
%%
load data/rainfall_RD_LC3.mat
goodRainStation = struct('src_id',0,'x',0,'y',0,'recLength',0,'maxRecV',0,'dailyV',[0 0]);
n = 1;
for i=1:length(rainfall_RD_LC3)
    if ~isempty(rainfall_RD_LC3{i})
        dataTable = rainfall_RD_LC3{i};
        d = dataTable.OB_DATE;
        p = dataTable.PRCP_AMT;
        ind = p<1000&p>=0;
        d = d(ind);
        p = p(ind);
        if length(p)>365*30
        M = max(p);
        goodRainStation(n).src_id = stationInfor.src_id(i);
        goodRainStation(n).x = stationInfor.X_Coor(i);
        goodRainStation(n).y = stationInfor.Y_Coor(i);
        goodRainStation(n).recLength = length(p);
        goodRainStation(n).maxRecV = M;
        goodRainStation(n).dailyV = table(d,p);
        n = n+1;
        end
    end
end
save data/goodRainStation goodRainStation
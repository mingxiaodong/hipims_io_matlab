%*****this script is to generate daily rainfall files for each station from
%yearly daily rainfall files
clear,clc
yearFileAddress = 'E:\Data\CEDA\RD\yearly_files\'; % address of the yearly files
% stationFileAddress = 'E:\Data\CEDA\RD\station_files_LC4\'; %address of the station files to be written
cd(yearFileAddress)
% fileList = dir('*.txt'); % information struct of the yearly files
stationInfor = readtable('C:\Users\b4042552\Google Drive\MyResearch\London\MIDIAS_Station_LC4_D.txt');
% dailyRec = false(height(stationInfor),1);
% hourlyRec = false(height(stationInfor),1);
% stationInfor.dailyRec = dailyRec;
% stationInfor.hourlyRec = hourlyRec;
%%
shpLC4 = shaperead('C:/Users/b4042552/Google Drive/MyResearch/London/LondonCatchment4.shp');
StationsInterested = stationInfor;
figure
mapshow(shpLC4,'FaceColor','none','EdgeColor','b')
mapshow(StationsInterested.X,StationsInterested.Y,'DisplayType','Point')
axis image
%%
clc
yearFileAddress = 'E:\Data\CEDA\RH\yearly_files\'; % address of the yearly files
stationFileAddress = 'E:\Data\CEDA\RH\station_files_LC4\'; %address of the station files to be written
cd(yearFileAddress)
fileList = dir('*.txt'); % information struct of the yearly files
srcID_StationInterested = StationsInterested.src_id;
recFound = false(height(stationInfor),1);
srcID_yearfile = cell(length(fileList) ,1); % to store srcID of each yearly file
column_src_yearlyfile = 7; % column of the src_ID 7: hourly rain
% column_src_yearlyfile = 8; % column of the src_ID 8: daily rain;
for i = 1:length(fileList) 
    tic
    fileName = fileList(i).name; % name of the yearly file   
    %%***read the src_ID from the yearly file
    
        % all the other columns are omitted with this format
    formatSpec = [repmat('%*s',[1,column_src_yearlyfile-1]) '%f%*[^\n]'];
    delimiter = ',';
    fileID = fopen([yearFileAddress,fileName]);
        % read the column of src_id
    srcID_all = textscan(fileID, formatSpec, 'Delimiter', delimiter);
    fclose(fileID);
    srcID_all = srcID_all{1}; % all the station IDs record in the yearly file
    srcID_unique = unique(srcID_all); % unique station IDs in the yearly file              
    %%***read data into lines string
    formatSpec = '%s'; delimiter = '';
    fileID = fopen([yearFileAddress,fileName]);
        % each line is stored as a long string
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter);
    fclose(fileID);
    dataArray = dataArray{1};
    %%***add data to text files for each station
    ind_StationFound = false(size(srcID_unique));
    for j = 1:length(srcID_unique)
        srcID = srcID_unique(j);
        ind_StationInterested = srcID_StationInterested==srcID;
        if sum(ind_StationInterested)>=1
            recFound(ind_StationInterested) = true;
            ind_StationFound(j) = 1; % means the srcID is whant I am interested in
            stationFileName = [num2str(srcID) '.txt'];
            dataArrayToWrite = dataArray(srcID_all==srcID);
            % open an exist staion file or cread to a new one
            fileID = fopen([stationFileAddress stationFileName],'a');
            % attach the
            fprintf(fileID,'%s\n',dataArrayToWrite{:});
            fclose(fileID);
        end
    end
    srcID_interested_1year = srcID_unique(ind_StationFound);
    srcID_yearfile{i} = srcID_interested_1year;
    disp(['Year' fileName(16:19) ':' num2str(numel(srcID_interested_1year)) ' stations have been found.'])
    toc
end
stationInfor.hourlyRec = recFound;
%%
clc
figure
mapshow(shpLC4,'FaceColor','none','EdgeColor','b')
mapshow(stationInfor.X_Coor,stationInfor.Y_Coor,'DisplayType','Point')
% symbolspec = makesymbolspec('Point',{'Default',});
mapshow(stationInfor.X_Coor(stationInfor.hourlyRec),stationInfor.Y_Coor(stationInfor.hourlyRec),...
    'DisplayType','point','Marker','o','MarkerEdgeColor','k')
axis image
%%
writetable(stationInfor,'C:\Users\b4042552\Google Drive\MyResearch\London\MIDIAS_StationFound_LC4.txt')
%*****this script is to generate daily rainfall files for each station from
%yearly daily rainfall files
clear,clc
yearFileAddress = 'F:\Data\CEDA\RH\yearly_files\'; % address of the yearly files
stationFileAddress = 'F:\Data\CEDA\RH\station_files_London\'; %address of the station files to be written
cd(yearFileAddress)
fileList = dir('*.txt'); % information struct of the yearly files
%%
load C:\Users\b4042552\Dropbox\Matlab\London\MIDIAS_StationLondonVicinity.mat
StationsInterested = MIDIAS_StationLondonVicinity;
srcID_StationInterested = StationsInterested.src_id;
%%
srcID_yearfile = cell(length(fileList) ,1); % to store srcID of each yearly file
column_src_yearlyfile = 7; % column of the src_ID 8: daily rain; 7: hourly rain
for i = 54%1:length(fileList) 
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
    ind = false(size(srcID_unique));
    for j = 1:length(srcID_unique)
        srcID = srcID_unique(j);
        if sum(srcID_StationInterested==srcID)
            ind(j) = 1; % means the srcID is whant I am interested in
            stationFileName = [num2str(srcID) '.txt'];
            dataArrayToWrite = dataArray(srcID_all==srcID);
            % open an exist staion file or cread to a new one
            fileID = fopen([stationFileAddress stationFileName],'a');
            % attach the
            fprintf(fileID,'%s\n',dataArrayToWrite{:});
            fclose(fileID);
        end
    end
    srcID_interested_1year = srcID_unique(ind);
    srcID_yearfile{i} = srcID_interested_1year;
    disp([fileName ' has been divided to ' num2str(numel(srcID_interested_1year)) ' station files'])
    toc
end
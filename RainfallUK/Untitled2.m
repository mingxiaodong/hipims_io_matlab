clear,clc
yearFileAddress = 'F:\Data\CEDA\RD\yearly_files\';
stationFileAddress = 'F:\Data\CEDA\RD\station_files\';
cd(yearFileAddress)
fileList = dir('*.txt');
%%
for i = 1:2%length(fileList)
    fileName = fileList(i).name;    
    %% read the src_ID
    column = 8; % column of the src_ID
    formatSpec = [repmat('%*s',[1,column-1]) '%f%*[^\n]']; delimiter = ',';
    fileID = fopen([yearFileAddress,fileName]);
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter);
    fclose(fileID);
    srcIDall = dataArray{1}; % all the station record in the yearly file
    srcIDunique = unique(srcIDall); % unique station ID in the yearly file
    %% read data in line
    formatSpec = '%s'; delimiter = '';
    fileID = fopen([yearFileAddress,fileName]);
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter);
    fclose(fileID);
    dataArray = dataArray{1};
    %% write data to text files for each station
    for j = 1:length(srcIDunique)
        srcID = srcIDunique(j);
        stationFileName = [num2str(srcID) '.txt'];
        % A = exist([stationFileAddress stationFileName],'file');
        dataArrayToWrite = dataArray(srcIDall==srcID);
        if ~isempty(dataArrayToWrite)
            fileID = fopen([stationFileAddress stationFileName],'a');
            fprintf(fileID,'%s\n',dataArrayToWrite{:});
            fclose(fileID);
        end
    end
    disp([fileName 'has been divided to station files'])
end
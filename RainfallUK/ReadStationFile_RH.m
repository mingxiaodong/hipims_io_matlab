function RH_T = ReadStationFile_RH(filelocation,srcID)
%% importRDf_T   import data from station daily rainfall txt files and transfer the data to a table
%% Initialize variables.
filename = [filelocation,'/',num2str(srcID),'.txt'];
delimiter = ',';
% 15 columns, the colum name
VariableNames = {'OB_END_TIME','ID','ID_TYPE','OB_HOUR_COUNT','VERSION_NUM',...
    'MET_DOMAIN_NAME','SRC_ID','REC_ST_IND','PRCP_AMT','PRCP_DUR',...
    'PRCP_AMT_Q','PRCP_DUR_Q','METO_STMP_TIME','MIDAS_STMP_ETIME','PRCP_AMT_J'};
fileID = fopen(filename);
C = textscan(fileID,'%s %d %s %d %d %s %d %d %d %d %d %d %s %s %d',...
    'Delimiter',delimiter);
fclose(fileID);
% T = readtable(filename,'Delimiter',delimiter,'ReadVariableNames',false);
T = table;
for i=1:length(C)
    if length(C{i})==length(C{1})
        T.(i) = C{i};
        T.Properties.VariableNames{i} = VariableNames{i};
    end
end
% T.Properties.VariableNames = VariableNames;
T.OB_END_TIME = datetime(T.OB_END_TIME,'InputFormat','yyyy-MM-dd HH:mm');
T.METO_STMP_TIME = categorical(T.METO_STMP_TIME);
% T.PRCP_AMT_J = categorical(T.PRCP_AMT_J);
RH_T = T;
end

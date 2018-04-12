function RD_T = ReadStationFile_RD(filelocation,srcID)
%% importRDf_T   import data from station daily rainfall txt files and transfer the data to a table
%% Initialize variables.
filename = [filelocation,'\',num2str(srcID),'.txt'];
delimiter = ',';
% 15 columns, the colum name
VariableNames = {'ID','ID_TYPE','OB_DATE','VERSION_NUM','MET_DOMAIN_NAME',...
 'OB_END_CTIME','OB_DAY_CNT','SRC_ID','REC_ST_IND','PRCP_AMT','OB_DAY_CNT_Q'...
 ,'PRCP_AMT_Q','METO_STMP_TIME','MIDAS_STMP_ETIME','PRCP_AMT_J'};
formatSpec = '%d%s%{yyyy-MM-dd HH:mm}D%d%s%d%d%d%d%f%d%d%{yyyy-MM-dd HH:mm}D%C%d%*[^\n]';
T = readtable(filename,'Format',formatSpec,'Delimiter',delimiter,'ReadVariableNames',false);
T.Properties.VariableNames = VariableNames;
T.OB_DATE = datetime(T.OB_DATE,'InputFormat','yyyy-MM-dd HH:mm');
T.METO_STMP_TIME = categorical(T.METO_STMP_TIME);
T.PRCP_AMT_J = categorical(T.PRCP_AMT_J);
RD_T = T;
end

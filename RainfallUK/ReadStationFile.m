function RainTable = ReadStationFile(filelocation,srcID,filetype)
%% importRDf_T   import data from station daily rainfall txt files and transfer the data to a table
%% Initialize variables.
filename = [filelocation,'/',num2str(srcID),'.txt'];
delimiter = ',';
% 15 columns, the colum name
if strcmp(filetype,'day')
    VariableNames = {'ID','ID_TYPE','OB_DATE','VERSION_NUM','MET_DOMAIN_NAME',...
        'OB_END_CTIME','OB_DAY_CNT','SRC_ID','REC_ST_IND','PRCP_AMT','OB_DAY_CNT_Q'...
        ,'PRCP_AMT_Q','METO_STMP_TIME','MIDAS_STMP_ETIME','PRCP_AMT_J'};
    formatSpec = '%d%s%{yyyy-MM-dd HH:mm}D%d%s%d%d%d%d%f%d%d%s%C%C%*[^\n]';
    T = readtable(filename,'Format',formatSpec,'Delimiter',delimiter,'ReadVariableNames',false);
    T.Properties.VariableNames = VariableNames;
elseif strcmp(filetype,'hour')
    VariableNames = {'OB_END_TIME','ID','ID_TYPE','OB_HOUR_COUNT','VERSION_NUM',...
    'MET_DOMAIN_NAME','SRC_ID','REC_ST_IND','PRCP_AMT','PRCP_DUR',...
    'PRCP_AMT_Q','PRCP_DUR_Q','METO_STMP_TIME','MIDAS_STMP_ETIME','PRCP_AMT_J'};
    formatSpec = ['%{yyyy-MM-dd HH:mm}D%d%s%d%d',...
        '%s%d%d%f%f%d',...
        '%f%f%s%C%C%*[^\n]'];
    T = readtable(filename,'Format',formatSpec,'Delimiter',delimiter,'ReadVariableNames',false);
    T.Properties.VariableNames(1:15) = VariableNames;
end

% T.OB_DATE = datetime(T.OB_DATE,'InputFormat','yyyy-MM-dd HH:mm');
% T.METO_STMP_TIME = categorical(T.METO_STMP_TIME);
% T.PRCP_AMT_J = categorical(T.PRCP_AMT_J);
RainTable = T;
end

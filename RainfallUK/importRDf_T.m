function RD_T = importRDf_T(year,filelocation,stationNum)
%% importRDf_T   import data from daily rainfall txt files and transfer the data to a table
%   RD_T = importRDf_T(year,fileaddress,stationNum) import daily rianfall data as a
%   table, year is the recorded year, filelocation is the location of txt
%   file, stationNum is the SRC_ID of Met station
% See also importRDf_table
%
% created by Xiaodong Ming on 2015/12/7
%% Initialize variables.
filehead = 'midas_raindrnl_';
filetail = '.txt';
filename = [filelocation, filehead, num2str(year),'01-',num2str(year),'12',filetail];
delimiter = ',';
% 15 columns, the colum name
VariableNames = {'ID','ID_TYPE','OB_DATE','VERSION_NUM','MET_DOMAIN_NAME',...
 'OB_END_CTIME','OB_DAY_CNT','SRC_ID','REC_ST_IND','PRCP_AMT','OB_DAY_CNT_Q'...
 ,'PRCP_AMT_Q','METO_STMP_TIME','MIDAS_STMP_ETIME','PRCP_AMT_J'};
T = readtable(filename,'Delimiter',delimiter,'ReadVariableNames',false);
T.Properties.VariableNames = VariableNames;
ind = T.SRC_ID == stationNum;
T = T(ind,:);
T.OB_DATE = datetime(T.OB_DATE,'InputFormat','yyyy-MM-dd HH:mm');
T.METO_STMP_TIME = categorical(T.METO_STMP_TIME);
T.PRCP_AMT_J = categorical(T.PRCP_AMT_J);
RD_T = T;
end

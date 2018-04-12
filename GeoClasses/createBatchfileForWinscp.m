%% creat cmd file for winscp
clear,clc
filename = 'G:\cudaSWEsSolver\Eden\cmdScript.cmd';
% folder in the server
numGPU = 8;%for MultiGPU
remoteDir = '/home/b1055010/GeoClasses/release/bin/Eden10mMGPUs/';% folder in the server
loacalDir = 'G:\cudaSWEsSolver\Eden\MultiGPU\';% folder in local
%% ** upload boundary files to remote
winscpStatements = cell(numGPU,1);
fileLocation = 'G:\cudaSWEsSolver\DefenceFailure\input\field';
filemask = 'failure_monitor.dat *_BC_?.dat precipitation_source_*.dat';
winscpStatements{1} = ['lcd ' fileLocation];
for i = 0:7
    winscpStatements{i+2} = sprintf(['put ',...
        filemask,...
        ' %d/input/field/'],i);
end
%% ** sychronise input folders to remote
winscpStatements = cell(numGPU,1);
for i = 0:7
    winscpStatements{i+1} = sprintf('synchronize remote %d\\input',i);
end
%% ** download output files from remote
winscpStatements = cell(numGPU,1);
% filemask = '*_gauges.dat';
% filemask = 'h_max_*.asc';
filemask = 'h_max_*.asc; *_gauges.dat';

for i=0:numGPU-1
    winscpStatements{i+1} = sprintf(['synchronize local',...
        ' -filemask=""',filemask,'"" '...
        '%d\\output %d/output'],i,i);
end
%% ** clean remote output folders
winscpStatements = cell(numGPU,1);
for i=0:7
    winscpStatements{i+1} = sprintf('mkdir %d/output rdir %d/output',i,i);
end
%%
%********************write file******************
fileID = fopen(filename,'w');
cmdFileHeadstr = {'cd C:\Program Files (x86)\WinSCP','@echo off'...
    'winscp.com /command ^',...
    '  "open mysession" ^',...
    ['  "cd ' remoteDir ' " ^'],...
    ['  "lcd ' loacalDir ' " ^']};
for i=1:length(cmdFileHeadstr)
    fprintf(fileID,'%s\n',cmdFileHeadstr{i});
end
for i=1:length(winscpStatements)
    fprintf(fileID,'  %s\n',['"' winscpStatements{i} '"^']);
end
fprintf(fileID,'  %s\n','"exit"');
fclose(fileID);
%*************************************************
%% run batchfile
clc
cd('C:\Program Files (x86)\WinSCP')
[status,cmdout] = system(filename,'-echo');
cd G:\cudaSWEsSolver\DefenceFailure
%% clean local file
for i=0:7
    cd(['G:\cudaSWEsSolver\DefenceFailure\MultiGPUs\' num2str(i)])
%     rmdir('output','s')
     mkdir('output')
end
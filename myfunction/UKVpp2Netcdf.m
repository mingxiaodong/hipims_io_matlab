function [status,cmdout] = UKVpp2Netcdf(ppFileAddress,serverLocation)
%UKVpp2Netcdf convert UKV pp file to netcdf file via xconv in a Linux
%server
% [status,cmdout] = UKVpp2Netcdf(ppFileAddress,serverLocation)
% ppFileAddress is the folder of pp files. All the files in that folder
% will be converted to nc files.
% serverLocation is the location of xconv installed, a file named
% script.tcl is also required.
% Created by Xiaodong Ming on 2017-11-08

%% *********create and execute cmd to upload pp file to server*************
winscpDir = 'C:\Program Files (x86)\WinSCP';
cmdStatements = {...
'@echo off',...
'cd C:\Program Files (x86)\WinSCP',...
'winscp.com /command ^',...
'    "open mysession" ^',...
['	  "cd ' serverLocation '" ^'],...
['	  "lcd ""' ppFileAddress '""" ^'],...
'    "put *.pp" ^',...
'    "exit"'};
filename = 'tempScript.cmd';
fileID = fopen(filename,'w');
for i=1:length(cmdStatements)
    fprintf(fileID,'  %s\n',cmdStatements{i});
end
fclose(fileID);
currentDir = pwd;
clc
cd(winscpDir)
system([currentDir '\' filename]);
cd(currentDir)
delete(filename)
%**************************************************************************

%% *********create and execute cmd to convert pp file to netcdf************
cmdStatements = {...
['cd ' serverLocation],...
'./convsh1.93 script.tcl *.pp'};
filename0 = 'tempScript.txt';
fileID = fopen(filename0,'w');
for i=1:length(cmdStatements)
    fprintf(fileID,'%s\n',cmdStatements{i});
end
fclose(fileID);
cmdStatements = {'@echo off','C:','cd "C:\Program Files (x86)\PuTTY"',...
['plink.exe -ssh ceg-gpu01.ncl.ac.uk -l b1055010 -pw xiaxilingpu -m "',...
currentDir '\' filename0 '"']};
filename1 = 'tempScript.cmd';
fileID = fopen(filename1,'w');
for i=1:length(cmdStatements)
    fprintf(fileID,'%s\n',cmdStatements{i});
end
fclose(fileID);
system([currentDir '\' filename1]);
delete(filename0)
delete(filename1)
%**************************************************************************

%% *********create and execute cmd to convert pp file to netcdf************
winscpDir = 'C:\Program Files (x86)\WinSCP';
cmdStatements = {...
'@echo on',...
'cd C:\Program Files (x86)\WinSCP',...
'winscp.com /command ^',...
'    "open mysession" ^',...
['	  "cd ' serverLocation '" ^'],...
['	  "lcd ""' ppFileAddress '""" ^'],...
'    "get *.nc" ^',...
'    "rm *.nc" ^',...
'    "rm *.pp" ^',...
'    "exit"'};
filename = 'tempScript.cmd';
fileID = fopen(filename,'w');
for i=1:length(cmdStatements)
    fprintf(fileID,'  %s\n',cmdStatements{i});
end
fclose(fileID);
currentDir = pwd;
cd(winscpDir)
[status,cmdout] = system([currentDir '\' filename],'-echo');
cd(currentDir)
delete(filename)
end


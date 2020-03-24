function UpdatePrecipitationAll(rainSourceMat,remoteDir,numGPU)
%% to create and upload precipitation_source_all.dat file to server
% create file
filename = 'precipitation_source_all.dat';
numOfRainSource = size(rainSourceMat,2)-1;
formatSpec = repmat('%.15e ',[1 numOfRainSource] );
formatSpec = ['%09.2f ' formatSpec(1:end-1) '\n'];
fileID = fopen(filename,'w');
fprintf(fileID,'%u\n',numOfRainSource);
fprintf(fileID,formatSpec,rainSourceMat');
fclose(fileID);
% upload file
OperationID = 2; 
filemask = filename;
loacalDir = pwd;
UploadOneFile(remoteDir,loacalDir,OperationID,numGPU,filemask);
cd(loacalDir)
% delete(filename)
end

function [status,cmdout] = UploadOneFile(remoteDir,loacalDir,operationType,numGPU,varargin)
%WINSCPOPERATION Execute four types of Winscp file operation as the following list: 
%%*Operation List
%   NO. Operation Name          Description
%   --  ------                  -----------
%   1   SychronizeInput         Sychronize the local 'Input' folder to remote
%   2   UploadFieldFiles        Upload field files(eg. rainfall; boundary) to remote field folder
%   3   DownloadOutput          Download output files to local
%   4   CleanRemoteOutput       Clean remote output folder
%%
currentDir = pwd;
winscpDir = 'C:\Program Files (x86)\WinSCP'; %the location of winscp 
OperationNames = {'SychronizeInput','UploadFeildFiles','DownloadOutput','CleanRemoteOutput'};
filemask = '*.dat';
if ischar(operationType)
    OperationID = strcmp(operationType,OperationNames);
    OperationID = find(OperationID == 1);
else 
    OperationID = operationType;
end
if numel(OperationID)~=1
   error('operationType should be an integer 1~4 or a string')
end
if isempty(loacalDir)
    loacalDir = pwd;
end
cd(loacalDir)
folderMask = 'output';
if length(varargin) == 1
    filemask = varargin{1}; %'failure_monitor.dat *_BC_?.dat precipitation_source_*.dat'
elseif length(varargin) == 2
    filemask = varargin{1};
    folderMask = varargin{2};
end
if strcmp(folderMask,'field')
    folderMask = 'input/field/';
elseif strcmp(folderMask,'mesh')
    folderMask = 'input/mesh/';
end
%% ********************creat cmd file for winscp
winscpStatements = cell(numGPU,1);
switch OperationID
    case 1 % SychronizeInput
        if numGPU==1 %single GPU
            winscpStatements = {'synchronize remote input input',...
                'put times_setup.dat input/*'};
        else %multi-GPU
            winscpStatements = {sprintf('synchronize remote')}; 
        end
    case 2 % UploadFieldFiles
        if numGPU==1 %single GPU
            % upload from local single field folder to remote single field folder
            winscpStatements = {['synchronize remote -filemask=""',filemask,...
                '"" input\field input/field/']};
        else %multi GPU
            if exist('0','dir') == 7 % '0' is a folder
                % upload from local MultiFolder to remote multiFolder
                for i = 1:numGPU
                    % go to local field folder
%                     winscpStatements{2*i-1} = sprintf('lcd %s\\%d\\input\\field',loacalDir,i-1);
                    % synchronize to remote field folder
                    winscpStatements{i} = sprintf(['synchronize remote -filemask=""',...
                        filemask,'"" %d\\input\\field %d/input/field'],i-1,i-1);
                end
            else % upload from local single-folder to remote multi-folder
%                 winscpStatements{1} = ['lcd ' loacalDir '\input\field'];
                for i = 1:numGPU
                    winscpStatements{i} = sprintf(['put ',...
                        filemask,'  %d/input/field/*'],i-1);
                end
            end
        end
    case 3 % Download
        if numGPU==1 %single GPD
           winscpStatements = {sprintf(['synchronize local',...
                    ' -filemask=""',filemask,'"" '...
                    folderMask ' ' folderMask])};
        else %multi GPU
            for i=1:numGPU
                if strcmp(filemask,'DEM.txt')
                    lcdName = [num2str(i-1) '\input\mesh' ];
                    winscpStatements{i} = ['synchronize local',' -filemask=""',filemask,'"" '...
                    lcdName ' ' num2str(i-1) '/input/mesh'];
                else
                    lcdName = [num2str(i-1) '\' folderMask];
                    winscpStatements{i} = ['synchronize local',' -filemask=""',filemask,'"" '...
                    lcdName ' ' num2str(i-1) '/' folderMask];
                end
                delete([loacalDir '\' lcdName '\' filemask])
            end
        end
    case 4 % CleanRemoteOutput
        if numGPU==1 %single GPD
            winscpStatements{1} = 'rmdir output';
            winscpStatements{2} = 'mkdir output';
        else %multi GPU
            for i=1:numGPU
                winscpStatements{2*i-1} = sprintf('rmdir %d/output',i-1);
                winscpStatements{2*i} = sprintf('mkdir %d/output',i-1);
            end
        end
end

%********************write cmd file******************
filename = [currentDir '/tempScript.cmd'];
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
cd(winscpDir)
[status,cmdout] = system(filename,'-echo');
delete(filename)
cd(currentDir)

end


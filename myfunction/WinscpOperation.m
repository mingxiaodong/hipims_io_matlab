function [status,cmdout] = WinscpOperation(remoteDir,loacalDir,operationType,numGPU,varargin)
%WINSCPOPERATION Execute four types of Winscp file operation as the following list: 
%%*Operation List
%   NO. Operation Name          Description
%   --  ------                  -----------
%   1   SychronizeInput         Sychronize the local 'Input' folder to remote
%   2   UploadFieldFiles        Upload field files(eg. rainfall; boundary) to remote field folder
%   3   DownloadOutput          Download output files to local
%   4   CleanRemoteOutput       Clean remote output folder
%  [status,cmdout] = WinscpOperation(remoteDir,loacalDir,operationType,numGPU,varargin)
%   remoteDir is the folder directory in server. loacalDir is the folder
%   directory in local PC. operationType is the type of operation, its
%   value could be a string or a integer from 1~4. numGPU is the number of
%   GPUs will be used to run model. varargin is an optional parameter to
%   give the value of filemask, for example: '*.asc; *_BC_?.dat'
%   
% Created by Xiaodong Ming 2017-10-18
% Updated by Xiaodong Ming 2018-01-09
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
    case 1 % Upload the Input Folder and files
        if numGPU==1 %single GPU
            winscpStatements = {'synchronize remote input input'};
        else %multi-GPU
            for i=1:numGPU
                loaclFolderName = [loacalDir '\' num2str(i-1) '\output'];
                delete([loaclFolderName '\*'])
            end
            winscpStatements = {sprintf('synchronize remote')};
        end
    case 2 % Upload Files
        if numGPU==1 % single GPU
            % upload from local single field folder to remote single field folder
            winscpStatements = {['synchronize remote -filemask=""',filemask,...
                '"" input\field input/field/']};
        else % multi GPU
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
                for i = 1:numGPU
                    winscpStatements{i} = sprintf(['synchronize remote -filemask=""',...
                        filemask,'"" input\\field %d/input/field/'],i-1);
                end
            end
        end
    case 3 % Download
        if numGPU==1 %single GPD
           winscpStatements = {sprintf(['synchronize local',...
                    ' -filemask=""',filemask,'"" '...
                    folderMask ' ' folderMask])};
        else %multi GPU
            n = 1;
            for i=1:numGPU
                if strcmp(filemask,'DEM.txt')
                    lcdName = [num2str(i-1) '\input\mesh' ];
                    winscpStatements{i} = ['get ',num2str(i-1),'/input/mesh/', filemask, ' ', lcdName '\*' ];
                else
                    lcdName = [num2str(i-1) '\' folderMask];
                    if contains(filemask,'gauges') %if download gauges data
                        winscpStatements{n} = ['get ',num2str(i-1),'/input/field/gauges* ', num2str(i-1),'\input\field\*' ];
                        n = n+1;
                    end                   
                    filemask_Str = split(filemask,';');
                    for nMask = 1:length(filemask_Str)
                        filemask_Str_1 = filemask_Str{nMask};
                        winscpStatements{n} = ['get ',num2str(i-1),'/output/', filemask_Str_1, ' ', lcdName '\*' ];
                        n = n+1;
                    end
                end
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
cmdFileHeadstr = {['cd ' winscpDir],'@echo off'...
    'winscp.com /command ^',...
    '  "open mysession" ^',...
    ['  "cd ""' remoteDir '"" " ^'],...
    ['  "lcd ""' loacalDir '"" " ^']};
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


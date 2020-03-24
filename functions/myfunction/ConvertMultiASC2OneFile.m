function [Z_all,R_all] = ConvertMultiASC2OneFile(caseFolder,resultFileName,outputName)
%% CombineMultiGPUResults Summary of this function goes here
%   Detailed explanation goes here
% Section number must be sequenced from bottom to top as (numGPU-1:-1:0)
% in the boundary between each section, there are two lines overlayed
list = dir(caseFolder);
list(~[list.isdir])=[];
folderNum = nan(1,length(list));
for i=1:length(list)
    folderNum(i) = str2double(list(i).name);
end
folderNum(isnan(folderNum))=[];
folderNum_max = max(folderNum);

numSec = numel(folderNum);
x11Vec = nan(numSec,1);
y11Vec = x11Vec;
cellsizeVec = x11Vec;
nrowsVec =  x11Vec;
ncolsVec = x11Vec;
% read the reference data of all asc files
for i=folderNum
    oneFileName = [caseFolder '\'  num2str(i) '\output\' resultFileName];
    [x11,y11,cellsize,nrows,ncols] = ReadAscHead(oneFileName);
    x11Vec(i+1) = x11;
    y11Vec(i+1) = y11;
    cellsizeVec(i+1) = cellsize;
    nrowsVec(i+1) = nrows;
    ncolsVec(i+1) = ncols;
end
x11_all =  unique(x11Vec);
y11_all = max(y11Vec);
cellsize_all = unique(cellsizeVec);
ncols_all = unique(ncolsVec);
nrows_all = sum(nrowsVec)-2*(numSec-1);
xllcorner = x11_all - cellsize_all/2;
yllcorner = y11_all - cellsize_all*nrows_all + cellsize_all/2;
% write asc head rows
WriteAscHead(outputName,ncols_all,nrows_all,xllcorner,yllcorner,cellsize_all)
R_all = makerefmat(x11_all,y11_all,cellsize_all,-cellsize_all);
[~,indSeq] = sort(y11Vec,'descend');
sectionSeq = folderNum(indSeq);
%% wirte matrix
for i=sectionSeq
        oneFileName = [caseFolder '\'  num2str(i) '\output\' resultFileName];
        [Z,~] = ArcgridreadM(oneFileName);
        if i==folderNum_max
            dlmwrite(outputName,Z,'-append','delimiter','\t')
        else
            dlmwrite(outputName,Z(3:end,:),'-append','delimiter','\t')
        end
end
%%
if nargout==2
    [Z_all,~] = ArcgridreadM(outputName);
end
end
%%
function [x11,y11,cellsize,nrows,ncols] = ReadAscHead(filename)
fileID = fopen(filename);
if fileID == -1
    error([filename ' cannot be openned'])
end
C = textscan(fileID,'%s %f',6);
ncols = C{2}(1);
nrows = C{2}(2);
xllcorner = C{2}(3);
yllcorner = C{2}(4);
cellsize = C{2}(5);
x11 = xllcorner+0.5*cellsize;
y11 = yllcorner+(nrows-0.5)*cellsize;
fclose(fileID);
end

function WriteAscHead(filename,ncols,nrows,xllcorner,yllcorner,cellsize)
fileID = fopen(filename,'w');
fprintf(fileID,'ncols    %d\n', ncols);
fprintf(fileID,'nrows    %d\n', nrows);
fprintf(fileID,'xllcorner    %f\n', xllcorner);
fprintf(fileID,'yllcorner    %f\n', yllcorner);
fprintf(fileID,'cellsize    %f\n', cellsize);
fprintf(fileID,'NODATA_value    -9999\n');
fclose(fileID);
end
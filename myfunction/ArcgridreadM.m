function [Z,R] = ArcgridreadM(filename)
%ArcgridreadM Simplified version of arcgridread. asc file must have 6 head lines
% Created by Xiaodong Ming on 2017-11-27
% See also arcgridread
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
nanValue = C{2}(6);
formatSpec = [repmat('%f',[1,ncols]) '%[^\n\r]'];
Z = textscan(fileID,formatSpec,nrows);
fclose(fileID);
Z = [Z{1:end-1}];
Z(Z==nanValue) = nan;
x11 = xllcorner+0.5*cellsize;
y11 = yllcorner+(nrows-0.5)*cellsize;
R = makerefmat(x11,y11,cellsize,-cellsize);
end
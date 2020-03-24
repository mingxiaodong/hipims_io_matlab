function Arcgridwrite(filename,Z,R)
%ARCGRIDWRITE Write gridded data set to Arc ASCII Grid Format
%
%   ARCGRIDWRITE(filename,Z,R) write a matlab format grid with its 
%   reference into an Arc ASCII file Grid format.  
%   Z is a 2D array containing the data values.  R is a
%   referencing matrix (see MAKEREFMAT).  -9999 is assigned to elements 
%   of Z corresponding to null data values in the grid file.
%
%   See also ARCGRIDREAD.

% Created by Ming 2017-06-09.

% x11 and y11 specify the map location of the center of the first (1,1) 
% pixel in the image or the first element of the data grid
% x00 and y00 specify the centre of the upper-left fake pixel of (x11,y11)
x00 = R(3,1); 
y00 = R(3,2);
gridSize = abs(R(2));
nrows = size(Z,1);
ncols = size(Z,2);

% XLLCORNER and YLLCORNER are the coordinates of the
% lower left corner of the lower left cell of Z.
xllcorner = x00 + gridSize/2;
yllcorner = y00 - gridSize*nrows - gridSize/2;
NODATA_value = -9999;
Z(isnan(Z)) = NODATA_value;

fileID = fopen(filename,'w');
fprintf(fileID,'ncols    %d\n', ncols);
fprintf(fileID,'nrows    %d\n', nrows);
fprintf(fileID,'xllcorner    %f\n', xllcorner);
fprintf(fileID,'yllcorner    %f\n', yllcorner);
fprintf(fileID,'cellsize    %f\n', gridSize);
fprintf(fileID,'NODATA_value    -9999\n');
dlmwrite(filename,Z,'-append','delimiter','\t')
fclose(fileID);
end
function [X, Y] = Ind2Map(ind,siz,R)
% convert index of an array with reference mat to its map cooridates X,Y
% ind: index of the array
% siz: size of the array
% R: reference matrix of the array
[row,col] = ind2sub(siz,ind);
[x11,y11] = refMat2firstPixel(R);
cellsize = abs(R(2));
X = x11+(col-1)*cellsize;
Y = y11-(row-1)*cellsize;
end

function [x11,y11] = refMat2firstPixel(R)
% convert refmat to the coordinate at the centre of the first pixel
x00 = R(3,1);
y00 = R(3,2);
dx = R(2,1);
dy = R(1,2);
x11 = x00+dx;
y11 = y00+dy;
end
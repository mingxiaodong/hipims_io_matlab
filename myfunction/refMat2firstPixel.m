function [x11,y11] = refMat2firstPixel(R)
% convert refmat to the coordinate of the first pixel
x00 = R(3,1);
y00 = R(3,2);
dx = R(2,1);
dy = R(1,2);
x11 = x00+dx;
y11 = y00+dy;
end
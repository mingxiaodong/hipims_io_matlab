function [X,Y] = Raster2FeaturePoints(Z,R)
%[X,Y] = RASTER2FEATUREPOINTS(Z,R) Conver raster grid to feature points
%
%   [X,Y] = Raster2FeaturePoints(Z,R) rerurns the X and Y coordinates of
%   the feature points. Z is the raster grid matrix and R is the reference
%   matrix.
%   Created by Xiaodong Ming on 2017-08-02
%   See also MAKEREFMAT, ARCGRIDREAD.
dx = R(2,1); dy = R(1,2);
[m,n] = size(Z);
x11 = R(3,1)+dx;
y11 = R(3,2)+dy;
x = x11:dx:x11+dx*(n-1);
y = y11:dy:y11+dy*(m-1);
[X,Y] = meshgrid(x,y);
end
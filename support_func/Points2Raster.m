function newZ = Points2Raster(sizeZ,R,PointsVec,method)
%% convert point values to a raster with interpolating methods
% sizeZ: the size of targetted raster array
% R: the reference mat of targetted raster
% PointsVec: (x,y,v) a 3-col vector gives the x,y coordiantes and values
% method: interpolate methods (nearest, linear,)
% return an array with the same size of raster
% Created by Xiaodong Ming on 17-Oct-2019
[X,Y] = Raster2FeaturePoints(sizeZ,R);
F = scatteredInterpolant(PointsVec(:,1),PointsVec(:,2),PointsVec(:,3),method);
newZ = F(X,Y);
end
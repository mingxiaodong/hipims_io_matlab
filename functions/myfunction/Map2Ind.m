function outputInd = Map2Ind(X,Y,sizeZ,R)
% OutputInd = Map2Ind(X,Y,[nrows,ncols],R) convert map coordiantes XY 
% to the index of the map matrix
%       X and Y: points map coordiantes
%       [nrows,ncols]: 2-integer vector, size of domain matrix
%       R: reference mat of the domain matrix
% Created by Xiaodong Ming on 11/3/2018
[row,col] = map2pix(R,X(:),Y(:));
row = round(row); col = round(col);
row(row>sizeZ(1)|row<1) = nan;
col(col>sizeZ(2)|col<1) = nan;
outputInd = sub2ind(sizeZ,row,col);
outputInd = reshape(outputInd,size(X));
end


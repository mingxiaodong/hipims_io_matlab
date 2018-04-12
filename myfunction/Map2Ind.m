function outputInd = Map2Ind(X,Y,sizeZ,R)
%outputInd = Map2Ind(X,Y,sizeNewZ,newR) convert map coordiantes XY to index of Z matrix of the map
%   xy: 2col array X and Y map coordiantes
%   sizeZ 2-integer vector, size of Z mat [rows cols]
%   R: reference mat of map
%   Detailed explanation goes here
% Created by Xiaodong Ming on 11/3/2018
[row,col] = map2pix(R,X(:),Y(:));
row = round(row); col = round(col);
row(row>sizeZ(1)|row<1) = nan;
col(col>sizeZ(2)|col<1) = nan;
outputInd = sub2ind(sizeZ,row,col);
outputInd = reshape(outputInd,size(X));
end


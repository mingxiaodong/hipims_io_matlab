function [z_New,r_New] = ClipDEM(Z,R,frameCoords)
%[z_New,r_New] = ClipDEM(Z,R,frameCoords) clip DEM based on the frame
%   defined by frameCoords. Z is the data mat and R is the reference mat.
%   frameCoords is a two column vector with x and y coords of the frame
% Created by Xiaodong Ming on 2017-08-04
x = frameCoords(:,1);
y = frameCoords(:,2);
[row,col] = map2pix(R,x,y);
min_row = min(round(row)); max_row = max(round(row));
min_col = min(round(col)); max_col = max(round(col));
if min_row<1
    min_row=1;
end
if max_row>size(Z,1)
    max_row=size(Z,1);
end
if min_col<1
    min_col=1;
end
if max_col>size(Z,2)
    max_col=size(Z,2);
end
row = min_row:max_row;
col = min_col:max_col;
z_New = Z(row,col);
r_New = R;
r_New(3,1) = r_New(3,1)+(min_col-1)*r_New(2,1);
r_New(3,2) = r_New(3,2)+(min_row-1)*r_New(1,2);
end
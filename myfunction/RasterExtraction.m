function [Z_new,R_new] = RasterExtraction(Z,R,mask)
%RASTEREXTRACTION Extract a raster Z by a mask array
%   [Z_new,R_new] = RasterExtraction(Z,R,mask) Extract a raster Z by a mask
%       vectorand return new data matrix Z_new and reference matrix R_new.
%   Z and R: the value and spatial reference matrix of the raster to be extracted. 
%   mask: a two column matrix giving the X and Y coordinates of points that
%           define the extent of extraction. If mask has only two points 
%           then their diagonal rectangle will be the extraction mask. 
X = mask(:,1); Y = mask(:,2);
R_new = R;
if size(mask,1) == 2
    rectangle_xy = [min(X),min(Y);...
                    max(X),min(Y);...
                    max(X),max(Y);...
                    min(X),max(Y);
                    min(X),min(Y)];
    X = rectangle_xy(:,1);
    Y = rectangle_xy(:,2);
    [row_range,col_range] = map2pix(R,X([4,2]),Y([4,2]));%top-left&bottom-right
    row_range = round(row_range); row_range(row_range<1)=1;row_range(row_range>size(Z,1)) = size(Z,1);
    col_range = round(col_range); col_range(col_range<1)=1;col_range(col_range>size(Z,2)) = size(Z,2);
    row_min = row_range(1); col_min = col_range(1);
    row_max = row_range(2); col_max = col_range(2);
    Z_new = Z(row_min:row_max,col_min:col_max);
    R_new(3,1) = R(3,1)+ R(2,1)*(col_min-1);
    R_new(3,2) = R(3,2)+ R(1,2)*(row_min-1);
else
    [XX_Z,YY_Z] = Raster2FeaturePoints(Z,R);
    in = inpolygon(XX_Z(:),YY_Z(:),X,Y);
    ind_new = 1:length(in);
    Z(~in) = nan;
    [row_mask,col_mask] = ind2sub(size(Z),ind_new(in));
    row_mask = round(row_mask); col_mask = round(col_mask);
    row_mask(row_mask<1)=1; row_mask(row_mask>size(Z,1)) = size(Z,1);
    col_mask(col_mask<1)=1; col_mask(col_mask>size(Z,2)) = size(Z,2);
    Z_new = Z(min(row_mask):max(row_mask),min(col_mask):max(col_mask));
    R_new(3,1) = R(3,1)+ R(2,1)*(min(col_mask)-1);
    R_new(3,2) = R(3,2)+ R(1,2)*(min(row_mask)-1);
end
end

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

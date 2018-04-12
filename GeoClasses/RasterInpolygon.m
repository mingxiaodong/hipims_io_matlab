function  in = RasterInpolygon(Z,R,xv,yv)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
ind = isnan(xv);
x_vertex = xv(ind);
y_vertex = yv(ind);
x_all = x_vertex;
y_all = y_vertex;
n = numel(x_all);
gridSize = abs(R(2));
indZ = nan(size(Z));
for j =1:numel(x_vertex)-1 %number of vertex in one polygon
    x1 = x_vertex(j); x2 = x_vertex(j+1);
    y1 = y_vertex(j); y2 = y_vertex(j+1);
    if x1 ~= x2
        xy_slope = (y2-y1)/(x2-x1);
        xy_cosin = 1/sqrt(1+xy_slope^2);
        gridGap = xy_cosin*gridSize*(x2-x1)/abs(x1-x2);
        x_segment = x1:gridGap:x2;
        y_segment = y1 + (x_segment-x1)*xy_slope;
    else %x1==x2 means a vertical line segment
        if y1==y2 %means they are duplicate points, not a line segment
            continue
        else
            gridGap = gridSize*(y2-y1)/abs(y1-y2);
            y_segment = y1:gridGap:y2;
            x_segment = x1 + y_segment*0;
        end
    end
    x_all(n+1:n+numel(x_segment)) = x_segment;
    y_all(n+1:n+numel(x_segment)) = y_segment;
    n = numel(x_all);
end
pix_RowCol = map2pix(R,x_all, y_all); % cols and rows
pix_RowCol = round(pix_RowCol); 
pix_RowCol(pix_RowCol==0)=1; % check the rows and cols out of Z range
pix_RowCol(pix_RowCol(:,1)>size(Z,1),1)= size(Z,1);
pix_RowCol(pix_RowCol(:,2)>size(Z,2),2)= size(Z,2);
pix_Ind = sub2ind(size(Z),pix_RowCol(:,1),pix_RowCol(:,2));
pix_Ind = unique(pix_Ind,'stable'); % indice converted from cols and rows and sorted
indZ(pix_Ind) = 1;
for i=1:size(indZ,1)
    rowValue = indZ(i,:);
    indcol = find(indZ==1);
    nPts = numel(indcol);
    if nPts<=1
        continue
    elseif round(nPts/2) == nPts/2 %even number
        ind1row = zeros(1,2);
        n = 1;
        for j = 1:nPts/2
            ind1 = nPts(2*j-1):nPts(2*j);
            ind1row(n:n+numel(ind1)-1) = ind1;
            n = n+numel(ind1);
        end
    end
    rowValue(ind1row) = 1;
    indZ(i,:) = rowValue;
end
in = indZ==1;
end
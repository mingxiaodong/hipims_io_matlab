function ind = RasterizeLine(xq,yq,R,sizeZ)
% rasterize a line vector on a raster grid
% xq,yq: coordinates of the line vertex
% sizeZ: the size of the grid
% R: the reference matrix of the grid
 
n = sum(~isnan(xq))-1; %number of vertex of the line
subs_Cell = cell(n,2);
[rows,cols] = map2pix(R,xq,yq);
for i=1:n
    rowsDiff = abs(rows(i)-rows(i+1));
    colsDiff = abs(cols(i)-cols(i+1));
    nPoints = max(rowsDiff,colsDiff)+1;
    nPoints = ceil(nPoints);
    X = linspace(xq(i),xq(i+1),nPoints);
    Y = linspace(yq(i),yq(i+1),nPoints);
    [row,col] = map2pix(R,X',Y');
    row = round(row);
    col = round(col);
    row(row>sizeZ(1))=nan;
    col(col>sizeZ(2))=nan;
    row(row<1)=nan;
    col(col<1)=nan;
    subs_Cell{i,1} = row;
    subs_Cell{i,2} = col;
end
subs_Cell = cell2mat(subs_Cell);
ind = sub2ind(sizeZ,subs_Cell(:,1),subs_Cell(:,2));
ind = unique(ind);
ind(isnan(ind))=[];
end
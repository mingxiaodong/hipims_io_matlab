%% build one uniform km Grid represented by XY points in the polygon
clear,clc
YRD = shaperead('/Users/Sheldon1/Filr/My Files/Data/YRD bound/YRD_county.shp');
MinX = YRD(1).BoundingBox(1);
MaxX = YRD(1).BoundingBox(2);
MinY = YRD(1).BoundingBox(3);
MaxY = YRD(1).BoundingBox(4);
for i = 1:140
    if MinX > YRD(i).BoundingBox(1);
        MinX = YRD(i).BoundingBox(1);
    end
    if MaxX < YRD(i).BoundingBox(2);
        MaxX = YRD(i).BoundingBox(2);
    end
    if MinY > YRD(i).BoundingBox(3);
        MinY = YRD(i).BoundingBox(3);
    end
    if MaxY < YRD(i).BoundingBox(4);
        MaxY = YRD(i).BoundingBox(4);
    end
end
YRD_bound = [MinX,MinY;MaxX,MaxY];
Xran = MinX:1000:MaxX;
Yran = MinY:1000:MaxY;
[XkmG,YkmG] = meshgrid(Xran,Yran);
XkmV = XkmG(:);
YkmV = YkmG(:);
CkmV = zeros(size(XkmV)); % record which polygon the point belong to
% verify each point is inside the polygon or not
for i = 1:140
    xv = YRD(i).X;
    yv = YRD(i).Y;
    IN = inpolygon(XkmV,YkmV,xv,yv); 
    %zero means the point is not inside any of the polygon
    CkmV(IN) = i;
end
CkmG = reshape(CkmV,size(XkmG)); % Grid data representing points belonging
YRD_KmsGrid = [XkmV,YkmV,CkmV];
% A uniform grid of YRD, column 3 represents which polygon points belong to
YRD_KmsGrid = YRD_KmsGrid(YRD_KmsGrid(:,3)>0,:); 
mapshow(YRD_KmsGrid(:,1),YRD_KmsGrid(:,2),'DisplayType','point')
%mapshow(YRD,'DisplayType','polygon','FaceColor','none');
%mapshow(XkmG,YkmG,CkmG,'DisplayType', 'mesh')
%save YRD_KmsGrid YRD_KmsGrid XkmG YkmG CkmG
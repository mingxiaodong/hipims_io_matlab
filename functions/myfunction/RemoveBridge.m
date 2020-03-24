function Z0 = RemoveBridge(Z,R,shp)
% remove bridges on DEM via interpolation methods
% Z,R: height and reference matrix of DEM
% shp: polygons define the position of bridges
%     each feature has four vertex, vertex 1,vertex 2 are in the same side of the
%     bridge,  vertex 3 and vertex 4 are in the other side
% return New DEM height matrix
% functions required:
%     RasterizeLine, GenerateBedElevationRaster,
%     Raster2FeaturePoints, Map2Ind
% Created by Xiaodong Ming on 2019-5-24
Z0 = Z;
for i=1:length(shp)
    XY = unique([shp(i).X',shp(i).Y'],'rows','stable');
    XY(isnan(XY(:,1)),:)=[]; 
    if length(XY)~=4
        warning(['Feature' num2str(i) 'ignored: not a 4-vertext polygon'])
        continue
    end        
    % identify the lines along bridge and lines across bridge
    % bankLine: flat,higher, lines across the bridge
    % crossLine: curve, lower, lines along outside edge of bridge
    X = XY(:,1);
    Y = XY(:,2);
    lineL1 = [X([1,2]),Y([1,2])]; 
    lineL2 = [X([4,3]),Y([4,3])]; 
    lineW1 = [X([1,4]),Y([1,4])];%lines across the bridge
    lineW2 = [X([2,3]),Y([2,3])];
    cellsize = R(2);
    indL1 = RasterizeLine(lineL1(:,1),lineL1(:,2),R,size(Z));
    indL2 = RasterizeLine(lineL2(:,1),lineL2(:,2),R,size(Z));
    indW1 = RasterizeLine(lineW1(:,1),lineW1(:,2),R,size(Z));
    indW2 = RasterizeLine(lineW2(:,1),lineW2(:,2),R,size(Z));
    zL1 = Z(indL1);
    zL2 = Z(indL2);
    zW1 = Z(indW1);
    zW2 = Z(indW2);
    if mean(zL1)+mean(zL2) > mean(zW1)+mean(zW2)
        bankLine1 = lineL1;
        bankLine2 = lineL2;
        crossInd1 = indW1;
        crossInd2 = indW2;
    else
        bankLine1 = lineW1;
        bankLine2 = lineW2;
        crossInd1 = indL1;
        crossInd2 = indL2;
    end
    % lines along the outside edge of the bridge
    [I,J]=ind2sub(size(Z),crossInd1);
    [crossX1,crossY1] = pix2map(R,I,J);
    crossLine1 = [crossX1,crossY1,Z(crossInd1)];    
    [I,J]=ind2sub(size(Z),crossInd2);
    [crossX2,crossY2] = pix2map(R,I,J);
    crossLine2 = [crossX2,crossY2,Z(crossInd2)];
    disp(i)
    crossSectionLine = {crossLine1,crossLine2};
    errorDistance = cellsize/2;
    [Z_new,R_new] = GenerateBedElevationRaster(bankLine1,bankLine2,crossSectionLine,cellsize,errorDistance);
    [X_new,Y_new] = Raster2FeaturePoints(Z_new,R_new);
    X_new(isnan(Z_new))=[];
    Y_new(isnan(Z_new))=[];
    Z_new(isnan(Z_new))=[];
    ind_new = Map2Ind(X_new,Y_new,size(Z),R);
    Z0(ind_new) = Z_new;
%     disp(i)
end
end
%% example
clear,clc
cd('/Users/ming/OneDrive - Newcastle University/Newcastle/SharedData')
DEMName = 'DEM2m_NewcastleFrame2.asc';
[Z,R] = ArcgridreadM(DEMName);%'DEM2m_NewcastleFrame2.asc'
shp = shaperead('Bridges.shp');
Z0 = RemoveBridge(Z,R,shp);
% mapshow(Z0,R,'DisplayType','surface')
% Arcgridwrite('DEM2m_NewcastleFrame2_BridgeRemoved.asc',Z0,R)
%% read a river DEM and create a bridge
clear,clc
[Z,R] = ArcgridreadM('fakeriver.asc');
bridgeValues = Z(:,130:135);
bridgeValues(~isnan(bridgeValues)) = max(bridgeValues(:));
% Z(:,130:135) = bridgeValues;
% figure('units','normalized','outerposition',[0 0 1 1])
figure('Color','w')
mapshow(Z,R,'DisplayType','mesh')
view([23 26])
box on
grid on
ax = gca;
ax.BoxStyle = 'full';
% export_fig riverchannel_BR.jpg -r300
%% create a four-vertex polygon covering the bridge
vert1 = [39.4,61.7];
vert2 = [39.4,52.1];
vert3 = [41.8,52.1];
vert4 = [41.8,61.7];
polygonXY = [vert1;vert2;vert3;vert4];
X = [polygonXY(:,1);polygonXY(1,1)];
Y = [polygonXY(:,2);polygonXY(1,2)];
ind = RasterizeLine(X,Y,R,size(Z));
Z0 = Z;
Z0(ind) = 20;
mapshow(Z0,R,'DisplayType','mesh')
view([23 26])
box on
grid on
ax = gca;
ax.BoxStyle = 'full';
%%
bankLine1 = [X([1,4]),Y([1,4])];
bankLine2 = [X([2,3]),Y([2,3])];
cellsize = R(2);
crossInd1 = RasterizeLine(X([1,2]),Y([1,2]),R,size(Z));
[I,J]=ind2sub(size(Z),crossInd1);
[crossX1,crossY1] = pix2map(R,I,J);
crossLine1 = [crossX1,crossY1,Z(crossInd1)];
crossInd2 = RasterizeLine(X([4,3]),Y([4,3]),R,size(Z));
[I,J]=ind2sub(size(Z),crossInd2);
[crossX2,crossY2] = pix2map(R,I,J);
crossLine2 = [crossX2,crossY2,Z(crossInd2)];
crossSectionLine = {crossLine1,crossLine2};
%%
errorDistance = 0.2;
[Z_new,R_new] = GenerateBedElevationRaster(bankLine1,bankLine2,cellsize,crossSectionLine,errorDistance);
[X_new,Y_new] = Raster2FeaturePoints(Z_new,R_new);
X_new(isnan(Z_new))=[];
Y_new(isnan(Z_new))=[];
Z_new(isnan(Z_new))=[];
ind_new = Map2Ind(X_new,Y_new,size(Z),R);
Z0 = Z;
Z0(ind_new) = Z_new;

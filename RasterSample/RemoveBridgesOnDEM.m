clear,clc
shpPolygon = shaperead('tempContour.shp');
[z1,r1] = arcgridread('dem5m_extract.asc');
%%
[X1,Y1] = Raster2FeaturePoints(z1,r1);
xv = [shpPolygon(1).X];
yv = [shpPolygon(1).Y];
in = inpolygon(X1,Y1,xv,yv);
X_interp = X1(~in);
Y_interp = Y1(~in);
V_interp = z1(~in);
xq = X1(in);
yq = Y1(in);
ind = isnan(V_interp);
X_interp(ind) = []; Y_interp(ind) = []; V_interp(ind) = [];
F = scatteredInterpolant(X_interp,Y_interp,V_interp);
vq = F(xq,yq);
z1_new = z1;
z1_new(in) = vq;
%%
figure
mapshow(shpPolygon,'FaceColor','none')
%%
figure
mapshow(z1_new,r1,'DisplayType','Surface')
mask = [540000,175000;550000, 185000; 540000,185000];
[Z_new,R_new] = RasterExtraction(Z,R,mask);
%%
figure
% mapshow(Z,R,'DisplayType','Surface')
hold on
mapshow(Z_new,R_new,'DisplayType','Surface')
axis equal
hold off
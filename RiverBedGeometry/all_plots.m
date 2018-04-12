%% show river bed points and river boundary
figure;
hold on
plot(riverBedPoints_xy(:,1),riverBedPoints_xy(:,2),'.b');
mapshow(channelBoundPoly(:,1),channelBoundPoly(:,2),'DisplayType','Polygon','FaceColor','none','EdgeColor','r');
hold off
xlabel('x');ylabel('y')
axis image
%% plot scatter points river polygon outline
figure(1);
mapshow(channalBound_section(:,1),channalBound_section(:,2),'DisplayType','Polygon','FaceColor','none');
xlabel('x');ylabel('y')
hold on;
scatter(crossLineStart(:,1),crossLineStart(:,2),'k*'); 
scatter(crossLineEnd(:,1),crossLineEnd(:,2),'b*');
scatter(bankLineTop(:,1),bankLineTop(:,2),'ro'); 
scatter(bankLineBottom(:,1),bankLineBottom(:,2),'go'); 
scatter(bankLineBottom(1,1),bankLineBottom(1,2),'rd','filled')
hold off
%% 3D plot for cross line
figure(2)
hold on
scatter3(crossLineStart(:,1),crossLineStart(:,2),crossLineStart(:,3),'k*'); 
scatter3(crossLineEnd(:,1),crossLineEnd(:,2),crossLineEnd(:,3),'b*');
hold off
xlabel('x');ylabel('y')
%% scatter plot for relative coords of outline
figure(3);
hold on;
scatter(crossLeft_lwh(:,1),crossLeft_lwh(:,2),'k*');
scatter(crossRight_lwh(:,1),crossRight_lwh(:,2),'b*');
scatter(bankTop_lw(:,1),bankTop_lw(:,2),'ro');
scatter(bankBottom_lw(:,1),bankBottom_lw(:,2),'go');
scatter(bankBottom_lw(1,1),bankBottom_lw(1,2),'r+')
hold off
xlabel('l');ylabel('w')
%% scatter plot for absolute coords of all points
figure(4);
mapshow(channalBound_section(:,1),channalBound_section(:,2),'DisplayType','Polygon','FaceColor','none');
xlabel('x');ylabel('y')
hold on;
scatter(crossLineStart(:,1),crossLineStart(:,2),'k*'); 
scatter(crossLineEnd(:,1),crossLineEnd(:,2),'b*');
scatter(bankLineTop(:,1),bankLineTop(:,2),'ro'); 
scatter(bankLineBottom(:,1),bankLineBottom(:,2),'go'); 
scatter(bedPoints(:,1),bedPoints(:,2)); 
hold off
%% scatter plot for relative coords of all points
figure(5)
hold on;
scatter(crossLeft_lwh(:,1),crossLeft_lwh(:,2),'k*');
scatter(crossRight_lwh(:,1),crossRight_lwh(:,2),'b*');
scatter(bankTop_lw(:,1),bankTop_lw(:,2),'ro');
scatter(bankBottom_lw(:,1),bankBottom_lw(:,2),'go');
scatter(bedPoints_l,bedPoints_w);
hold off
xlabel('l');ylabel('w')
%% plot scatter3 points for L, W, H
figure(6);
scatter3(bedPoints_l,bedPoints_w,bedPoints_h,'m'); 
hold on;
scatter3(crossLeft_lwh(:,1),crossLeft_lwh(:,2),crossLeft_lwh(:,3),'k*');
scatter3(crossRight_lwh(:,1),crossRight_lwh(:,2),crossRight_lwh(:,3),'b*');
xlabel('L');ylabel('W');zlabel('H')
hold off
%% plot absolute surface
i = 10;
channalBound_section = riverPolySections(i).lineAll;
crossLineStart= riverPolySections(i).line1;
bankLineBottom = riverPolySections(i).line2;
crossLineEnd= riverPolySections(i).line3;
bankLineTop = riverPolySections(i).line4;
bedPoints_section = bedPoints_cell{i};  
figure(7)
% bedPoints_all = [bedPoints, bedPoints_h;crossLineStart;crossLineEnd];
% mapshow(channalBound_section(:,1),channalBound_section(:,2),'DisplayType','polygon','Facecolor','none')
tri = delaunay(bedPoints_section(:,1),bedPoints_section(:,2));
trisurf(tri,bedPoints_section(:,1),bedPoints_section(:,2),bedPoints_section(:,3),'EdgeColor','none')
hold on
% scatter3(bedPoints(:,1),bedPoints(:,2),bedPoints_h,'m.'); %axis image
plot3(crossLineStart(:,1),crossLineStart(:,2),crossLineStart(:,3))
scatter3(crossLineStart(:,1),crossLineStart(:,2),crossLineStart(:,3),'k*');
scatter3(crossLineEnd(:,1),crossLineEnd(:,2),crossLineEnd(:,3),'b*');
xlabel('x towards East')
ylabel('y towards North')
zlabel('z elevation')
hold off
ax = gca;
ax.DataAspectRatio = [1,1,0.2];
%% plot river bed 3D points cell
figure(8)
hold on
for i=1:length(bedPoints_cell)
    data = bedPoints_cell{i};
%     mapshow(riverPolySections(i).lineAll(:,1),riverPolySections(i).lineAll(:,2),'DisplayType','Polygon','FaceColor','none');
    scatter3(data(:,1),data(:,2),data(:,3),'.')
end
hold off
xlabel('x towards East')
ylabel('y towards North')
zlabel('z elevation')
ax = gca;
ax.DataAspectRatio = [1,1,0.2];

%% plot scatter points river polygon outline
figure;
mapshow(channalBound_section(:,1),channalBound_section(:,2),'DisplayType','Polygon','FaceColor','none');
xlabel('x');ylabel('y')
hold on;
scatter(crossLineStart(:,1),crossLineStart(:,2),'k*'); 
scatter(crossLineEnd(:,1),crossLineEnd(:,2),'b*');
scatter(bankLineTop(:,1),bankLineTop(:,2),'ro'); 
scatter(bankLineBottom(:,1),bankLineBottom(:,2),'go'); 
hold off
% %% plot DEM and river poly
% z_data = z_riverBed;
% figure
% map_h = mapshow(z_data,r_dem,'DisplayType','surface');
% % zdatam(map_h,z_data-500)
% hold on
% % mapshow(channelBoundary,'FaceColor','y')
% mapshow(channelBoundaryLine,'Color','y')
% hold off
% axis image
% % hold on;scatter(x_ToBeInterp,y_ToBeInterp,'r*');hold off
% % hold on;scatter(bedPoints(i,1),bedPoints(i,2),'y*'); hold off
%% plot river bed elevation map
figure;
mapshow(z_demNewRiverBed,r_dem,'DisplayType','surface')
axis image
colorbar
figure(1) 
%show result of function BankCrossLinePreprocessriver
%river channel polygon and end points of bank line
mapshow(channelPoly,'FaceColor','none');
hold on
plot(channelPoly.X,channelPoly.Y)
scatter(bLine1(1,1),bLine1(1,2),'r*');
scatter(bLine2(1,1),bLine1(2,2),'r*');
scatter(bLine1(end,1),bLine1(end,2),'bs');
scatter(bLine2(end,1),bLine2(end,2),'bs');

for i=1:length(cLines)
%     if ~isempty(cLines{i})
    xyzdata = cLines{i};
    plot(xyzdata(end,1),xyzdata(end,2),'ko')
    text(xyzdata(1,1),xyzdata(1,2),num2str(i));
%     scatter(xyzdata(1,1),xyzdata(1,2),'go')
%     end
end
hold off
axis image
%% show result of function Channel2Sections
figure(2)
mapshow(sectionPoly,'FaceColor','none')
hold on
for i=1:length(cLines)
    xyzdata = cLines{i};
    plot(xyzdata(end,1),xyzdata(end,2),'bs')
    plot(xyzdata(1,1),xyzdata(1,2),'rs');
end
hold off
axis equal
%% show result of function DiscretizeChannel2Quadrangle
figure(3)
hold on
for i=1:length(sectionPoly)
    subSectionPoly = subSectionPolyCell{i};
    splitPoints1 = splitPoints1Cell{i};
    splitPoints2 = splitPoints2Cell{i};
    mapshow(subSectionPoly,'FaceColor','none')
    plot(sectionPoly(i).X,sectionPoly(i).Y,'g-','LineWidth',2)
    scatter(splitPoints1(:,1),splitPoints1(:,2),'r*')
    scatter(splitPoints2(:,1),splitPoints2(:,2),'b*')
end
clear subSectionPoly splitPoints1 splitPoints2
axis equal
%% show result of function LineInterpolation
figure(4)
hold on
for i=3%:length(sectionPoly)
%     subSectionPoly = subSectionPolyCell{i};
%     crossLines = subSectionLineCell{i};
%     splitPoints1 = splitPoints1Cell{i};
%     splitPoints2 = splitPoints2Cell{i};
%     mapshow(subSectionPoly,'FaceColor','none')
%     plot(sectionPoly(i).X,sectionPoly(i).Y,'g-','LineWidth',2)
%     scatter(splitPoints1(:,1),splitPoints1(:,2),'r*')
%     scatter(splitPoints2(:,1),splitPoints2(:,2),'b*')
%     scatter3(crossLines(:,1),crossLines(:,2),crossLines(:,3));
%     scatter3(crossLines(:,4),crossLines(:,5),crossLines(:,6));
    channelLines = channelLinesCell{i};
    for j = 1:length(channelLines)
        plot3(channelLines(j).X,channelLines(j).Y,channelLines(j).Z,'b-');
    end
end
% clear subSectionPoly splitPoints1 splitPoints2 crossLines
function sectionPoly = Channel2Sections(bLine1,bLine2,cLines)
%Channel2Sections Summary of this function goes here
% split channel polygon to sections based on cross section lines
sectionPoly = struct('X',0,'Y',0,'Geometry','Polygon','BoundingBox',0,'bank1',0,'bank2',0);
for i=1:length(cLines)-1
    cross0 = cLines{i};
    cross1 = cLines{i+1};
    [~,ind0] = PrjToBankLine(cross0(1,1:2),bLine1);
    [~,ind1] = PrjToBankLine(cross1(1,1:2),bLine1);
    bank1 = bLine1(ind0(2):ind1(1),:);
    [~,ind0] = PrjToBankLine(cross0(end,1:2),bLine2);
    [~,ind1] = PrjToBankLine(cross1(end,1:2),bLine2);
    bank2 = bLine2(ind0(2):ind1(1),:);
    sectionPoly(i).bank1 = [cross0(1,1:2);   bank1; cross1(1,1:2)];
    sectionPoly(i).bank2 = [cross0(end,1:2); bank2; cross1(end,1:2)];
    XY =  [bank1;cross1(:,1:2);bank2(end:-1:1,1:2);cross0(end:-1:1,1:2)];
    sectionPoly(i).X = XY(:,1);
    sectionPoly(i).Y = XY(:,2);
    sectionPoly(i).Geometry = 'Polygon';
    sectionPoly(i).BoundingBox = [min(sectionPoly(i).X),min(sectionPoly(i).Y);...
                                  max(sectionPoly(i).X),max(sectionPoly(i).Y) ];
end
end

function line_xy_all = ConnectLines2Polygon(lineB1,lineB2,crossPoints_xyzCell)
% lineB1,lineB2: column vectors represent bank lines
% lineC1,lineC2: column vectors represent bank lines
line_xy_all = struct('lineAll',[],'line1',[],'line2',[],'line3',[],'line4',[]);

for i = 1:length(crossPoints_xyzCell)-1
    lineC1 = crossPoints_xyzCell{i};
    lineC2 = crossPoints_xyzCell{i+1};    
    xy_c1s = lineC1(1,1:2);
    xy_c1e = lineC1(end,1:2);
    xy_c2s = lineC2(1,1:2);
    xy_c2e = lineC2(end,1:2);
    %lineC1 is the first line
    [minV_b1c1e,~] = nearestDistance(lineB1,xy_c1e);
    [minV_b2c1e,~] = nearestDistance(lineB2,xy_c1e);
    if minV_b1c1e < minV_b2c1e %closer to bank1, ban1 is the 2nd line
        secondLine = lineB1;
        fourthLine = lineB2;
    else               %closer to bank2, ban2 is the 2rd line
        secondLine = lineB2;
        fourthLine = lineB1;
    end
    [~,ind_2ndLs] = nearestDistance(secondLine,xy_c1e);
    % distance between the 2nd line section and lineC2 end points
    [minV_2ndLC2s,minI_2ndLC2s] = nearestDistance(secondLine,xy_c2s);
    [minV_2ndLC2e,minI_2ndLC2e] = nearestDistance(secondLine,xy_c2e);
    if minV_2ndLC2s<minV_2ndLC2e %closer to the start point of lineC2
        ind_2ndLe = minI_2ndLC2s;
    else                     %closer to the end point of lineC2
        ind_2ndLe = minI_2ndLC2e;
        lineC2 = lineC2(end:-1:1,:); % change the points' sequence of lineC2        
        xy_c2s = lineC2(1,1:2);
        xy_c2e = lineC2(end,1:2);
    end
    if ind_2ndLs>ind_2ndLe
        signGap = -1;
    else
        signGap = 1;
    end
    
    xy_c1e_Two = [xy_c1e;lineC1(end-1,1:2)];
    xy_intersecLine12 = ExtendToBankLine(xy_c1e_Two,secondLine,ind_2ndLs);
    
    secondLine_section = secondLine(ind_2ndLs:signGap:ind_2ndLe,:);
    if size(secondLine_section,1) == 1
        inFlag = 0;
    else
        inFlag = IsItBetweenTwoPoints(secondLine_section(1,1:2),secondLine_section(2,1:2),...
        xy_intersecLine12);
    end
    if inFlag == 1
        secondLine_section(1,:) =[];
    end
    xy_c2s_Two = [xy_c2s;lineC2(2,1:2)];
    if size(secondLine_section,1) == 1
        inFlag = 0;
    else
        xy_intersecLine23= ExtendToBankLine(xy_c2s_Two,secondLine,ind_2ndLe);
        inFlag = IsItBetweenTwoPoints(secondLine_section(end,1:2),secondLine_section(end-1,1:2),...
        xy_intersecLine23);
    end
    if inFlag == 1
        secondLine_section(end,:) =[];
    end
    secondLine_Final = [xy_intersecLine12;secondLine_section;xy_intersecLine23];
    [~,ind_4thLs] = nearestDistance(fourthLine,xy_c2e); % nearest point to c2e
    [~,ind_4thLe] = nearestDistance(fourthLine,xy_c1s); % nearest point to c1s
    if ind_4thLs>ind_4thLe
        signGap = -1;
    else
        signGap = 1;
    end
    fourthLine_section = fourthLine(ind_4thLs:signGap:ind_4thLe,:);
    xy_c2e_Two = [xy_c2e;lineC2(end-1,1:2)];
    xy_intersecLine34= ExtendToBankLine(xy_c2e_Two,fourthLine,ind_4thLs);
    if size(fourthLine_section,1) == 1
        inFlag = 0;
    else
        inFlag = IsItBetweenTwoPoints(fourthLine_section(1,1:2),fourthLine_section(2,1:2),...
        xy_intersecLine34);
    end
    if inFlag == 1
        fourthLine_section(1,:) =[];
    end
    xy_c1s_Two = [xy_c1s;lineC1(2,1:2)];
    if size(fourthLine_section,1) == 1
        inFlag = 0;
    else
        xy_intersecLine41= ExtendToBankLine(xy_c1s_Two,fourthLine,ind_4thLe);
        inFlag = IsItBetweenTwoPoints(fourthLine_section(end,1:2),fourthLine_section(end-1,1:2),...
        xy_intersecLine41);
    end
    if inFlag == 1
        fourthLine_section(end,:) =[];
    end
    fourthLine_Final = [xy_intersecLine34;fourthLine_section;xy_intersecLine41];
    line_xy_all(i).lineAll = [lineC1(:,1:2);secondLine_Final;lineC2(:,1:2);fourthLine_Final];
    line_xy_all(i).line2 = secondLine_Final;    
    line_xy_all(i).line4 = fourthLine_Final;
    % add end points to lineC1
%     F_scatter = scatteredInterpolant(lineC1,'linear','linear');
    z_intersecLine41 = lineC1(1,3);%F_scatter(xy_intersecLine41(1),xy_intersecLine41(2));
    z_intersecLine12 = lineC1(end,3);%F_scatter(xy_intersecLine12(1),xy_intersecLine12(2));
    xyz_Line12 = [xy_intersecLine12 z_intersecLine12];
    xyz_Line41 = [xy_intersecLine41 z_intersecLine41];
    line_xy_all(i).line1 = [xyz_Line41;lineC1;xyz_Line12];
    % add end points to lineC2
%     F_scatter = scatteredInterpolant(lineC2,'linear','linear');
    z_intersecLine23 = lineC2(1,3);%F_scatter(xy_intersecLine23(1),xy_intersecLine23(2));
    z_intersecLine34 = lineC2(end,3);%F_scatter(xy_intersecLine34(1),xy_intersecLine34(2));
    xyz_Line23 = [xy_intersecLine23 z_intersecLine23];
    xyz_Line34 = [xy_intersecLine34 z_intersecLine34];
    line_xy_all(i).line3 = [xyz_Line23;lineC2;xyz_Line34];
    
end
end

function xy_intersec= ExtendToBankLine(xy_TwoCrossPoints,xyBankLine,nearestInd)
% calculate the the projection of xy_Point to line
% (nearest and second nearest points) 
xy11 = xy_TwoCrossPoints(1,:);
xy12 = xy_TwoCrossPoints(2,:);
xy21 = xyBankLine(nearestInd,:); %point1 on bank line
xyBankLine(nearestInd) = nan;
[~,secondNearestInd] = nearestDistance(xyBankLine,xy11);
xy22 = xyBankLine(secondNearestInd,:); %point2 on bank line
x11 = xy11(1); y11 = xy11(2);
x12 = xy12(1); y12 = xy12(2);
x21 = xy21(1); y21 = xy21(2);
x22 = xy22(1); y22 = xy22(2);
syms x y
eqn1 = (y-y11)*(x11-x12)==(x-x11)*(y11-y12);
eqn2 = (y-y21)*(x21-x22)==(x-x21)*(y21-y22);
[A,B] = equationsToMatrix([eqn1, eqn2], [x, y]);
A = double(A); B = double(B);
sol = linsolve(A,B);
xsol = sol(1);
ysol = sol(2);
xy_intersec = [xsol ysol];
end

function inFlag = IsItBetweenTwoPoints(point1,point2,mypoint)% check whether the the intersection point is between the two bank points
if (mypoint(1)-point1(1))*(mypoint(1)-point2(1))<0&&...
        (mypoint(2)-point1(2))*(mypoint(2)-point2(2))<0
    inFlag = 1;
else
    inFlag = 0;
end
end
% function xy_prj = PrjToBankLine(xy_Point,xyLine,nearestInd)
% % calculate the the projection of xy_Point to line
% % (nearest and second nearest points) 
% xy0 = xy_Point;
% xy1 = xyLine(nearestInd,:);
% xyLine(nearestInd) = nan;
% [~,secondNearestInd] = nearestDistance(xyLine,xy_Point);
% xy2 = xyLine(secondNearestInd,:);
% vec_a = xy0 - xy1; % vector from nearest Point to xyPoint
% vec_b = xy2 - xy1; % vector from nearest Point to second nearest point
% vec_b_e = vec_b/norm(vec_b);
% vec_a_prj = dot(vec_a,vec_b_e)*vec_b_e; %vector projected on line
% xy_prj = vec_a_prj+xy1; %vector to coordinates
% end

function [minV,minInd] = nearestDistance(xyLine,xy_Point)
% distance between line points and single cross point
D = @(xyLine,xy_Point) sum((xyLine-repmat(xy_Point,[length(xyLine) 1])).^2,2);
squareD = D(xyLine,xy_Point);
[minV,minInd] = min(squareD.^0.5);
end
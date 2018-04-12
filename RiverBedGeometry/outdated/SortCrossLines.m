function crossPointsSortedCell= SortCrossLines...
    (crossPoints_xyzCell,bank1Points,bank2Points)
%sort cross lines based on their relative location along the river bed
crossPointsSortedCell = cell(size(crossPoints_xyzCell));

% check whether the points of the two river banks towards to the same direction 
if abs(bank1Points(1,1)-bank2Points(1,1))>abs(bank1Points(1,1)-bank2Points(end,1))
    bank2Points = bank2Points(end:-1:1,:);
end

% define the bottom bank line
meanL_sort = zeros(size(crossPoints_xyzCell));
for i = 1:length(crossPointsSortedCell)
    crossPoints = crossPoints_xyzCell{i};
    twoCrossP1 = crossPoints(1:2,1:2);
    twoCrossP2 = crossPoints(end-1:end,1:2);
    crossP1_xy = twoCrossP1(1,1:2);
    [minV_b1,~] = nearestPointDistance(bank1Points,crossP1_xy);
    [minV_b2,~] = nearestPointDistance(bank2Points,crossP1_xy);
    if minV_b1 < minV_b2 % start point of crossLine is closer to bank1
        bankLineBottom = bank1Points;
        bankLineTop = bank2Points; 
    else % start point of crossLine is closer to bank2
        bankLineBottom = bank1Points;
        bankLineTop = bank2Points; 
    end
    [~,indC1] = nearestPointDistance(bankLineBottom,twoCrossP1(1,1:2));
    [~,indC2] = nearestPointDistance(bankLineTop,twoCrossP2(end,1:2));
    intersecP1 = ExtendToBankLine(twoCrossP1,bankLineBottom,indC1);
    intersecP2 = ExtendToBankLine(twoCrossP2,bankLineTop,indC2);
    [~,P1Ind1] = nearestPointDistance(bankLineBottom,intersecP1);
    bankLineBottom_temp= bankLineBottom; 
    bankLineBottom_temp(P1Ind1,:) = nan;
    [~,P1Ind2] = nearestPointDistance(bankLineBottom_temp,intersecP1);
    [~,P2Ind1] = nearestPointDistance(bankLineTop,intersecP2);
    bankLineTop_temp = bankLineTop;
    bankLineTop_temp(P2Ind1,:) = nan;
    [~,P2Ind2] = nearestPointDistance(bankLineTop_temp,intersecP2);
    p1Inds = sort([P1Ind1,P1Ind2]);
    p2Inds = sort([P2Ind1,P2Ind2]);
    bankLineBottom_add = [bankLineBottom(1:p1Inds(1),:);intersecP1;bankLineBottom(p1Inds(2):end,:)];
    bankLineTop_add = [bankLineTop(1:p2Inds(1),:);intersecP2;bankLineTop(p2Inds(2):end,:)];
    bankLineBottom_lw = RelativeCoordsBankPoints(bankLineBottom_add,0);
    bankLineTop_lw = RelativeCoordsBankPoints(bankLineTop_add,1);
    meanL = (bankLineBottom_lw(p1Inds(1)+1,1)+bankLineTop_lw(p2Inds(1)+1,1))/2;
    meanL_sort(i) = meanL;
end
meanL_sort = [meanL_sort,(1:length(meanL_sort))'];
meanL_sort = sortrows(meanL_sort,1);
crossPointsSortedCell = crossPoints_xyzCell(meanL_sort(:,2));
end

function xy_intersec= ExtendToBankLine(xy_TwoCrossPoints,xyBankLine,nearestInd)
% calculate the the intersection of the line defined by xy_TwoCrossPoints
% to line defined by two bank points(nearest and second nearest points) 
xy11 = xy_TwoCrossPoints(1,:);
xy12 = xy_TwoCrossPoints(2,:);
xy21 = xyBankLine(nearestInd,:); %point1 on bank line
xyBankLine(nearestInd) = nan;
[~,secondNearestInd] = nearestPointDistance(xyBankLine,xy11);
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

function [minV,minInd] = nearestPointDistance(xyLine,xy_Point)
% distance between line points and single cross point
D = @(xyLine,xy_Point) sum((xyLine-repmat(xy_Point,[length(xyLine) 1])).^2,2);
squareD = D(xyLine,xy_Point);
[minV,minInd] = min(squareD.^0.5);
end
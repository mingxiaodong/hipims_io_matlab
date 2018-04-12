function xy_intersec= ExtendToBankLine(xy_TwoCrossPoints,xyBankLine,varargin)
% calculate the the intersection of the line defined by xy_TwoCrossPoints
% to line defined by two bank points(nearest and second nearest points) 
xy11 = xy_TwoCrossPoints(1,:); 
xy12 = xy_TwoCrossPoints(2,:);
[~, nearestInd] = nearestPointDistance(xyBankLine,xy11);
if ~isempty(varargin)
    nearestInd = varargin{1};
end
xy21 = xyBankLine(nearestInd,:); %point1 on bank line
xyBankPoints = xyBankLine;
xyBankPoints(nearestInd) = nan;
betweenFlag = false; % whether the intersection points is between point1 and point2 on bank line
while(~betweenFlag)
    [~,secondNearestInd] = nearestPointDistance(xyBankPoints,xy11);
    xy22 = xyBankLine(secondNearestInd,:); %point2 on bank line
    xyBankPoints(secondNearestInd) = nan;
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
    betweenFlag = (xsol-x21)*(xsol-x22)<=0&(ysol-y21)*(ysol-y22)<=0;
end
xy_intersec = [xsol ysol];
end
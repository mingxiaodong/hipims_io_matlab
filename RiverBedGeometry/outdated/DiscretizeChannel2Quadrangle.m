function [splitPointsAll_1,splitPointsAll_2,sectionPoly] = DiscretizeChannel2Quadrangle(bank1Points,bank2Points,errorDistance)
%DiscretizeChannel2Quadrangle Split the river channel to quadrangular or
%triangular polygon according to a given error threshold.
%   bank1Points and bank2Points are the coordiantes of points in each bank,
%   the direction of the bank line should be the same. errorDistance is the
%   maximum distance between bank points and the discretized bank line,
%   which indicates the error of discretization. splitPointsAll_1 and 
%   splitPointsAll_2 are the split points in each bank. sectionPoly is the
%   struct containing information of the divided polygons of river channel.
%   Created by Xiaodong Ming on 2017-10-30
%   See also DiscretizeChannel2Points
splitInd1 = SplitSingleBank(bank1Points,errorDistance); splitInd1 = [1 splitInd1]';
splitInd2 = SplitSingleBank(bank2Points,errorDistance); splitInd2 = [1 splitInd2]';
split_Points1 = bank1Points(splitInd1,:);
split_Points2 = bank2Points(splitInd2,:);
if numel(splitInd1)==2&&numel(splitInd2)==2
    splitPointsAll_1 = split_Points1;
    splitPointsAll_2 = split_Points2;
else
    [xy_2prj1, indBetween1] = PrjToBankLine(split_Points2,bank1Points);
    [xy_1prj2, indBetween2] = PrjToBankLine(split_Points1,bank2Points);
    splitPointsAll_1 = CombineSplitPointsOfTwoBank(bank1Points,splitInd1,xy_2prj1,indBetween1);
    splitPointsAll_2 = CombineSplitPointsOfTwoBank(bank2Points,splitInd2,xy_1prj2,indBetween2);
end
sectionPoly = struct('X',0,'Y',0,'Geometry','Polygon','BoundingBox',0,'CutPoints0',0,'CutPoints1',0);
for n=1:length(splitPointsAll_1)-1
    x1 = splitPointsAll_1(n,1);   y1 = splitPointsAll_1(n,2);
    x2 = splitPointsAll_1(n+1,1); y2 = splitPointsAll_1(n+1,2);
    x3 = splitPointsAll_2(n+1,1); y3 = splitPointsAll_2(n+1,2);
    x4 = splitPointsAll_2(n,1);   y4 = splitPointsAll_2(n,2);
    sectionPoly(n).X = [x1 x2 x3 x4]';
    sectionPoly(n).Y = [y1 y2 y3 y4]';
    sectionPoly(n).CutPoints0 = [x1 y1; x4 y4];
    sectionPoly(n).CutPoints1 = [x2 y2; x3 y3];
    sectionPoly(n).Geometry = 'Polygon';
    sectionPoly(n).BoundingBox = [min(sectionPoly(n).X),min(sectionPoly(n).Y);...
                                  max(sectionPoly(n).X),max(sectionPoly(n).Y)];
end

end

function splitPointsAll = CombineSplitPointsOfTwoBank(bankPoints,splitIndOrig,prjPoints,prjIndBetween)
% prjPoints = xy_1prj2;
% splitIndOrig = splitInd2;
% bankPoints = bank2Points;
% indBetween = indBetween2;
indBetween = prjIndBetween;
split_Points_all = [prjPoints;bankPoints]*0;
splitIndPrj = zeros(length(indBetween),1);
n_ind = 1;
for i=1:length(indBetween)
    row_e = indBetween(i,1);
    if i==1
        row_s = 1;
    else
        row_s = indBetween(i-1,2);
    end
    split_Points_all(n_ind:row_e-row_s+n_ind,:) = bankPoints(row_s:row_e,:);
    n_ind = row_e-row_s+n_ind+1;
    split_Points_all(n_ind,:) = prjPoints(i,:);
    splitIndPrj(i) = n_ind;
    splitIndOrig(splitIndOrig>=n_ind) = splitIndOrig(splitIndOrig>=n_ind)+1;
    n_ind = n_ind+1;
    if i==length(indBetween)
        split_Points_all(n_ind:end,:) = bankPoints(row_e+1:length(bankPoints),:);
    end
end
splitIndAll = sort([splitIndPrj; splitIndOrig]);
splitPointsAll = split_Points_all(splitIndAll,:);
end

function outputArg1 = SplitSingleBank(bankPoints_l,errorDistance)
%SplitRiverChannel Summary of this function goes here
%   Detailed explanation goes here
% bankPoints_l: relative length of bank points from the first point.
splitInd = length(bankPoints_l);
n_section = 1;
indPs = 1; %point index at start of section
indPe = length(bankPoints_l); %point index at end of section
while(indPs<indPe)
    x1 = bankPoints_l(indPs,1); y1 = bankPoints_l(indPs,2);
    x2 = bankPoints_l(indPe,1); y2 = bankPoints_l(indPe,2);
    A = y2-y1;
    B = x1-x2;
    C = (y2-y1)*(-x1)+(x1-x2)*(-y1);
    % lineEqns = (y2-y1)*(x-x1)+(x1-x2)*(y-y1)==0;
    distanceEqn = @(x,y) abs(A*x+B*y+C)/sqrt(A^2+B^2);
    distanceValues = distanceEqn(bankPoints_l(indPs:indPe,1),bankPoints_l(indPs:indPe,2));
    [M,I] = max(distanceValues);
    if M>errorDistance
        indPe = I+indPs-1;
    else
        splitInd(n_section) = indPe;
        n_section = n_section+1;
        indPs = indPe+1;
        indPe = length(bankPoints_l);
    end      
end
outputArg1 = splitInd;
end

function [xy_prj, indBetween]= PrjToBankLine(xy_Point,xyLine)
% calculate the the projection of xy_Point to line
% (nearest and second nearest points) 
xy_prj = xy_Point*0;
indBetween = xy_prj;
for n=1:size(xy_prj,1)
    xy0 = xy_Point(n,1:2);
    xyBankLine = xyLine;
    [~,nearestInd] = nearestPointDistance(xyBankLine,xy0);
    xy_1near = xyBankLine(nearestInd,:); % the nearest point to xy0
    if nearestInd==1 % the previous point of xy_1near
        xy_1near_pre = xy_1near;
    else
        xy_1near_pre = xyBankLine(nearestInd-1,:);
    end
    if nearestInd==length(xyBankLine) % the next point of xy_1near
        xy_1near_nxt = xy_1near;
    else
        xy_1near_nxt = xyBankLine(nearestInd+1,:);
    end
%     xyBankLine(nearestInd,:) = nan;
    dotProduct_pre = dot(xy0-xy_1near,xy_1near_pre - xy_1near);
    dotProduct_nxt = dot(xy0-xy_1near,xy_1near_nxt - xy_1near);
    if max(dotProduct_pre,dotProduct_nxt)<=0
        xy_prj0 = xy_1near;
    else
        if dotProduct_pre>=dotProduct_nxt
            xy_2near = xy_1near_pre;
        else
            xy_2near = xy_1near_nxt;
        end
        vec_a = xy0 - xy_1near; % vector from xy_1near to xy0
        vec_b = xy_2near - xy_1near; % vector from xy_1near to xy_2near
        if norm(vec_b)==0
            xy_prj0 = xy_1near;
        else
            vec_b_e = vec_b/norm(vec_b);
            vec_a_prj = dot(vec_a,vec_b_e)*vec_b_e; %vector projected on line
            xy_prj0 = vec_a_prj+xy_1near; %vector to coordinates
        end
    end

    % find the location of xy_prj in xyLine
    [~,ind1] = nearestPointDistance(xyBankLine,xy_prj0);
    xyBankLine(ind1,:) = nan;
    [~,ind2] = nearestPointDistance(xyBankLine,xy_prj0);
    indBetween0 = [ind1,ind2];
    xy_prj(n,1:2) = xy_prj0;
    indBetween(n,1:2) = sort(indBetween0);
end
end
function [minV,minInd] = nearestPointDistance(xyLine,xy_Point)
% distance between line points and single cross point
D = @(xyLine,xy_Point) sum((xyLine-repmat(xy_Point,[length(xyLine) 1])).^2,2);
squareD = D(xyLine,xy_Point);
[minV,minInd] = min(squareD.^0.5);
end

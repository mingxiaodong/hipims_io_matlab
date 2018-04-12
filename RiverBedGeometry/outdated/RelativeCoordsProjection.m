function [output_l,output_w] = RelativeCoordsProjection(bankLineTop,...
    bankLineBottom,inputPoints)
% find the nearest two points for input each point
output_l = zeros(length(inputPoints(:,1)),1);
output_w = output_l;
for i=1:length(inputPoints(:,1))
    x0 = inputPoints(i,1); y0 = inputPoints(i,2);
    [xy_prj1, indBetween1] = PrjToBankLine([x0,y0],bankLineBottom);
    [xy_prj2, indBetween2] = PrjToBankLine([x0,y0],bankLineTop);
    indBetween1 = sort(indBetween1);
    indBetween2 = sort(indBetween2);
    bankLineBottom_add = [bankLineBottom(1:indBetween1(1),:);xy_prj1;bankLineBottom(indBetween1(2):end,:)];
    bankLineTop_add = [bankLineTop(1:indBetween2(1),:);xy_prj2;bankLineTop(indBetween2(2):end,:)];
    bankLineBottom_lw = RelativeCoordsBankPoints(bankLineBottom_add,0);
    bankLineTop_lw = RelativeCoordsBankPoints(bankLineTop_add,1);
    meanL = (bankLineBottom_lw(indBetween1(1)+1,1)+bankLineTop_lw(indBetween2(1)+1,1))/2;
    output_l(i) = meanL;
    distance1 = norm(xy_prj1-[x0,y0]);
    distance2 = norm(xy_prj2-[x0,y0]);
    w_aboveBottom = distance1/(distance1+distance2);
    output_w(i) = w_aboveBottom;
end
end


function [xy_prj, indBetween]= PrjToBankLine(xy_Point,xyLine)
% calculate the the projection of xy_Point to line
% (nearest and second nearest points) 
xy0 = xy_Point;
[~,nearestInd] = nearestDistance(xyLine,xy0);
xy1 = xyLine(nearestInd,:);
xyLine(nearestInd,:) = nan;
[~,secondNearestInd] = nearestDistance(xyLine,xy_Point);
xy2 = xyLine(secondNearestInd,:);
xyLine(nearestInd,:) = xy1;
vec_a = xy0 - xy1; % vector from nearest Point to xyPoint
vec_b = xy2 - xy1; % vector from nearest Point to second nearest point
vec_b_e = vec_b/norm(vec_b);
vec_a_prj = dot(vec_a,vec_b_e)*vec_b_e; %vector projected on line
xy_prj = vec_a_prj+xy1; %vector to coordinates
% find the location of xy_prj in xyLine
[~,ind1] = nearestDistance(xyLine,xy_prj);
xyLine(ind1,:) = nan;
[~,ind2] = nearestDistance(xyLine,xy_prj);
indBetween = [ind1,ind2];
end

function [minV,minInd] = nearestDistance(xyLine,xy_Point)
% distance between line points and single cross point
D = @(xyLine,xy_Point) sum((xyLine-repmat(xy_Point,[length(xyLine) 1])).^2,2);
squareD = D(xyLine,xy_Point);
[minV,minInd] = min(squareD.^0.5);
end
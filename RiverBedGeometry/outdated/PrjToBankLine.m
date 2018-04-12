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
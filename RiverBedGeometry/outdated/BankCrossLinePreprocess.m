function [bLine1,bLine2,cLines,channelPoly] = BankCrossLinePreprocess(bankLine1, bankLine2, crossLines)
%BankCrossLinePreprocess
%   1. unify the direction of points on the two river banks lines
%   2. unify the direction of points on each cross line from bank1 to bank2
%   3. sort	the cross-section lines towards the same direction of bank lines
%   4. cut part of cross-section lines outside channel
%% 1.unify the direction of points on the two river banks lines
if isnan(bankLine1(end,1))
    bankLine1(end,:)=[];
end
if isnan(bankLine2(end,1))
    bankLine2(end,:)=[];
end
if abs(bankLine1(1,1)-bankLine2(1,1))>abs(bankLine1(1,1)-bankLine2(end,1))
    bankLine2 = bankLine2(end:-1:1,:);
end
bLine1 = bankLine1;
bLine2 = bankLine2;
channelPoly = struct('X',0,'Y',0,'Geometry','Polygon','BoundingBox',0);
channelPoly.X = [bankLine1(:,1);bankLine2(end:-1:1,1);nan];
channelPoly.Y = [bankLine1(:,2);bankLine2(end:-1:1,2);nan];
channelPoly.BoundingBox = [min(channelPoly.X),min(channelPoly.Y);...
    max(channelPoly.X),max(channelPoly.Y) ];
%% 2.unify the direction of points on each cross line from bank1 to bank2
for i=1:length(crossLines)
    xyzdata = crossLines{i};
    xy_Point = xyzdata(1,1:2);
    [xy_prj1, ~] = PrjToBankLine(xy_Point,bankLine1);
    [xy_prj2, ~] = PrjToBankLine(xy_Point,bankLine1);
    if norm(xy_Point-xy_prj1)>norm(xy_Point-xy_prj2)
        xyzdata = xyzdata(end:-1:1,:);
    end
    crossLines{i} = xyzdata;
end
%% 3 sort the cross-section lines towards the same direction of bank lines
%% 4 cut the cross-section lines part outside channel polygon
cLines = crossLines;
bankLine1_infill = InfillPoints(bankLine1,4);
bankLine2_infill = InfillPoints(bankLine2,4);
% firstCrossPoints = nan(length(crossLines),2);
prjPoints = nan(length(crossLines),2);
prjInd = prjPoints;
for i=1:length(crossLines)
    xyzdata = crossLines{i};
    % 1. find the position point of each cross line inside channel
    in = inpolygon(xyzdata(:,1),xyzdata(:,2),channelPoly.X,channelPoly.Y);
    if sum(in)==0 % cross-section line is outside the channel polygon
        if norm(xyzdata(1,1:2)-bankLine1(1,1:2))<norm(xyzdata(1,1:2)-bankLine1(end,1:2))%the first
            [xy0,ind0] = PrjToBankLine(bankLine1(1,1:2),xyzdata(:,1:2));
            [xy1,ind1] = PrjToBankLine(bankLine2(1,1:2),xyzdata(:,1:2));
            prjInd(i,:) = [1 1];
        else %the last cross section line
            [xy0,ind0] = PrjToBankLine(bankLine1(end,1:2),xyzdata(:,1:2));
            [xy1,ind1] = PrjToBankLine(bankLine2(end,1:2),xyzdata(:,1:2));
            prjInd(i,:) = [numel(bankLine1_infill) numel(bankLine1_infill)];
        end
        prjPoints(i,:) = xy0;
        xyzdata_cut = xyzdata(ind0(2):ind1(1),:);
    else
        xyzdata_cut = xyzdata(in,:);
        firstCrossPoint = xyzdata_cut(1,1:2);
        lastCrossPoint = xyzdata_cut(end,1:2);
    % 2. project the main position points to bank line1
        [xy_prj0, indbetween] = PrjToBankLine(firstCrossPoint,bankLine1_infill);
        [xy_prj1, ~] = PrjToBankLine(lastCrossPoint,bankLine2_infill);
        prjPoints(i,:) = xy_prj0;
        prjInd(i,:) = indbetween;
        [xy0,ind0] = PrjToBankLine(xy_prj0,xyzdata(:,1:2));
        [xy1,ind1] = PrjToBankLine(xy_prj1, xyzdata(:,1:2));
    end
    
    z0 = mean(xyzdata(ind0,3));
    z1 = mean(xyzdata(ind1,3));
    xyzdata_exp = [xy0,z0;xyzdata_cut;xy1,z1];
    cLines{i} = xyzdata_exp;
end
    % 3. compare the position of the projected points
meanPrjInd = mean(prjInd,2);
[~,sortInd] = sort(meanPrjInd);
cLines = cLines(sortInd);
end


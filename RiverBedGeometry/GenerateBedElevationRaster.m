function [Z,R] = GenerateBedElevationRaster(bankLine1,bankLine2,cellsize,varargin)
% GENERATEBEDELEVATIONRASTER
% [Z,R] = GenerateBedElevationRaster(bankLine1,bankLine2,cellsize,crossSectionLine,errorDistance)
% bankLine1,bankLine2: xy coordinates of points on the two banks of river
% cellsize: is the output raster cell size
% crossSectionLine: a cell of xyz value of points on cross section lines. 
%       If it is empty, then default two parabola lines will be given.
% errorDistance: the error value when linearize the bank line
% Created by Xiaodong Ming on 2017-3-11
if isempty(varargin)
    errorDistance = cellsize*2;
    crossSectionLine = [];
elseif length(varargin)==1
    if isnumeric(varargin{1})
        errorDistance = varargin{1};
        crossSectionLine = [];
    else
        crossSectionLine = varargin{1};
        errorDistance = cellsize*2;
    end
elseif length(varargin)==2
    if isnumeric(varargin{1})
        errorDistance = varargin{1};
        crossSectionLine = varargin{2};
    else
        errorDistance = varargin{2};
        crossSectionLine = varargin{1};
    end
end
%% BankCrossLinePreprocess
% 1	Points on the two bank lines should be listed towards the same direction (upstream or downstream)
% 2	Points on each cross-section line should be listed towards the same direction (bank1 to bank2)
% 3	Cross-section lines should be listed towards the same direction of bank lines (upstream or downstream)
%*****************subfunction: BankCrossLinePreprocess******************
[bLine1, bLine2,cLines,channelPoly] = BankCrossLinePreprocess(bankLine1, bankLine2, crossSectionLine);
%% River channel segmentation
%*****************subfunction: Channel2Sections**********************   
sectionPoly = Channel2Sections(bLine1,bLine2,cLines);
%% Linearization of channel sections (for each section)
splitPoints1Cell = cell(length(sectionPoly),1);
splitPoints2Cell = cell(length(sectionPoly),1);
subSectionPolyCell = cell(length(sectionPoly),1);
for i = 1:length(sectionPoly)
    bank1Points = sectionPoly(i).bank1;
    bank2Points = sectionPoly(i).bank2;
%*****************subfunction: DiscretizeChannel2Quadrangle**************** 
    [sP1,sP2,allSubSections] = DiscretizeChannel2Quadrangle(bank1Points,bank2Points,errorDistance);
    splitPoints1Cell{i} = sP1;
    splitPoints2Cell{i} = sP2;
    subSectionPolyCell{i} = allSubSections;
end
%% Discretization of channel subsections (for each subsection)
interpGap = (0:0.01:1)';
%*****************subfunction: LineInterpolation***************************
[~, cutLineCell] = LineInterpolation(cLines,subSectionPolyCell,splitPoints1Cell,splitPoints2Cell,interpGap);
resolution = cellsize;
%% perform interpolation for each riverPolySection generated based on cross section lines
wholeRaster_r = makerefmat(min(channelPoly.X),max(channelPoly.Y),cellsize,-cellsize);
wholeCol = round(range(channelPoly.X)/cellsize)+1;
wholeRow = round(range(channelPoly.Y)/cellsize)+1;
wholeRaster_z = nan(wholeRow,wholeCol);
%*****************subfunction: Raster2FeaturePoints************************
[wholeRaster_x,wholeRaster_y] = Raster2FeaturePoints(wholeRaster_z,wholeRaster_r);
for i = 1:length(subSectionPolyCell)
    cutLines = cutLineCell{i};
    allSubSections   = subSectionPolyCell{i};
    for j=1:length(allSubSections)
        crossLine0 = cutLines(j).XYZ;
        crossLine1 = cutLines(j+1).XYZ;
        subSecPoly1 = allSubSections(j);
        in_subsection = inpolygon(wholeRaster_x,wholeRaster_y,subSecPoly1.X,subSecPoly1.Y);
        if sum(in_subsection(:))>0
%*****************subfunction: DiscretizeChannel2Points********************
            [X_grid,Y_grid,Z_grid] = DiscretizeChannel2Points(crossLine0,crossLine1,resolution);
            X_grid = X_grid(:); Y_grid = Y_grid(:); p = [X_grid,Y_grid];
            [p,ia,~]= unique(p,'rows','stable');
            v = Z_grid(ia);
            F = scatteredInterpolant(p,v,'linear','nearest');
            xq = wholeRaster_x(in_subsection);
            yq = wholeRaster_y(in_subsection);
            zq = F(xq,yq);
            wholeRaster_z(in_subsection) = zq;
        end
    end
end
%%
Z = wholeRaster_z;
R = wholeRaster_r;
end

%% ************************subfunctions**************************************
% with subfunctions: PrjToBankLine, InfillPoints
function [bLine1,bLine2,cLines,channelPoly] = BankCrossLinePreprocess(bankLine1, bankLine2, crossLines)
%BankCrossLinePreprocess
%   1. unify the direction of points on the two river banks lines
%   2. unify the direction of points on each cross line from bank1 to bank2
%   3. sort	the cross-section lines towards the same direction of bank lines
%   4. cut part of cross-section lines outside channel
%%1.unify the direction of points on the two river banks lines
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
%%2.unify the direction of points on each cross line from bank1 to bank2
if isempty(crossLines)
    x1 = linspace(bLine1(1,1),bLine2(1,1),10);
    y1 = linspace(bLine1(1,2),bLine2(1,2),10);
    z1 = 40*(linspace(0,1,10)-0.5).^2;
    x2 = linspace(bLine1(end,1),bLine2(end,1),10);
    y2 = linspace(bLine1(end,2),bLine2(end,2),10);
    z2 = 40*(linspace(0,1,10)-0.5).^2;
    crossLines{1} = [x1',y1',z1']; clear x1 y1 z1
    crossLines{2} = [x2',y2',z2']; clear x2 y2 z2
end

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
%%3 sort the cross-section lines towards the same direction of bank lines
%%4 cut the cross-section lines part outside channel polygon
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

% with subfunctions: PrjToBankLine
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

% with subfunctions: SplitSingleBank, PrjToBankLine, CombineSplitPointsOfTwoBank
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
% with subfunctions: RelativeCoordsCrossPoints
function [X_grid,Y_grid,Z_grid] = DiscretizeChannel2Points(crossLine0,crossLine1,resolution)
%DiscretizeChannel2Points Discretize quadrangular/triangular channel
%section to feature points based on the elevation of the points in the 
%start and end cross lines of the river channel section. The banks line of
%this channel section must be straight!
%   crossLine1,crossLine2: x,y,z coordinates and z value of cross line points
%   resolution: resolution of the output points
%   Created by Xiaodong Ming on 2017-10-30.
%   See also DiscretizeChannel2Quadrangle
    DistanceTwoPoints = @(xy1,xy2) (sum((xy1-xy2).^2)).^0.5;
    crossLine0_D = DistanceTwoPoints(crossLine0(1,1:2),crossLine0(end,1:2));
    crossLine1_D = DistanceTwoPoints(crossLine1(1,1:2),crossLine1(end,1:2));
    bankLine0_D = DistanceTwoPoints(crossLine0(1,1:2),crossLine1(1,1:2));
    bankLine1_D = DistanceTwoPoints(crossLine0(end,1:2),crossLine1(end,1:2));
    n_bank  = round(max( bankLine0_D, bankLine1_D)/resolution)+1;
    n_cross = round(max(crossLine0_D,crossLine1_D)/resolution)+1;
    W = linspace(0,1, n_cross)';
    L = linspace(0,1,  n_bank)';
    X_grid = nan(n_cross,n_bank);
    Y_grid = X_grid;
    Z_grid = X_grid;
    crossLines0_W = RelativeCoordsCrossPoints(crossLine0,0);
    crossLines1_W = RelativeCoordsCrossPoints(crossLine0,1);
    crossLine0_mW = interp1(crossLines0_W(:,2),crossLine0(:,3),W);
    crossLine1_mW = interp1(crossLines1_W(:,2),crossLine1(:,3),W);
    lineraInterp = @(xyz0,xyz1,r) xyz0+r*(xyz1-xyz0);
    crossLine0_XYZ = lineraInterp(crossLine0(  1,:),crossLine0(end,:),W);
    crossLine1_XYZ = lineraInterp(crossLine1(  1,:),crossLine1(end,:),W);
    crossLine0_XYZ(:,3) = crossLine0_mW;
    crossLine1_XYZ(:,3) = crossLine1_mW;

    for i=1:n_cross
        bankLineXYZ  = lineraInterp(crossLine0_XYZ(i,:),crossLine1_XYZ(i,:),L);
        X_grid(i,:) = bankLineXYZ(:,1)';
        Y_grid(i,:) = bankLineXYZ(:,2)';
        Z_grid(i,:) = bankLineXYZ(:,3)';
    end

end
% with subfunctions: nearestPointDistance
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

% with subfunctions: RelativeCoordsCrossPoints, RelativeCoordsBankPoints
function [channelLinesCell, cutLineCell]= LineInterpolation(cLines,subSectionPolyCell,splitPoints1Cell,splitPoints2Cell,interpGap)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
subSectionLineCell = cell(length(subSectionPolyCell),1);
cutLineCell = cell(length(subSectionPolyCell),1);
channelLinesCell = cell(length(subSectionPolyCell),1);

xyInterp = @(xy0,xy1,interpRate) [linspace(xy0(1),xy1(1),numel(interpRate));...
                               linspace(xy0(2),xy1(2),numel(interpRate))];
for i=1:length(subSectionPolyCell)
    crossLine0 = cLines{i};
    xy0 = crossLine0(1,:); xy1 = crossLine0(end,:);
    crossLine0_interp_XY = xyInterp(xy0,xy1,interpGap);
    crossLine1 = cLines{i+1};
    xy0 = crossLine1(1,:); xy1 = crossLine1(end,:);
    crossLine1_interp_XY = xyInterp(xy0,xy1,interpGap);
    
    crossLine0_W = RelativeCoordsCrossPoints(crossLine0(:,1:2),0);
    xv0 = [crossLine0_W(:,2),crossLine0(:,3)];
    [~,ia,~] = unique(xv0(:,1),'stable'); xv0 = xv0(ia,:);
    crossLine1_W = RelativeCoordsCrossPoints(crossLine1(:,1:2),1);
    xv1 = [crossLine1_W(:,2),crossLine1(:,3)];
    [~,ia,~] = unique(xv1(:,1),'stable'); xv1 = xv1(ia,:);
    crossLine0_interp_Z = interp1(xv0(:,1),xv0(:,2),interpGap);
    crossLine1_interp_Z = interp1(xv1(:,1),xv1(:,2),interpGap);
    crossLine0_interp = [crossLine0_interp_XY',crossLine0_interp_Z];
    crossLine1_interp = [crossLine1_interp_XY',crossLine1_interp_Z];
%     sSPoly = subSectionPolyCell{i};
    sPts1 = splitPoints1Cell{i};
    sPts2 = splitPoints2Cell{i};
    sPtsMd = (sPts1+sPts2)/2;
    sPtsMd_L = RelativeCoordsBankPoints(sPtsMd,0);sPtsMd_L = sPtsMd_L(:,1);
    cutLineStc = struct('X',0,'Y',0,'Z',0);
    for j=1:length(sPtsMd)
        xy0 = sPts1(j,:); xy1 = sPts2(j,:);
        cutLine_XY = xyInterp(xy0,xy1,interpGap);
        k = sPtsMd_L(j);
        cutLine_Z = (1-k)*crossLine0_interp_Z+k*crossLine1_interp_Z;
        cutLine = [cutLine_XY',cutLine_Z];
        cutLineStc(j).X = cutLine(:,1);
        cutLineStc(j).Y = cutLine(:,2);
        cutLineStc(j).Z = cutLine(:,3);
        cutLineStc(j).XYZ = cutLine;
    end
    channelLines = struct('X',0,'Y',0,'Z',0);
    Xall = [cutLineStc.X];
    Yall = [cutLineStc.Y];
    Zall = [cutLineStc.Z];
    for n=1:length(interpGap)
        channelLines(n).X = Xall(n:length(interpGap):end)';
        channelLines(n).Y = Yall(n:length(interpGap):end)';
        channelLines(n).Z = Zall(n:length(interpGap):end)';
    end
    subSectionLineCell{i} = [crossLine0_interp,crossLine1_interp];
    cutLineCell{i} = cutLineStc;
    channelLinesCell{i}  = channelLines;
end
end

% without subfunctions
function [minV,minInd] = nearestPointDistance(xyLine,xy_Point)
% distance between line points and single cross point
D = @(xyLine,xy_Point) sum((xyLine-repmat(xy_Point,[length(xyLine) 1])).^2,2);
squareD = D(xyLine,xy_Point);
[minV,minInd] = min(squareD.^0.5);
end
% without subfunctions
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
% without subfunctions
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
% without subfunctions
function bank1Points_lw = RelativeCoordsBankPoints(bank1Points,w_Value)
%%bank1Points_lw = RelativeCoordsBankPoints(bank1Points,w_Value)
% build the ralative coordinates for river bank points
% L coordinate: length along the river bank (ratio)
% W coordinate: width across the river bank (ratio)
% H coordinate: height of the river bank    (absolute)
    % function to calculate the distance between two neigbour points in one river bank
Distance_InterPoints = @(A)((A(1:end-1,1)-A(2:end,1)).^2+(A(1:end-1,2)-A(2:end,2)).^2).^0.5;
%%*****for bank 1
interDis_Bank1Points = Distance_InterPoints(bank1Points);
cumuDist_Bank1Points = cumsum(interDis_Bank1Points);
% L, W, and H coordinate of river bank 1
coor_L_Bank1Points = [0; cumuDist_Bank1Points/sum(interDis_Bank1Points)];
coor_W_Bank1Points = coor_L_Bank1Points*0+w_Value;
bank1Points_lw = [coor_L_Bank1Points,coor_W_Bank1Points];
end
% without subfunctions
function cross1Points_lw = RelativeCoordsCrossPoints(cross1Points,l_Value)
%%build the ralative coordinates for river bank points
% L coordinate: length along the river bank (ratio)
% W coordinate: width across the river bank (ratio)
    % function to calculate the distance between two neigbour points in one river bank
cross1Points = cross1Points(:,1:2);
Distance_InterPoints = @(A)((A(1:end-1,1)-A(2:end,1)).^2+(A(1:end-1,2)-A(2:end,2)).^2).^0.5;
%%*****for bank 1
interDis_Cross1Points = Distance_InterPoints(cross1Points);
cumuDist_Cross1Points = cumsum(interDis_Cross1Points);
% L, W, and H coordinate of river bank 1
coor_W_Cross1Points = [0; cumuDist_Cross1Points/sum(interDis_Cross1Points)];
coor_L_Cross1Points = coor_W_Cross1Points*0+l_Value;
cross1Points_lw = [coor_L_Cross1Points,coor_W_Cross1Points];
end
% without subfunctions
function [X,Y] = Raster2FeaturePoints(Z,R)
%[X,Y] = RASTER2FEATUREPOINTS(Z,R) Conver raster grid to feature points
%
%   [X,Y] = Raster2FeaturePoints(Z,R) rerurns the X and Y coordinates of
%   the feature points. Z is the raster grid matrix and R is the reference
%   matrix.
%   Created by Xiaodong Ming on 2017-08-02
%   See also MAKEREFMAT, ARCGRIDREAD.
dx = R(2,1); dy = R(1,2);
[m,n] = size(Z);
x11 = R(3,1)+dx;
y11 = R(3,2)+dy;
x = x11:dx:x11+dx*(n-1);
y = y11:dy:y11+dy*(m-1);
[X,Y] = meshgrid(x,y);
end
% without subfunctions
function outPutPoints = InfillPoints(inputPoints,infillTimes)
% linear infill points(centre)
if nargin == 1
    infillTimes=1;
end
for i = 1:infillTimes
    inputPoints_half = (inputPoints(1:end-1,:)+inputPoints(2:end,:))/2;
    newX = [inputPoints(1:end-1,1),inputPoints_half(:,1)];
    newY = [inputPoints(1:end-1,2),inputPoints_half(:,2)];
    newX = newX'; newX = newX(:);
    newY = newY'; newY = newY(:);
    outPutPoints = [[newX,newY];inputPoints(end,:)];
    inputPoints = outPutPoints;
end
end
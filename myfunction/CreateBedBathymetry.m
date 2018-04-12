function [bedRaster,bedPoints]= CreateBedBathymetry(bedPoly_xy,bankLine1,bankLine2,crossLines,varargin) 
%CrossSection2Bathymetry Generate river bed bathymetry based on its
% 3D river cross section lines, river bound polygon and bank lines
%
%   bedPoints = CrossSection2Bathymetry(bedPoly_xy,
%   bankLine1,bankLine2,crossLines,riverBedPoints_xy)
%   Generate a vector giving the elevation of riverBedPoints_xy. 
%   bedPoly_xy is a 3-column matrix representing river bound polygon.
%   bankLine1 and bankLine2 are 2-column matrice representing river 
%   bound lines. crossLines is a cell consisting of 3-column 
%   matrices giving the x,y,z coordinates of 3D points from the cross 
%   section lines
%   See also ConnectLines2Polygon, SortCrossLines

% Created by Xiaodong Ming on 24 Jun 2017.
%% check the format of inputs
if size(bedPoly_xy,2)~=2
    error('Format of input channelBoundPoly should be a 2-column matrix')
end
if size(bankLine1,2)~=2||size(bankLine2,2)~=2
    error('Format of input riveBankLine should be a 2-column matrix')
end
if ~iscell(crossLines)
    error('Format of input crossSection3DLines should be a cell')
end
if size(varargin{1},1)==3&&size(varargin{1},2)==2 % a reference matrix
    bedRaster_r = varargin{1};
    nrows = range(bedPoly_xy(:,2))/bedRaster_r(2); nrows = ceil(nrows);
    ncols = range(bedPoly_xy(:,1))/bedRaster_r(2); ncols = ceil(ncols);
    bedRaster_z = nan(nrows,ncols);
    [x_bedPoints,y_bedPoints] = Raster2FeaturePoints(bedRaster_z,bedRaster_r);
    ind_raster = inpolygon(x_bedPoints,y_bedPoints,bedPoly_xy(:,1),bedPoly_xy(:,2));
    xy_bedPoints = [x_bedPoints(ind_raster),y_bedPoints(ind_raster)];
else
    xy_bedPoints = varargin{1};
    bedRaster_z = [];
    bedRaster_r = [];
end
    
%polyarea
%% check the channelBoundaryLine data
bank1Points = bankLine1;
bank2Points = bankLine2;

% check whether the points of the two river banks towards to the same direction 
if abs(bank1Points(1,1)-bank2Points(1,1))>abs(bank1Points(1,1)-bank2Points(end,1))
    bank2Points = bank2Points(end:-1:1,:);
end
%% deal with Cross Section Lines data
%****cut the line points outside river polygon****
crossPoints_xyzCell = crossLines;
for n = 1:length(crossPoints_xyzCell)
    xyz_data = crossPoints_xyzCell{n};
    ind_withinRiver = inpolygon(xyz_data(:,1),xyz_data(:,2),...
        bedPoly_xy(:,1), bedPoly_xy(:,2));
    xyz_data = xyz_data(ind_withinRiver,:);
    crossPoints_xyzCell{n} = xyz_data;    
end
clear xyz_data ind_withinRiver n
%*****sort the Cross Section Lines based on their relative location
% L coordinate: length along the river bank (ratio)
% W coordinate: width across the river bank (ratio)
bank1Points_Infill = InfillPoints(bank1Points,4); %increase the density of points
bank2Points_Infill = InfillPoints(bank2Points,4);
crossPointsSortedCell= SortCrossLines(crossPoints_xyzCell,bank1Points_Infill,bank2Points_Infill);
%% split river bed polygon(points inside) based on the cross lines
riverPolySections = ConnectLines2Polygon(bank1Points_Infill,bank2Points_Infill,...
    crossPointsSortedCell);
x_bedPoints = xy_bedPoints(:,1);
y_bedPoints = xy_bedPoints(:,2);
z_bedPoints = nan(size(x_bedPoints));
for i=1:length(riverPolySections)
    crossLineStart= riverPolySections(i).line1;
    crossLineStart = crossLineStart(end:-1:1,:);
    bankLineBottom = riverPolySections(i).line2;
    crossLineEnd= riverPolySections(i).line3;
    bankLineTop = riverPolySections(i).line4;
    bankLineTop = bankLineTop(end:-1:1,:);
    % unify the directions
    channalBound_section = riverPolySections(i).lineAll;
    %%figure(1);figure(2)
%%***build the ralative coordinates for river bank points
    % L coordinate: length along the river bank (ratio)
    % W coordinate: width across the river bank (ratio)
%     bankLineTop = InfillPoints(bankLineTop,4);
%     bankLineBottom = InfillPoints(bankLineBottom,4);
    crossLeft_lw = RelativeCoordsCrossPoints(crossLineStart,0);
    crossRight_lw = RelativeCoordsCrossPoints(crossLineEnd,1);
    crossLeft_lwh = [crossLeft_lw crossLineStart(:,3)];
    crossRight_lwh = [crossRight_lw crossLineEnd(:,3)];
%%figure(3);
%%*** get the matrix index of bed points inside the current polygon
    ind_withinRiver = inpolygon(x_bedPoints,y_bedPoints,...
        channalBound_section(:,1), channalBound_section(:,2));
    bedPoints_section = [x_bedPoints(ind_withinRiver,:),y_bedPoints(ind_withinRiver,:)];
%%figure(4);
%%***calculate the relative coords of river bed points
                % use external function
    [bedPoints_l,bedPoints_w] = RelativeCoordsProjection(bankLineTop,...
    bankLineBottom,bedPoints_section);
%     [bedPoints_l,bedPoints_w] = RelativeCoordsInterp(bankLineTop,bankTop_lw,...
%     bankLineBottom,bankBottom_lw,bedPoints);
%%figure(5)  
    bedPoints_h = interpRatio2Map(crossLeft_lwh,crossRight_lwh,bedPoints_l,bedPoints_w);
    z_bedPoints(ind_withinRiver) = bedPoints_h;
end
%%
bedPoints = z_bedPoints;
bedRaster_z(ind_raster) = bedPoints;
bedRaster = {bedRaster_z,bedRaster_r};
end

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
    [minV_b1c1e,~] = nearestPointDistance(lineB1,xy_c1e);
    [minV_b2c1e,~] = nearestPointDistance(lineB2,xy_c1e);
    if minV_b1c1e < minV_b2c1e %closer to bank1, ban1 is the 2nd line
        secondLine = lineB1;
        fourthLine = lineB2;
    else               %closer to bank2, ban2 is the 2rd line
        secondLine = lineB2;
        fourthLine = lineB1;
    end
    [~,ind_2ndLs] = nearestPointDistance(secondLine,xy_c1e);
    % distance between the 2nd line section and lineC2 end points
    [minV_2ndLC2s,minI_2ndLC2s] = nearestPointDistance(secondLine,xy_c2s);
    [minV_2ndLC2e,minI_2ndLC2e] = nearestPointDistance(secondLine,xy_c2e);
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
    [~,ind_4thLs] = nearestPointDistance(fourthLine,xy_c2e); % nearest point to c2e
    [~,ind_4thLe] = nearestPointDistance(fourthLine,xy_c1s); % nearest point to c1s
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

function inFlag = IsItBetweenTwoPoints(point1,point2,mypoint)% check whether the the intersection point is between the two bank points
if (mypoint(1)-point1(1))*(mypoint(1)-point2(1))<0&&...
        (mypoint(2)-point1(2))*(mypoint(2)-point2(2))<0
    inFlag = 1;
else
    inFlag = 0;
end
end

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

function hq=interpRatio2Map(lwh_line1,lwh_line2, lq,wq)
hq = nan(size(lq));
for i=1:length(lq)
    l_value = lq(i);
    w_value = wq(i);
    % 1D interp h based on w coords
    [~,ia,~]=unique(lwh_line1(:,2));
    h1 = interp1(lwh_line1(ia,2),lwh_line1(ia,3),w_value,'PCHIP');
    [~,ia,~]=unique(lwh_line2(:,2));
    h2 = interp1(lwh_line2(ia,2),lwh_line2(ia,3),w_value,'PCHIP');
    hq(i) = h1+(h2-h1)*l_value;
    
end
end

function cross1Points_lw = RelativeCoordsCrossPoints(cross1Points,l_Value)
%% build the ralative coordinates for river bank points
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
function bank1Points_lw = RelativeCoordsBankPoints(bank1Points,w_Value)
%% build the ralative coordinates for river bank points
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
[~,nearestInd] = nearestPointDistance(xyLine,xy0);
xy1 = xyLine(nearestInd,:);
xyLine(nearestInd,:) = nan;
[~,secondNearestInd] = nearestPointDistance(xyLine,xy_Point);
xy2 = xyLine(secondNearestInd,:);
xyLine(nearestInd,:) = xy1;
vec_a = xy0 - xy1; % vector from nearest Point to xyPoint
vec_b = xy2 - xy1; % vector from nearest Point to second nearest point
vec_b_e = vec_b/norm(vec_b);
vec_a_prj = dot(vec_a,vec_b_e)*vec_b_e; %vector projected on line
xy_prj = vec_a_prj+xy1; %vector to coordinates
% find the location of xy_prj in xyLine
[~,ind1] = nearestPointDistance(xyLine,xy_prj);
xyLine(ind1,:) = nan;
[~,ind2] = nearestPointDistance(xyLine,xy_prj);
indBetween = [ind1,ind2];
end

function crossPointsSortedCell= SortCrossLines...
    (crossPoints_xyzCell,bank1Points,bank2Points)
%sort cross lines based on their relative location along the river bed
crossPointsSortedCell = cell(size(crossPoints_xyzCell));

% check whether the points of the two river banks towards to the same direction 
if abs(bank1Points(1,1)-bank2Points(1,1))>abs(bank1Points(1,1)-bank2Points(end,1))
    bank2Points = bank2Points(end:-1:1,:);
end

% define the bottom bank line
meanL_sort = zeros(length(crossPoints_xyzCell),1);
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
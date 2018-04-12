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


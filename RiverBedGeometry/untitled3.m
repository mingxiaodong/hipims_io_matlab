clc,clear
tic
addpath /Users/b4042552/Dropbox/Matlab/myfunction
load('CarlisleRiverChannelData.mat')
bankLine1 = [channelBoundaryLine(1).X', channelBoundaryLine(1).Y'];
bankLine2 = [channelBoundaryLine(2).X', channelBoundaryLine(2).Y'];
cellsize = 5;
%% BankCrossLinePreprocess
% ?	Points on the two bank lines should be listed towards the same direction (upstream or downstream)
% ?	Points on each cross-section line should be listed towards the same direction (bank1 to bank2)
% ?	Cross-section lines should be listed towards the same direction of bank lines (upstream or downstream)

[bLine1, bLine2,cLines,channelPoly] = BankCrossLinePreprocess(bankLine1, bankLine2, crossSectionLine);
%% River channel segmentation
sectionPoly = Channel2Sections(bLine1,bLine2,cLines);
%% Linearization of channel sections (for each section)
errorDistance = 10;
splitPoints1Cell = cell(length(sectionPoly),1);
splitPoints2Cell = cell(length(sectionPoly),1);
subSectionPolyCell = cell(length(sectionPoly),1);
for i = 1:length(sectionPoly)
    bank1Points = sectionPoly(i).bank1;
    bank2Points = sectionPoly(i).bank2;
    [sP1,sP2,allSubSections] = DiscretizeChannel2Quadrangle(bank1Points,bank2Points,errorDistance);
    splitPoints1Cell{i} = sP1;
    splitPoints2Cell{i} = sP2;
    subSectionPolyCell{i} = allSubSections;
    clear sP1 sP2 sSPoly
end
%% Discretization of channel subsections (for each subsection)
interpGap = (0:0.01:1)';
[channelLinesCell, cutLineCell] = LineInterpolation(cLines,subSectionPolyCell,splitPoints1Cell,splitPoints2Cell,interpGap);
resolution = cellsize;

%% perform interpolation for each riverPolySection generated based on cross section lines

%%

wholeR = makerefmat(min(channelPoly.X),max(channelPoly.Y),cellsize,-cellsize);
wholeCol = round(range(channelPoly.X)/cellsize)+1;
wholeRow = round(range(channelPoly.Y)/cellsize)+1;
wholeRaster_z = nan(wholeRow,wholeCol);
[wholeRaster_x,wholeRaster_y] = Raster2FeaturePoints(wholeRaster_z,wholeR);
rasterCell = cell(length(sectionPoly),1);

for i = 1:length(subSectionPolyCell)
    cutLines = cutLineCell{i};
    allSubSections   = subSectionPolyCell{i};
    sectionPoly1 = sectionPoly(i);
%     bedRaster_r = makerefmat(min(sectionPoly1.X),max(sectionPoly1.Y),cellsize,-cellsize);
%     ncols = round(range(sectionPoly1.X)/cellsize)+1;
%     nrows = round(range(sectionPoly1.Y)/cellsize)+1;
%     secRaster_z = nan(nrows,ncols);
%     [secRaster_x,secRaster_y] = Raster2FeaturePoints(secRaster_z,bedRaster_r);
    for j=1:length(allSubSections)
        crossLine0 = cutLines(j).XYZ;
        crossLine1 = cutLines(j+1).XYZ;
        subSecPoly1 = allSubSections(j);
        in_subsection = inpolygon(wholeRaster_x,wholeRaster_y,subSecPoly1.X,subSecPoly1.Y);
        if sum(in_subsection(:))>0
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
%     ind = ~isnan(secRaster_z);
%     bedPoint_z = secRaster_z(ind);
%     bedPoint_x = secRaster_x(ind);
%     bedPoint_y = secRaster_y(ind);
%     p = map2pix(wholeR,bedPoint_x,bedPoint_y);
%     p = round(p); 
%     p(p<1)=1; p(p(:,1)>wholeRow,1)=wholeRow;
%     p(p(:,2)>wholeCol,2)=wholeCol;
%     [p,ia,~]= unique(p,'rows','stable');
%     bedPoint_z = bedPoint_z(ia);
%     inds = sub2ind(size(wholeRaster_z),p(:,1),p(:,2));
%     wholeRaster_z(inds) = bedPoint_z;
end
toc
%%
% Arcgridwrite('CarlisleRiver.asc',wholeRaster_z,wholeR)
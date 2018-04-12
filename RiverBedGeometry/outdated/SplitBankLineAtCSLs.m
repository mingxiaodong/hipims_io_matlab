function sectionLines = SplitBankLineAtCSLs(bankLine1,bankLine2,CSLs)
% SplitBankLineAtCSLs split bank line into section banks based on the cross
% section lines. bankLine1 and bankLine2 are sorted bank line points. CSLs
% is the sorted cross line cell.
sectionLines = struct('crossLine1',0,'crossLine2',0,'bankline1',0,'bankline2',0);
for i=1:length(CSLs)-1
    if iscell(CSLs)
        csl_1 = CSLs{i};
        csl_2 = CSLs{i+1};
    elseif isstruct(CSLs)
        csl_1 = [CSLs(i).X',CSLs(i).Y',0*CSLs(i).Y'];
        csl_2 = [CSLs(i+1).X',CSLs(i+1).Y',0*CSLs(i+1).Y'];
    end
    if isnan(csl_1(end))
        csl_1(end,:)=[];
    end
    if isnan(csl_2(end))
        csl_2(end,:)=[];
    end
    xyz_Point = csl_1(1,:);
    [xy_prj, indBetween_b1c1] = PrjToBankLine(xyz_Point(1:2),bankLine1);
    z_prj = xyz_Point(3); xyz_prj_b1c1 = [xy_prj z_prj];
    xyz_Point = csl_1(end,:);
    [xy_prj, indBetween_b2c1] = PrjToBankLine(xyz_Point(1:2),bankLine2);
    z_prj = xyz_Point(3); xyz_prj_b2c1 = [xy_prj z_prj];
    xyz_Point = csl_2(1,:);
    [xy_prj, indBetween_b1c2] = PrjToBankLine(xyz_Point(1:2),bankLine1);
    z_prj = xyz_Point(3); xyz_prj_b1c2 = [xy_prj z_prj];
    xyz_Point = csl_2(end,:);
    [xy_prj, indBetween_b2c2] = PrjToBankLine(xyz_Point(1:2),bankLine2);
    z_prj = xyz_Point(3); xyz_prj_b2c2 = [xy_prj z_prj];
    sectionLines(i).crossLine1 = [xyz_prj_b1c1;csl_1;xyz_prj_b2c1];
    sectionLines(i).crossLine2 = [xyz_prj_b1c2;csl_2;xyz_prj_b2c2];
    ind = [indBetween_b1c1(2),indBetween_b1c2(1)];
    sectionLines(i).bankline1 = bankLine1(min(ind):max(ind),:);
    ind = [indBetween_b2c1(2),indBetween_b2c2(1)];
    sectionLines(i).bankline2 = bankLine2(min(ind):max(ind),:);
end
end
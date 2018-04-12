function [z_dem_changed,pix_IndValue] = AmendDEM(Z,R,lineShp,changeValue,varargin)
%Z_New = AmendDEM(Z,R,lineShp,changeValue) change DEM value based on polyline
%   Detailed explanation goes here
% Z is DEM matrix, R is the reference matrix
% lineShp is plolyline struct containing X and Y coords for line points
% changeValue is an array containing the change values for each line
% z_dem_changed = AmendDEM(Z,R,lineShp,changeValue,'Replace') exactly
% replace original DEM with the changeValue.  It is default calculation choice
% z_dem_changed = AmendDEM(Z,R,lineShp,changeValue,'Add') add a positive or
% negtive value to the original DEM value.
% z_dem_changed = AmendDEM(Z,R,lineShp,changeValue,'Positive') only
% increase DEM values, if changeValue is lower than original, ignore the change

% load defenceFailureData
% defenceFailureData_ind = cell(length(defenceFailureData),2);

% Created by Xiaodong Ming at 8/3/2018
z_dem = Z; r_dem = R;
changeLine_shp = lineShp;
if isempty(varargin)
    calMethod = 'Replace';
else
    calMethod = varargin{1};
end
%% find the points to be changed
gridSize = abs(r_dem(2));
z_dem_changed = z_dem;
pix_IndValue = cell(length(changeLine_shp),1);
for i = 1:length(changeLine_shp) % number of lines
    us_cre_lev = min(changeValue(i,:)); % changeValue upstream
    ds_cre_lev = max(changeValue(i,:)); % changeValue downstream
    x_vertex = [changeLine_shp(i).X]'; x_vertex(isnan(x_vertex)) = [];
    y_vertex = [changeLine_shp(i).Y]'; y_vertex(isnan(y_vertex)) = [];
    if size(x_vertex,1) == 1
        x_vertex = x_vertex';
    end
    if size(y_vertex,1) == 1
        y_vertex = y_vertex';
    end
    x_all = x_vertex;
    y_all = y_vertex;
    n = numel(x_all);
    % obtain xy coors from each line segment based on a unit of gridSize
    for j=1:numel(x_vertex)-1 %number of vertex in one polyline
        x1 = x_vertex(j); x2 = x_vertex(j+1);
        y1 = y_vertex(j); y2 = y_vertex(j+1);
        if x1 ~= x2
            xy_slope = (y2-y1)/(x2-x1);
            xy_cosin = 1/sqrt(1+xy_slope^2);
            gridGap = xy_cosin*gridSize*(x2-x1)/abs(x1-x2);
            x_segment = x1:gridGap:x2;
            y_segment = y1 + (x_segment-x1)*xy_slope;
        else %x1==x2 means a vertical line segment
            if y1==y2 %means they are duplicate points, not a line segment
                continue
            else
                gridGap = gridSize*(y2-y1)/abs(y1-y2);
                y_segment = y1:gridGap:y2;
                x_segment = x1 + y_segment*0;
            end
        end
        x_all(n+1:n+numel(x_segment)) = x_segment;
        y_all(n+1:n+numel(x_segment)) = y_segment;
        n = numel(x_all);
    end
    pix_RowCol = map2pix(r_dem,x_all, y_all); % cols and rows
    pix_RowCol = round(pix_RowCol);
    pix_RowCol(pix_RowCol<=0)=nan; % check the rows and cols out of Z range
    pix_RowCol(pix_RowCol(:,1)>size(z_dem,1),1)= nan;
    pix_RowCol(pix_RowCol(:,2)>size(z_dem,2),2)= nan;
    pix_Ind = sub2ind(size(z_dem),pix_RowCol(:,1),pix_RowCol(:,2));
    pix_Ind(isnan(pix_Ind)) = [];
    pix_Ind = unique(pix_Ind,'stable');
    %     pix_Ind = unique(pix_Ind,'stable'); % indice converted from cols and rows and sorted
    %**************change the value for each polyline**************************
    % from upstream to downstream
    if ~isempty(pix_Ind)
        if z_dem_changed(pix_Ind(1))>=z_dem_changed(pix_Ind(end))
            newValues = (linspace(us_cre_lev,ds_cre_lev,numel(pix_Ind)))';
        else
            newValues = (linspace(ds_cre_lev,us_cre_lev,numel(pix_Ind)))';
        end   
        if strcmpi(calMethod,'Positive')
            newValues = max(newValues,z_dem_changed(pix_Ind));
        elseif strcmpi(calMethod,'Add')
            newValues = newValues+z_dem_changed(pix_Ind);
        %else Replace
        end
        z_dem_changed(pix_Ind) = newValues;
        pix_IndValue{i} = [pix_Ind,newValues];
    end
end
end


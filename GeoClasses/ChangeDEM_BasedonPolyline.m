%% find grids passed through by polylines
%% read file
clear,clc
addpath C:/Users/b4042552/Dropbox/Matlab/myfunction
demFileName = 'K:\cudaSWEsSolver\Eden\dsm_10m_small.asc';
shpFileName = 'K:\Data\Cumbria\spatial_flood_defences_Cumm.shp';
changeLine_shp = shaperead(shpFileName);
[z_dem, r_dem] = arcgridread(demFileName);
% load defenceFailureData
% defenceFailureData_ind = cell(length(defenceFailureData),2);
%% plot DEM and Polylines
figure
fig_h1 = mapshow(z_dem, r_dem, 'DisplayType','Surface'); 
zdatam(fig_h1,z_dem-1000)
% demcmap(z_dem)
% axis manual %off
mapshow(changeLine_shp,'Color','red') % show catchment boundary
axis image 
caxis([-20,100])

%% find the points to be changed
gridSize = abs(r_dem(2));
z_dem_changed = z_dem;
for i = 1:length(changeLine_shp) % number of lines
    % only deal with wall and embankment
    if sum(strcmp(changeLine_shp(i).ASSE_TYPE,{'wall','embankment','flood_gate','high_ground'}))
        us_cre_lev = changeLine_shp(i).US_CRE_LEV; % ele upstream
        ds_cre_lev = changeLine_shp(i).DS_CRE_LEV; % ele downstream
    else
        us_cre_lev = changeLine_shp(i).US_CRE_LEV; % ele upstream
        ds_cre_lev = changeLine_shp(i).DS_CRE_LEV; % ele downstream
        %continue
    end
    x_vertex = [changeLine_shp(i).X]'; x_vertex(isnan(x_vertex)) = [];
    y_vertex = [changeLine_shp(i).Y]'; y_vertex(isnan(y_vertex)) = [];
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
    pix_RowCol(pix_RowCol==0)=1; % check the rows and cols out of Z range
    pix_RowCol(pix_RowCol(:,1)>size(z_dem,1),1)= size(z_dem,1);
    pix_RowCol(pix_RowCol(:,2)>size(z_dem,2),2)= size(z_dem,2);
    pix_Ind = sub2ind(size(z_dem),pix_RowCol(:,1),pix_RowCol(:,2));
    pix_Ind = unique(pix_Ind,'stable'); % indice converted from cols and rows and sorted
%**************change the value for each polyline
    % from upstream
    if z_dem_changed(pix_Ind(1))>=z_dem_changed(pix_Ind(end)) 
        changeValue = linspace(us_cre_lev,ds_cre_lev,numel(pix_Ind));
    else
        changeValue = linspace(ds_cre_lev,us_cre_lev,numel(pix_Ind));
    end
    if min(changeValue)>0
        z_dem_changed(pix_Ind) = changeValue+0.2;
    end
%     for n=1:length(defenceFailureData)
%         if changeLine_shp(i).OBJECTID == defenceFailureData(n).ObjectID
%             %z_dem_changed(pix_Ind) = 0;
%             defenceFailureData_ind{n,1} = pix_Ind;
%             defenceFailureData_ind{n,2} = [x_all y_all];
%         end
%     end
end
z_dem_changed = max(z_dem_changed,z_dem);
% sizeZmat = size(z_dem);
% save defenceFailureData_ind defenceFailureData_ind sizeZmat r_dem
%% show the changed DEM map with polylines
figure
fig_h2 = mapshow(z_dem_changed1, r_dem, 'DisplayType','mesh'); 
% zdatam(fig_h2,z_dem_changed-1000)
% demcmap(z_dem_changed)
% axis manual %off
hold on
mapshow(changeLine_shp,'Color','red') % show shp polyline
hold off
% axis image
% caxis([-20,100])
%% write changed DEM file
newDEM_name = 'dsm_10m_addDyke.asc';%dsm_DF4 defence failure
Arcgridwrite(newDEM_name,z_dem_changed,r_dem)
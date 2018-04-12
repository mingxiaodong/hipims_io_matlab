% to combine DEM asc files downloaded from digimap
clear,clc
gridShp = shaperead('C:\Users\b4042552\Google Drive\MyResearch\London\T10kmGrid_LC4.shp');
% filedir = 'F:\Data\Environment Agency (National 2m)\LC4\';
filedir = 'F:\Data\Environment Agency (National 2m)\LC4_Large\';
outputdir = filedir;
cd(filedir) % position of the asc files
%% get the file list
% cellsize of each asc file should be the same
% for n=1:length(gridShp)
    tic
%     writefileName = [outputdir gridShp(n).TILE_NAME '_2m.asc'];
    writefileName = [outputdir 'DSM10m_LC4.asc'];
%     subfileNames = [gridShp(n).TILE_NAME(1:3) '*' gridShp(n).TILE_NAME(4) '*_DSM_2m.asc'];
    subfileNames = '*_10m.asc';
    mylist = dir(subfileNames);% list of the files in current folder
%     for i=1:length(mylist)
%         if strcmpi(mylist(i).name([1:3,5]),gridShp(n).TILE_NAME)
%             mylist(i).isdir = true;
%         end
%     end
%     mylist = mylist([mylist.isdir]);
    %%*read data from DEM files
    delimiter = ' ';
    endRow = 6; % number of rows in the head of asc file
    formatSpec = '%s%f';
    RefValues = zeros(endRow,length(mylist)); % spatial reference value for each asc file
    fileExtent = nan(length(mylist),4); % for each asc file: extent left, right, bottom, top
    Z_Data_cell = cell(1,length(mylist));
    for i = 1:length(mylist)
        filename = mylist(i).name;
        % read the reference rows
        fileID = fopen(filename,'r');
        fileRef = textscan(fileID, formatSpec, endRow, 'Delimiter', delimiter,'MultipleDelimsAsOne', true);
        fclose(fileID);
        RefValues(:,i) = fileRef{2};
        ncols = RefValues(1,i);
        nrows = RefValues(2,i);
        xllcorner = RefValues(3,i);
        yllcorner = RefValues(4,i);
        cellsize = RefValues(5,i);
        extent_left = xllcorner+cellsize/2;
        extent_right = extent_left+cellsize*(ncols-1);
        extent_bottom = yllcorner+cellsize/2;
        extent_top = extent_bottom+cellsize*(nrows-1);
        fileExtent(i,:) = [extent_left extent_right extent_bottom extent_top];
        Z_Data = dlmread(filename,' ',endRow,0);
%         Z_Data = Z_Data(:,2:end); %********************as the first col is space!!
        Z_Data(Z_Data==RefValues(6,i))=nan; % nodata value
        if size(Z_Data,1)==nrows && size(Z_Data,2)==ncols
            Z_Data_cell{i} = Z_Data;
        else
            error('The size of Z data can not match ncols and nrows, please check the DEM file!')
        end     
    end
    % check cellsize
    if RefValues(5,:)==cellsize
        cellsize_new = cellsize;
    else
        error('cellsize of each asc file is not the same')
    end
    % define parameters of new DEM file
    extent_left_new = min(fileExtent(:,1));
    extent_right_new = max(fileExtent(:,2));
    extent_bottom_new = min(fileExtent(:,3));
    extent_top_new = max(fileExtent(:,4));
%     xllcorner_new = extent_left_new-cellsize_new/2; %lower-left corner
%     yllcorner_new = extent_bottom_new-cellsize_new/2;
    ncols_new = numel(extent_left_new:cellsize_new:extent_right_new);
    nrows_new = numel(extent_bottom_new:cellsize_new:extent_top_new);
    x11_new = extent_left_new; % centre of the fisrt pixel
    y11_new = extent_top_new;
    refMat_new = makerefmat(x11_new,y11_new,cellsize_new,-cellsize_new);
    z_new = nan(nrows_new,ncols_new);
    % pass sub Z data to combine new Z data
    for i = 1:length(mylist)
        extent_left = fileExtent(i,1);
        extent_top = fileExtent(i,4);
        [row_top,col_left] = map2pix(refMat_new,extent_left,extent_top);
        row_top = round(row_top); col_left = round(col_left);
        sub_zData = Z_Data_cell{i};
        [sub_nrows,sub_ncols] = size(sub_zData);
        z_new(row_top:row_top+sub_nrows-1, col_left:col_left+sub_ncols-1)...
            = sub_zData;
    end
    Arcgridwrite(writefileName,z_new,refMat_new)
%     disp([num2str(n) ' file created: ' gridShp(n).TILE_NAME '_2m.asc'])
    toc
% end
%% 
[Z,R] = arcgridread(writefileName);
figure; mapshow(Z,R,'DisplayType','Surface'); axis image
%%
clear,clc
cd('F:\Data\Environment Agency (National 2m)\LC4_Large')
mylist = dir('*_2m.asc');
newResolution = 10;
for n=1:length(mylist)
    [Z,R] = arcgridread(mylist(n).name);
    [z_New,r_new] = ResampleRaster(Z,R,newResolution);
    newFileName = [mylist(n).name(1:end-6) '10m.asc'];
    Arcgridwrite(newFileName,z_New,r_new)
    disp(['file ' num2str(n) ' resampled: ' newFileName])
end



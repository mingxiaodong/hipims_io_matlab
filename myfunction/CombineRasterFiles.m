function [Z_new,R_New]=CombineRasterFiles(inputFolder,outputFile)
% to combine DEM asc files downloaded from digimap
currentDir = pwd;
cd(inputFolder) % position of the asc files
%% get the file list
% cellsize of each asc file should be the same
subfileNames = '*.asc';
mylist = dir(subfileNames);% list of the files in current folder
%%*read data from DEM files
fileExtent = zeros(length(mylist),4); % spatial reference value for each asc file
for i = 1:length(mylist)
    filename = mylist(i).name;
    % read the reference rows
    fileID = fopen(filename,'r');
    headValues = textscan(fileID, '%s%f', 6, 'Delimiter', ' ','MultipleDelimsAsOne', true);
    fclose(fileID);
    headValues = headValues{2};
    ncols = headValues(1);
    nrows = headValues(2);
    xllcorner = headValues(3);
    yllcorner = headValues(4);
    cellsize = headValues(5);
    %     nonedataValue = headValues(6);
    left = xllcorner; % left edge
    right = xllcorner+cellsize*ncols; % right edge
    top = yllcorner+cellsize*nrows; % top edge
    bottom = yllcorner; % bottom edge
    fileExtent(i,:) = [left,right,bottom,top];
end

% define parameters of new DEM file
extent_left = min(fileExtent(:,1));
extent_right = max(fileExtent(:,2));
extent_bottom = min(fileExtent(:,3));
extent_top = max(fileExtent(:,4));
x11 = extent_left+cellsize/2;
y11 = extent_top-cellsize/2;
ncols = (extent_right-extent_left)/cellsize;
nrows = (extent_top-extent_bottom)/cellsize;
R_New = makerefmat(x11,y11,cellsize,-cellsize);
Z_new = nan(nrows,ncols);
% pass sub Z data to combine new Z data
for i = 1:length(mylist)
    [Z,R] = ArcgridreadM(mylist(i).name);
    x11_s = R(3,1)+R(2,1);
    y11_s = R(3,2)-R(2,1);
    [row1,col1] = map2pix(R_New,x11_s,y11_s);
    row1 = round(row1); col1 = round(col1);
    [sub_nrows,sub_ncols] = size(Z);
    disp([row1,row1+sub_nrows-1, col1,col1+sub_ncols-1])
    if row1<0
        disp(mylist(i).name)
        pause
    else
        Z_new(row1:row1+sub_nrows-1, col1:col1+sub_ncols-1)= Z;
    end
    %disp(i)
end
Arcgridwrite(outputFile,Z_new,R_New)
cd(currentDir)

end


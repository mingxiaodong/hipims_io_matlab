%% transfer result grid data to ascii file
clear,clc
curFolder = '/Users/b4042552/Google Drive/Data';
cd(curFolder)
% read dem.txt to obtain the asc title strings and values
filename = 'DEM.txt';
delimiter = ' ';
endRow = 6;
formatSpec = '%s%f';
fileID = fopen(filename,'r');
fileRef = textscan(fileID, formatSpec, endRow,'Delimiter', delimiter,'MultipleDelimsAsOne', true);
fclose(fileID);
RefValues = fileRef{2};
[Z,R] = arcgridread(filename);

%%
load('ResultCell_H72.mat', 'ResultCell_H','TimeSeries')
for n = 1:3%length(TimeSeries)
    T = TimeSeries(n);
    writeName = ['h_grid' num2str(T) 's.asc'];
    fileID = fopen(writeName,'w');
    for i=1:6
        forspec = '%-12s  %d\n';
        if i==3||i==4; forspec = '%-12s  %.4f\n'; end
        fprintf(fileID,'%-12s  %d\n', fileRef{1}{i}, RefValues(i));
    end
    Z_N = ResultCell_H{n}; Z_N(isnan(Z_N))= RefValues(end);
    dlmwrite(writeName,Z_N,'-append','delimiter',delimiter)
    fclose(fileID);
end

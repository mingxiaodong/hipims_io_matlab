%% Initialize variables.
clear,clc
T = 28800;
Filename = {'hU','eta'};

for i=1:2;
%%load backupfile
filpath1 = 'G:\cudaSWEsSolver\ThamesValley\output';
cd(filpath1)
BackupResult = dlmread([Filename{i} '_backup__' num2str(T) '.dat']); BackupResult(:,1:2)=[];

%load old initial files
filpath2 = 'G:\cudaSWEsSolver\ThamesValley1\input\field';
cd(filpath2)
startRow = 4;
endRow = startRow+length(BackupResult)-1;
if length(Filename{i})==2
    formatSpec = '%f%f%f%[^\n\r]'; formatSpec1 = '%-12d %.8f %.8f\n';
else
    formatSpec = '%f%f%[^\n\r]'; formatSpec1 = '%-12d %.8f\n'; 
end

fileID = fopen([Filename{i} '.dat'],'r');
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines', startRow-1, 'ReturnOnError', false);
fclose(fileID);
InitialValue = [dataArray{1:end-1}];
id_IV = [InitialValue(:,1) BackupResult];

fileID = fopen([Filename{i} '.dat'],'r');
filetail = textscan(fileID, '%s', 'Delimiter', '', 'MultipleDelimsAsOne', true, 'HeaderLines', length(BackupResult)+3, 'ReturnOnError', false);
filetail = filetail{1};
fclose(fileID);
%%write data
fileID = fopen([Filename{i} '1.dat'],'w');
% print valid cell and their initial value
fprintf(fileID,'$Element Number\n'); fprintf(fileID,'%d\n',length(id_IV));
fprintf(fileID,'$Element_id  Value\n'); fprintf(fileID,formatSpec1,id_IV');
for n=1:length(filetail)-1
    fprintf(fileID,'%s\n',filetail{n});
end
fprintf(fileID,'%s',filetail{end});
fclose(fileID);
end
%% create rainfall source file - single file with a matrix
clear,clc
cd('G:\Data\UK_EU\UK_Rainfall_Radar\201512')
mylist = dir('2015*.asc');
%%
[Z,~]=arcgridread(mylist(1).name);
m = length(mylist);
n = numel(Z);
T = 0:5*60:(m-1)*5*60;
RainData = zeros(m,n);
for i=1:m
    [Z,~]=arcgridread(mylist(i).name);
    Z(isnan(Z)) = 0; Z(Z<0) = 0;
    RainData(i,:) = Z(:);
end
T_RainMat = [T' RainData/3600/1000];
%%
T0 = datetime(mylist(1).name(1:12),'InputFormat','yyyyMMddHHmm');
T_datetime = T0 +T'/3600/24;
plot(T_datetime,T_RainMat(:,size(T_RainMat,2)/2))
%%

writeName = ['precipitation_source_all_' mylist(1).name(1:end-4) '.dat'];
dlmwrite(writeName,n);
dlmwrite(writeName,T_RainMat,'-append','delimiter',' ');
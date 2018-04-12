% download flood monitoring data from environment.data.gov.uk via URL
clc,clear
load CarlisleFlowStation % load Carlisle flow station information
url_head = 'http://environment.data.gov.uk/flood-monitoring/id/stations/';
url_middle = '/readings.csv?_sorted';
url_tail = '&_limit=10000';
dataCell = cell(height(CarlisleFlowStation),1);
for i = 1:height(CarlisleFlowStation)
    apiID = char(CarlisleFlowStation.APIID(i)); %'765013';
    startDate = ['&startdate=' datestr(datetime(2017,10,10),'yyyy-mm-dd')];
    endDate   = ['&enddate=',  datestr(datetime(2017,10,16),'yyyy-mm-dd')];
    urlstr = [url_head apiID url_middle startDate endDate url_tail];
    myreadtable = @(filename)readtable(filename,'HeaderLines',0, ...
        'Format','%s%s%f','Delimiter',',','MultipleDelimsAsOne',1);
    options = weboptions('ContentReader',myreadtable);
    data = webread(urlstr,options);
    dataCell{i} = data;
end
%%plot observations for each station
figure;
for i = 1:length(dataCell)
    data = dataCell{i};
    station_DateTime = datetime(data.dateTime,...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''');
    station_LevelStage = data.value;
    plot(station_DateTime,station_LevelStage)
    xlabel('data&time')
    ylabel('level stage')
    title(CarlisleFlowStation.Station(i))
    pause(0.5)
end
%%
urlstr = 'http://data.ceda.ac.uk/badc/ukmo-nimrod/data/composite/uk-5km/2017/';
options = weboptions('Timeout',60);
data = webread(urlstr,options);
splitStr = regexp(data,'\n','split');
expression = '<a href="\w{1,}\.tar';
startIndex = regexp(data,expression,'match');
%%
urlstr_newtail = startIndex{1}; urlstr_newtail = urlstr_newtail(10:end-1);
urlstr_new = [urlstr urlstr_newtail];
% expression = '\n';
% splitStr = regexp(data,expression,'split');
% filename = 'samples.txt';
%%
% myreadtable = @(filename)textscan(filename,'HeaderLines',0, ...
%         'Format','%s%s%f','Delimiter',',','MultipleDelimsAsOne',1);
options = weboptions('ContentType','text');
filename = 'samples.txt';
dataNew = webread(urlstr_new);
outfilename = websave(filename,urlstr_new,options);
% regexp
%%
URL = 'http://www.mathworks.com/matlabcentral/fileexchange';
filename = 'samples.html';
urlwrite(URL,filename,'get',{'term','urlwrite'});
%%
clear,clc
addpath C:\Users\b4042552\Dropbox\TempData\YRD
load CountyEvents153010_added %Events 140*1 cell added with date characters
% the code for adding date string to CountyEvent
% for i = 1:140
%     for j = 1:length(CountyEvent{i})
%         CountyEvent{i}(j).DateStartC = datestr(CountyEvent{i}(j).DateStart);
%         CountyEvent{i}(j).DateEndC = datestr(CountyEvent{i}(j).DateEnd);
%     end
% end
load CountyInfor

% load raw data of Loss records
[~, ~, raw, dates] = xlsread('Loss_Selected.xlsx','Sheet4','N2:Q1406','',@convertSpreadsheetExcelDates);
Code = [raw{:,1}]'; %county code
DateStart = datetime([dates{:,2}].', 'ConvertFrom', 'Excel', 'Format', 'dd/MM/yyyy');
DateEnd = datetime([dates{:,3}].', 'ConvertFrom', 'Excel', 'Format', 'dd/MM/yyyy');
EcoLoss = [raw{:,4}]';
LossRecord = table(Code,DateStart,DateEnd,EcoLoss);
clearvars -except CountyEvent CountyInfor LossRecord

%%
LossHazard = LossRecord;
Rain = zeros(height(LossHazard),1);
Wind = zeros(height(LossHazard),1);
for i = 1:height(LossHazard)    
    fid = find([CountyInfor.code]==LossRecord.Code(i)); % code of the county where the event happened
    if isempty(fid)
        continue
    end    
    CountyStruct = CountyEvent{fid};
    EventDate = datenum(LossRecord.DateStart(i));
    % to identify whether the Date of record loss is happened during an events
    DateIdnt1 = [CountyStruct.DateStart]'-EventDate;
    DateIdnt2 = [CountyStruct.DateEnd]'-EventDate;    
    [M,I] = min(DateIdnt1.*DateIdnt2);   
    if M>0
        [M1,I1] = min(abs(DateIdnt1));
        [M2,I2] = min(abs(DateIdnt2));
        erroN = 6;
        if M1<=M2 && M1<=erroN
            I = I1;    
        elseif M1 > M2 && M2<=erroN
            I = I2;          
        else
            I = NaN;
        end
    end
    if isnan(I)
       Rain(i) = NaN;
       Wind(i) = NaN;
    else
        Rain(i) = CountyStruct(I).MaxRain;
        Wind(i) = CountyStruct(I).MaxWind;%matched record
    end
end
LossHazard = [LossHazard,table(Rain,Wind)];
LossHazard_Fit = LossHazard;
LossHazard_Fit(isnan(LossHazard_Fit.Rain),:)=[];
LossHazard_Fit(LossHazard_Fit.Rain==0,:)=[];
%%
FitR = 538:564;
T = LossHazard_Fit(FitR,:);
T = sortrows(T,'DateStart','ascend');
X = T.Rain;
Y = T.Wind;
Z = T.EcoLoss;
cftool
%%
% load raw data of Loss records
[~, ~, raw, dates] = xlsread('EcoLoss_Case.xlsx','Sheet1','A2:I89','',@convertSpreadsheetExcelDates);
Code = [raw{:,1}]'; %county code
DateStart = datetime([dates{:,2}].', 'ConvertFrom', 'Excel', 'Format', 'dd/MM/yyyy');
DateEnd = datetime([dates{:,3}].', 'ConvertFrom', 'Excel', 'Format', 'dd/MM/yyyy');
LossRatio = [raw{:,9}]';
EcoLoss = [raw{:,4}]';
Rain = [raw{:,5}]';
Wind = [raw{:,6}]';
EcoLossCase = table(Code, DateStart, DateEnd, EcoLoss, Rain, Wind, LossRatio);
clearvars raw dates DateStart DateEnd LossRatio EcoLoss Rain Wind
%%
T2 = EcoLossCase(63:end,:);
X2 = T2.Rain;
Y2 = T2.Wind;
Z2 = T2.LossRatio;
%%
CountyCounts = struct2table(CountyInfor);
Counts = zeros(140,1);
for i = 1:140
    Counts(i) = sum(LossHazard_Fit.Code==CountyInfor(i).code);
end
CountyCounts = [CountyCounts, table(Counts)];
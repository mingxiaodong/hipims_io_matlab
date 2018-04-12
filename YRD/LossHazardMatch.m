%% 19-8-2015 match loss record with an multi-hazard event for each county
clear,clc
% load loss record data based on events
% [county code, year, month, day, loss ratio]
LossRecord = xlsread('LossRecord_Multi_disaster.xlsx','A2:E139');
% load multi-hazard events data for each county
load CountyEvents153010 %Events 140*1 cell
% add date string to CountyEvent
for i = 1:140
    for j = 1:length(CountyEvent{i})
        CountyEvent{i}(j).DateStartC = datestr(CountyEvent{i}(j).DateStart);
        CountyEvent{i}(j).DateEndC = datestr(CountyEvent{i}(j).DateEnd);
    end
end
load CountyInfor
%%
LossHazard = LossRecord;
for i = 1:length(LossRecord)
    
    fid1 = find([CountyInfor.code]==LossRecord(i,1)); % code of the county where the event happened
    CountyStruct = CountyEvent{fid1};     
    EventDate = datenum(LossRecord(i,2),LossRecord(i,3),LossRecord(i,4));
    
    % to identify whether the Date of record loss is happened during an events
    DateIdnt1 = [CountyStruct.DateStart]'-EventDate;
    DateIdnt2 = [CountyStruct.DateEnd]'-EventDate;    
    [M,I] = min(DateIdnt1.*DateIdnt2);
    
    if M<=0
        LossHazard(i,6) = CountyStruct(I).MaxRain;
        LossHazard(i,7) = CountyStruct(I).MaxWind;%matched record
    else
        [M1,I1] = min(abs(DateIdnt1));
        [M2,I2] = min(abs(DateIdnt2));
        erroN = 20;
        if M1<=M2 && M1<=erroN
            LossHazard(i,6) = CountyStruct(I1).MaxRain;
            LossHazard(i,7) = CountyStruct(I1).MaxWind;        
        elseif M1 > M2 && M2<=erroN
            LossHazard(i,6) = CountyStruct(I2).MaxRain;
            LossHazard(i,7) = CountyStruct(I2).MaxWind;            
        else
            LossHazard(i,6) = NaN;
            LossHazard(i,7) = NaN;
        end
    end
end
LossHazard_Fit=LossHazard;
LossHazard_Fit(isnan(LossHazard(:,6)),:)=[];
%%
%LossHazard_Fit(LossHazard_Fit(:,5)<0.02,:)=[];
FitR = 1:length(LossHazard_Fit);
X = LossHazard_Fit(FitR,6);
Y = LossHazard_Fit(FitR,7);
Z = LossHazard_Fit(FitR,5);
cftool
%% add date string to CountyEvent
for i = 1:140
    for j = 1:length(CountyEvent{i})
        CountyEvent{i}(j).DateStartC = datestr(CountyEvent{i}(j).DateStart);
        CountyEvent{i}(j).DateEndC = datestr(CountyEvent{i}(j).DateEnd);
    end
end
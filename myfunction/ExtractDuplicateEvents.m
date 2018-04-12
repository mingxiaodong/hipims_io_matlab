function PureEvent = ExtractDuplicateEvents(Dates,Observations,Daygap)
% PureEvent = ExtractDuplicateEvents remove the multiple events in the 
%   short duration defined by Daygap. The maximum value for each variable 
%   will be selected as the represent values in that duration.
% Dates: datetime array gives the time of observations
% Observations: martix with columns representing the variables
% Daygap: integer give the duration of extraction
% PureEvent: return value, a table of variable Date and Event
[Dates,I] = sort(Dates);
Observations = Observations(I,:);
gaps = [Daygap+1;Dates(2:end)-Dates(1:end-1)];
indGaps = gaps<=Daygap;
Date = Dates;
Event = Observations;
for i=1:length(Dates)
    if indGaps(i)
       Date(i) = Date(i-1);
       Event(i,:) = max(Event(i-1:i,:));
%        Event(i-1,:) = Event(i,:);
    end
end
PureEvent = table(Date,Event);
[~,ia,~] = unique(Date,'last');
PureEvent = PureEvent(ia,:);
end
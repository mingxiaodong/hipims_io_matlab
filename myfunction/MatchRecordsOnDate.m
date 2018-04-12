function matchedRec = MatchRecordsOnDate(rec1,rec2,varargin)
% matchedRec = MatchRecordsOnDate(rec1,rec2,daylag)
% first col of rec1 and rec2 are datetime vector, second cols are value to
% be matched
% daylag: 1 use the one day later obs from rec2 to match rec1
%        -1 use the one day before obs from rec2 to match rec1
% Created by X Ming on 2018-2-14
if isempty(varargin)
    daylag = 0;
    rec3 = [];
elseif length(varargin)==1
    if numel(varargin{1})==1
        daylag = varargin{1};
        rec3 = [];
    else
        daylag = 0;
        rec3 = varargin{1};
    end
elseif length(varargin)==2
    if numel(varargin{1})==1
        daylag = varargin{1};
        rec3 = varargin{2};
    else
        daylag = varargin{2};
        rec3 = varargin{1};
    end
end
if istable(rec1)
    rec1D = table2array(rec1(:,1));
    rec2D = table2array(rec2(:,1));
    rec1V = table2array(rec1(:,2));
    rec2V = table2array(rec2(:,2));

else
    rec1D = rec1(:,1);
    rec2D = rec2(:,1);
    rec1V = rec1(:,2);
    rec2V = rec2(:,2);
end
if ~isempty(rec3)
    rec3D = table2array(rec3(:,1));
    rec3V = table2array(rec3(:,2));
    if numel(unique(rec3D))<numel(rec3D)
        error('Arg3 has duplicate records')
    end
end
if numel(unique(rec1D))<numel(rec1D)
    error('Arg1 has duplicate records')
end
if numel(unique(rec2D))<numel(rec2D)
    error('Arg2 has duplicate records')
end
rec2D = rec2D+daylag;
if isempty(rec3)
    dateAll = [rec1D;rec2D];
else
    dateAll = [rec1D;rec2D;rec3D];
end
dateS = min(dateAll);
dateE = max(dateAll);
Date = (dateS:1:dateE)';
v1 = nan(size(Date));
v2 = v1;
v3 = v1;
ind = datenum(rec1D-dateS)+1;
v1(ind) = rec1V;
ind = datenum(rec2D-dateS)+1;
v2(ind) = rec2V;
if ~isempty(rec3)
    ind = datenum(rec3D-dateS)+1;
    v3(ind) = rec3V;
end
if isempty(rec3)
    matchedRec = table(Date,v1,v2);
    ind = isnan(v1)|isnan(v2);
else
    matchedRec = table(Date,v1,v2,v3);
    ind = isnan(v1)|isnan(v2)|isnan(v3);
end
matchedRec(ind,:) = [];
end
    
function [chi_u_sigV,chi_uB_sigV,u] = ChiDependenceSigTest(matchedRec,sampleTimes)
%% *significance test (in consideration to the seasonality)
% permutation test
% resample without replacement
% generate year block data;
% [chi_u,chiB_u]= QuantileAsymptoticalDependence(u,X,Y,sampleTimes);
recDate = [];
if istable(matchedRec)
    if width(matchedRec)==3
        recDate = table2array(matchedRec(:,1));% datetime
        X = table2array(matchedRec(:,2));
        Y = table2array(matchedRec(:,3));
    else
        X = table2array(matchedRec(:,1));
        Y = table2array(matchedRec(:,2));
    end
else
    X = matchedRec(:,1);
    Y = matchedRec(:,2);
end
[~,~,u] = ChiDependence([X,Y],100);
chi_u_re_all = zeros(length(u),sampleTimes);
chi_uB_re_all = chi_u_re_all;
if ~isempty(recDate) % consider seasonality
    yearAll = year(recDate);
    yearAllUniq = unique(yearAll);
    yearAllUniq(isnan(yearAllUniq))=[];
    dataYearBlock = cell(2,numel(yearAllUniq));
    for i=1:length(yearAllUniq)
        dateind = yearAll==yearAllUniq(i);
        dateOneYear = recDate(dateind);
        ind = datenum(dateOneYear-datetime(yearAllUniq(i),1,1))+1;
        data = [ind X(dateind) Y(dateind)];
        data(ind>365,:) = [];
        data1 = nan(365,1); data2 = data1;
        data1(data(:,1)) = data(:,2);
        data2(data(:,1)) = data(:,3);
        dataYearBlock{1,i} = data1;
        dataYearBlock{2,i} = data2;
        clear data1 data2 data
    end

    n = length(yearAllUniq);
    dataX = cell2mat(dataYearBlock(1,:)');
    for i=1:sampleTimes %repeat sampling sampleTimes
        % randomly select from X and Y
        yearInd = datasample(1:n,n,'Replace',false);
        datacell2 = dataYearBlock(2,yearInd);
        dataY = cell2mat(datacell2');
        [chi_u_re,chi_uB_re,u0] = ChiDependence([dataX,dataY],50);
        %     [chi_u_re,chi_uB_re]= QuantileAsymptoticalDependence(u,dataX,dataY);
        chi_u_re_all(:,i) = interp1(u0,chi_u_re(:,2),u,'linear');
        chi_uB_re_all(:,i) = interp1(u0,chi_uB_re(:,2),u,'linear');
    end
else % don't consider seasonality
    for i=1:sampleTimes %repeat sampling sampleTimes
        % randomly select from X and Y
        dataX = X;
        n = length(dataX);
        dayInd = datasample(1:n,n,'Replace',false);
        dataY = Y(dayInd);
        [chi_u_re,chi_uB_re,u0] = ChiDependence([dataX,dataY],50);
        %     [chi_u_re,chi_uB_re]= QuantileAsymptoticalDependence(u,dataX,dataY);
        chi_u_re_all(:,i) = interp1(u0,chi_u_re(:,2),u,'linear');
        chi_uB_re_all(:,i) = interp1(u0,chi_uB_re(:,2),u,'linear');
    end
end



chi_u_re_all = sort(chi_u_re_all,2,'descend');
chi_uB_re_all = sort(chi_uB_re_all,2,'descend');
% the 10th largest value is the 5% significance level
chi_u_sigV = chi_u_re_all(:,round(sampleTimes*0.05));
% chi_u_sigV(chi_u_sigV>=1|chi_u_sigV<=0)=nan;
chi_uB_sigV = chi_uB_re_all(:,round(sampleTimes*0.05));
% chi_uB_sigV(chi_uB_sigV>=1|chi_uB_sigV<=0)=nan;

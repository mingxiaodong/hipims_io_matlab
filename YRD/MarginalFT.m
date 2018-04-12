%% 19-8-2015 boundary fitting of rain and wind based on county events
clear,clc
load CountyEvents153010 %multi-hazard events in each county
FT_Result = struct('Rain',[0;0],'Wind',[0;0],'RainGEV',[0 0 0],'WindGEV',[0 0 0],...
    'RainKSTest',[0 0],'WindKSTest',[0 0]); %fitting and testing result of rain and wind

for i = 1:140 % county code    
    X = [CountyEvent{i,1}.MaxRain]';    
    Y = [CountyEvent{i,1}.MaxWind]';
    rcoef=gevfit(X); 
    wcoef=gevfit(Y);
    FT_Result(i).Rain = X;
    FT_Result(i).Wind = Y;
    FT_Result(i).RainGEV = rcoef; % rain distribution coef
    FT_Result(i).WindGEV = wcoef; % wind distribution coef
    % KS test
    test_cdf1 = [X,cdf('Generalized Extreme Value',X,rcoef(1),rcoef(2),rcoef(3))];
    test_cdf2 = [Y,cdf('Generalized Extreme Value',Y,wcoef(1),wcoef(2),wcoef(3))];
    [hx,px] = kstest(X,'CDF',test_cdf1); 
    FT_Result(i).RainKSTest = [hx,px];
    [hy,py] = kstest(Y,'CDF',test_cdf2);
    FT_Result(i).WindKSTest = [hy,py];
    % Chi2 test
    [hx,px] = chi2gof(X,'cdf',{@gevcdf,rcoef(1),rcoef(2),rcoef(3)},'nbins',10);
    FT_Result(i).RainChi2Test = [hx,px];
    [hy,py] = chi2gof(Y,'cdf',{@gevcdf,wcoef(1),wcoef(2),wcoef(3)},'nbins',10);
    FT_Result(i).WindChi2Test = [hy,py];
  
end
save FT_Result153010 FT_Result RainMin1 RainMin2 WindMin
%% to see the number of counties do not pass the test
RainKSTest = cell2mat({FT_Result.RainKSTest});
WindKSTest = cell2mat({FT_Result.WindKSTest});
disp('KS Test result for Rain')
disp(sum(RainKSTest(:,1)))
disp('KS Test result for Wind')
disp(sum(WindKSTest(:,1)))
RainChi2Test = cell2mat({FT_Result.RainChi2Test});
WindChi2Test = cell2mat({FT_Result.WindChi2Test});
disp('Chi2 Test result for Rain')
disp(sum(RainChi2Test(:,1)))
disp('Chi2 Test result for Wind')
disp(sum(WindChi2Test(:,1)))
%% Distribution fitting and testing for Rain
fid2 = find(WindKSTest==1);
Test2 = ones(size(fid2));
for n = 1:length(fid2)
    i = fid2(n);
    Y = [Events{i,1}.MaxWind]';
    wcoef=gevfit(Y);
    Test2(n) = chi2gof(Y,'cdf',{@gevcdf,wcoef(1),wcoef(2),wcoef(3)},'nbins',10);
end
disp('Chi2 Test results for Wind')
disp(sum(Test2))

%%
i = 112;
X = [CountyEvent{i,1}.MaxRain]';
rcoef = FT_Result(i).RainGEV;
[hx,px] = chi2gof(X,'cdf',{@gevcdf,rcoef(1),rcoef(2),rcoef(3)},'nbins',10); 
disp('The results of Chi-square Test for X:')
disp(['hx = ',num2str(hx)])


%% Plot the fitting result of Rain
i = 123;
X = [CountyEvent{i,1}.MaxRain]';
rcoef = FT_Result(i).RainGEV;
figure
subplot(2,1,1)  
[fx,xc] = ecdf(X);
ecdfhist(fx,xc)    % Histogram based on empirical cumulative distribution
xlabel('X - Rainfall (mm)');   
ylabel('f(x)');         
title('Histogram and Density Curve of Frequency Distribution for Rainfall')
hold on
xr=min(X):max(X);         % value of X (rainfall)
yrp=gevpdf(xr,rcoef(1),rcoef(2),rcoef(3));         % GEV probability density distribution
plot(xr,yrp,'r','LineWidth',2)  % plot probability density curve of rainfall
hold off
subplot(2,1,2)
cdfplot(X)
xlabel('X - Rainfall (mm)');   
ylabel('F(x)');
hold on
yrc=gevcdf(xr,rcoef(1),rcoef(2),rcoef(3));   % GEV cumulative distribution
plot(xr,yrc,'r','LineWidth',2)
legend('Empirical','Theoretical','Location','SE')
title('Empirical and Theoretical Distribution of Rainfall')
hold off

%% Plot the fitting result of Wind
i = 123;
Y = [CountyEvent{i,1}.MaxWind]';
wcoef = FT_Result(i).WindGEV; 
[hy,py]=chi2gof(Y,'cdf',{@gevcdf,wcoef(1),wcoef(2),wcoef(3)},'nbins',10,'Alpha',0.90);
disp('The results of Chi-square Test for Wind:')
disp(['hy = ',num2str(hy)])
figure
subplot(2,1,1) 
[fy, yc] = ecdf(Y);    %
ecdfhist(fy, yc, 10); %
xlabel('Y - Wind Speed (m/s)');  %
ylabel('f(y)');  %
title('Histogram and Density Curve of Frequency Distribution for Wind Speed')
hold on
xw=min(Y):0.1:max(Y);    %
ywp=gevpdf(xw,wcoef(1),wcoef(2),wcoef(3)); %
plot(xw,ywp,'r','LineWidth',2)  %
subplot(2,1,2)
cdfplot(Y)
hold on
ywc=gevcdf(xw,wcoef(1),wcoef(2),wcoef(3)); %
plot(xw,ywc,'r','LineWidth',2)
legend('Empirical','Theoretical','Location','SE')
title('Empirical and Theoretical Distribution of Wind Speed')
xlabel('Y - Wind Speed (m/s)');
ylabel('F(y)');
hold off
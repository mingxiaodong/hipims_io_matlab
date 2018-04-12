%% import paired observations
clc,clear
% create fake observations
% p = copularnd('Gumbel',2,1000); X = p(:,1)*50+100; Y = p(:,2)*2;
% X = rand(1000,1); Y = rand(1000,1);
load('/Users/b4042552/Google Drive/MyResearch/London/matlab/MatchedDailyData')
X = Matched_FS.Values(:,1);
Y = Matched_FS.Values(:,3);
%% view the data and calculate chi dependence and independence value
u = (0.01:0.02:0.95)';
[chi_u,chiB_u]= QuantileAsymptoticalDependence(u,X,Y);
%%
figure(1)
subplot(1,2,1)
scatter(X,Y,'.')
xlabel('X')
ylabel('Y')
axis square
subplot(1,2,2)
%%
plot(u, [chi_u,chiB_u])
lgd = legend({'$\chi(u)$','$\bar{\chi}(u)$'},'Interpreter','Latex');
lgd.FontSize = 12;
xlim([0 1])
ylim([0 1])
ylabel('\chi(u)')
xlabel('u')
axis square
%% significance test (without considering the seasonality)
% permutation test
% resample without replacement
sampleTimes = 199;
chi_u_re_all = zeros(length(u),sampleTimes);
n = length(X);
for i=1:sampleTimes %repeat sampling sampleTimes
    % randomly select from X and Y
    Y_resample = datasample(Y,n,'Replace',false);
    [chi_u_re,~]= QuantileAsymptoticalDependence(u,X,Y_resample);
    chi_u_re_all(:,i) = chi_u_re;
end
chi_u_re_all_sort = sort(chi_u_re_all,2,'descend');
% the 10th largest value is the 5% significance level
chi_u_sigV = chi_u_re_all_sort(:,round(sampleTimes*0.05));
figure(2) % 
plot(u,chi_u,'--b',u,chi_u_sigV,'--r')
legend({'\chi(u)','5% significance level'},'Location','best')
title('significance test for \chi(u)')
%% Confident Interval
% resample with replacement
sampleTimes = 200;
chi_u_re_all = zeros(length(u),sampleTimes);
n = length(X);
for i=1:sampleTimes %repeat sampling sampleTimes
    % randomly select from X and Y
    np = floor(rand(n,1)*n+1);
    X_resample = X(np);
    Y_resample = Y(np);
%     X_resample = datasample(X,n);
%     Y_resample = datasample(Y,n);
    [chi_u_re,~]= QuantileAsymptoticalDependence(u,X_resample,Y_resample);
    chi_u_re_all(:,i) = chi_u_re;
end
chi_u_re_all_sort = sort(chi_u_re_all,2,'ascend');
chi_u_CI05 = chi_u_re_all_sort(:,round(sampleTimes*0.025));
chi_u_CI95 = chi_u_re_all_sort(:,round(sampleTimes*0.975));
figure(3) 
plot(u,chi_u,'r') %
hold on
% plot(u,chi_u_re_all,'b.')
plot(u,[chi_u_CI05 chi_u_CI95],'b:')
hold off
legend({'\chi','95% CI'})
% xlim([0.3 1])
% ylim([0 1])
xlabel('X');ylabel('Y')
title('Confidence Interval of \chi(u)')
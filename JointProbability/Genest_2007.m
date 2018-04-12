%% Genest and Farve-2007-Journal of Hydrologic Engineering-
%  Everything You Always Wanted to Know about Copula Modeling but Were Afraid to Ask
% table 1: Learning Data Set
clear,clc
X = [-2.224 -1.538 -0.807 0.024 0.052 1.324]';
Y = [0.431 1.035 0.586 1.465 1.115 -1.847]';
Z = exp(X);
T = exp(3*Y);
% figure 1
figure(1)
subplot(2,2,1)
scatter(X,Y,'o')
title('Scatter plot of X,Y')
axis square
subplot(2,2,2)
scatter(Z,T,'o')
title('Scatter plot of Z,T')
axis square
% ranks
R = tiedrank(X);
S = tiedrank(Y);
R_star = tiedrank(Z);
S_star = tiedrank(T);
% figure 2
subplot(2,2,3)
scatter(R,S,'o')
title('Scatter plot of R,S')
axis square
subplot(2,2,4)
scatter(R_star,S_star,'o')
title('Scatter plot of R^*,S^*')
axis square
%% empirical copula
u = 0.1:0.1:0.9;
v = u;
[um, vm]= meshgrid(u,v);
U = [um(:),vm(:)];
Cn_uv = copulaEcdf([X,Y],U);
figure
surf(um,vm,reshape(Cn_uv,size(um)))
xlabel('Fx'); ylabel('Gy')
%% Measuring Dependence
% Pearson's Correlation coefficient
Rho_P = corr(X,Y,'Type','Pearson');
% Spearman’s rho
Rho_S = corr(X,Y,'Type','Spearman');
% Kendall's Tau
Rho_K = corr(X,Y,'Type','Kendall');
n = length(X);
% test statistic for Spearman’s rho
% rho~N(0,1/(n-1))
test_Z = (Rho_S-0)/sqrt(1/(n-1));
p = normcdf(test_Z,0,1);
%% Chi plot
[Chi,Lambda,K] = Chi_K_Plot_Fisher(X,Y);
figure
subplot(1,2,1)
plot(Lambda,Chi,'x')
xline = -1:0.1:1;
xlim([-1 1])
ylim([-1 1])
xlabel('\lambda'); ylabel('\chi')
axis square
subplot(1,2,2)
curvX = 0:0.01:1;
curvY = curvX-curvX.*log(curvX);
plot(K(:,1),K(:,2),'x',0:0.5:1,0:0.5:1,'k-',curvX,curvY,'k-')
xlim([0 1])
ylim([0 1])
axis square
xlabel('W_{i:n}'); ylabel('H_{(i)}')
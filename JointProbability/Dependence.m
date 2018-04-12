% calculate the dependence
%% *****load data and uniform the data
clc,clear
load('CountyEvents0921.mat', 'CountyEvent');
X = [CountyEvent{1}.MaxRain]';
Y = [CountyEvent{1}.MaxWind]';
n = numel(X);
[RankX,~] = tiedrank(X); % tied rank of X&Y 
[RankY,~] = tiedrank(Y); 
U = RankX/(n+1); % Uniform of X 
V = RankY/(n+1); % Uniform of Y
%% plot (X,Y); (Rank X, Rank Y) or (U,V)
td= 0.8;
ind = U>=td|V>=td;
figure
subplot(1,2,1)
scatter(X,Y), xlabel('X'), ylabel('Y'), axis square
hold on; scatter(X(ind),Y(ind),'*r'), hold off, title('Observations')
subplot(1,2,2)
scatter(U,V), xlabel('U'), ylabel('V'), axis square, title('CDF')

%% Dependence measure Chi (Coles,Heffernan,and Tawn, 1999)
td = 0:0.01:1; %threshold level (0,1) 
chi_u = dependence_chi(td,U,V); % my function
% X_t = prctile(X,t);
figure; plot(td,chi_u); xlabel('Threshold td'); ylabel('Dependence \chi')
%% calculate Spearman's Rho

%% plot empirical copula
[u,v] = meshgrid(0:0.05:1,0:0.05:1);
U = [u(:),v(:)];
Fxy = copulaEcdf([X Y],U); Fxy = reshape(Fxy,size(u)); %invoke Ecdf copula function
XX = quantile(X,u(:)); XX = reshape(XX,size(u));
YY = quantile(Y,v(:)); YY = reshape(YY,size(u));
figure
subplot(2,2,1)
surf(XX,YY,Fxy), title('Ecopula for variable value')
subplot(2,2,2)
surf(u,v,Fxy), title('Ecopula for variable cdf')
subplot(2,2,3)
scatter(X,Y), title('Scatters for variable')
subplot(2,2,4)
scatter(X,Y), title('Scatters for variable')
%%
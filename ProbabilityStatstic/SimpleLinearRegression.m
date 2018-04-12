%% Single Linear Regression
clc
load data
X = xy(:,1);
Y = xy(:,2);
Xbar = mean(X); 
Ybar = mean(Y);
n = numel(X);
anno_str = cell(1,3);

Sxy = sum((X-Xbar).*(Y-Ybar));
Sxx = sum((X-Xbar).^2);
Syy = sum((Y-Ybar).^2);
r_Correlation = Sxy/(sqrt(Sxx)*sqrt(Syy));
anno_str{3} = ['The sample Linear Correlation Coefficient: r_{corr} = ' num2str(r_Correlation)];

B1_hat = Sxy/Sxx;
B0_hat = Ybar-B1_hat*Xbar;
RegressionFormula = ['y = ' num2str(B0_hat,'%.3f') num2str(B1_hat,'%+.3f') 'x'];
anno_str{1} = ['Regression Equation: ' RegressionFormula];

Yhat = B0_hat + B1_hat*X;
SSE = sum((Y-Yhat).^2);
SST = sum((Y-Ybar).^2);
SSR = SST - SSE;
r_square = 1-SSE/SST;
anno_str{2} = ['R square: r^2 = ' num2str(r_square)];

s_square = SSE/(n-2);
s = sqrt(s_square);
s_B1hat = s/sqrt(Sxx);
al = 0.05;
t_al = tinv(1-al/2,n-2);
CI_B1_hat = [B1_hat-t_al*s_B1hat B1_hat + t_al*s_B1hat];
%%*plot
plot(X,Y,'*',X,Yhat,'-')
h_anno = annotation('textbox','String',anno_str);
h_anno.Position = [0.2 0.12 0.6 0.12];
h_anno.FontSize = 12;
h_anno.EdgeColor = 'none';
%% *bivariate probability distribution
clear,clc
X = randn(1000,1);
Y = randn(1000,1);
m1 = mean(X); s1 = std(X); 
m2 = mean(Y); s2 = std(Y);
C = cov(X,Y);
r = C(2)/(s1*s2);
x = linspace(min(X),max(X),20);
y = linspace(min(X),max(X),20);
[xx, yy] = meshgrid(x,y);
x = xx(:);
y = yy(:);
e_power = -( ((x-m1)/s1).^2 - 2*r*(x-m1).*(y-s2)/(s1*s2) + ((y-m2)/s2).^2 )/(2*(1-r)^2);
fxy = exp(e_power)/(2*pi*s1*s2*sqrt(1-r^2));
surf(xx,yy,reshape(fxy,size(xx)))
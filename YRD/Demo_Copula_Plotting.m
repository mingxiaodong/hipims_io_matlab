% This is an example to show how to illustrate the joint probability
% distribution of two random variables X and Y.
% Let's assume that X obeys weibull distribution with scale para a1= 1 and shape para b1 = 5;
% and Y obeys gamma distribution with scale para a2 = 2, b2 = 2; 
% and the joint distribution of X and Y obeys Gumbel Copula with alpha = 1.2 

% initial value of distribution parameters
a1 = 1; b1 = 1.5; % weibull parameters
a2 = 2; b2 = 2; % gamma parameters
alpha = 1.2; % Gumbel copula parameter

% generate X and Y values
x = 0:0.2:4;
y = 0:20; 

% generate the cdf value of x and y
Fx = wblcdf(x, a1,b1); % you can replace it with the distrution function of your data
Fy = gamcdf(y, a2,b2);

% plot the cdf figure
figure
subplot(1,2,1)
plot(x,Fx); xlabel('x'); ylabel('Fx'); title('Marginal Distributin of X')
subplot(1,2,2)
plot(y,Fy); xlabel('y'); ylabel('Fy'); title('Marginal Distributin of Y')

% generate the Fxy points corresponding to x and y values
[XX,YY] = meshgrid(x,y); %generate a mesh grid of X and Y for plotting
u = wblcdf(XX(:), a1,b1); % the cdf of all the values in XX
v = gamcdf(YY(:), a2,b2); % the cdf of all the values in YY
Fxy = copulacdf('Gumbel',[u v],alpha); % The joint cdf of XX and YY
Fxy = reshape(Fxy,size(XX));

% plot the joint cdf
figure; surf(XX,YY,Fxy)
xlabel('x'); ylabel('y'); zlabel('Fxy')
title('Joint distribution of X and Y')
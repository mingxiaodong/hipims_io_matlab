% solve equation: finite diference solver
% y' = x-y;
% If y0=1, i.e. y=1 when x=0, calculate the values of yi for the  
% domain 0 ?x? 3 with dx = 0.5
%% Euler Method
L = 3;
dx = 0.1;
jj = round(L/dx); %number of cells
x = 0:dx:L;
y = x+nan;
y(1) = 1;
for i=1:jj
    y(i+1) = y(i)+dx*(x(i)-y(i));
end
plot(x,y)
%% Improved Euler Method
L = 3;
figure
hold on
for dx = [0.25,0.5,1]
jj = round(L/dx); %number of cells
x = 0:dx:L;
y = x+nan;
y_p = y;
y(1) = 1;
for i=1:jj
    % Predictor
    xip05 = (x(i)+x(i+1))/2;
    y_p(i+1) = y(i)+dx*(xip05-y(i));
    % Corrector
    y_ip05 = (y(i)+y_p(i+1))/2;
    y(i+1) = y(i)+dx*(xip05-y_ip05);    
end
plot(x,y)
end
plot(x,y_p,'-o')
ylim([0 2.5])
grid on
hold off
%%
Qc = 5;
a = 45; b = 100; m = 0.56;
t = 0:5:50;
Q = [10 10 40 75 100 90 60 35 25 20 10];


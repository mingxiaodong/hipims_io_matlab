% QUICKEST scheme solution to the advection-diffusion equation
% with exponential initial distribution
% clear workspace
clear
clc

%define variables
xmin = 0;
xmax = 10;
u = 20;             % flow celerity
D = 0;              % diffusion coefficient
N = 100;            % no. nodes - 1
t = 0;              % initial time
tmax = 0.25;         % max value of time

% discretization of domain
dx = (xmax - xmin)/N;
x = xmin - dx : dx : xmax + dx;  
% x values go from xmin to xmax in steps of dx with two external nodes
dt = 0.001;         % time step

% set initial conditions
c0 = exp(-2*(x-2.5).^2);
c = c0;             % set initial conditions for c
cnew = c0;          % set initial conditions for c at next time step

% set up time loop
nsteps = tmax/dt;

% set numerical parameters
r = D*dt/(dx*dx);    % Fourier number
Cr = u*dt/dx;        % Courant number

% plot initial distribution
plot(x,c0,'g-');     
axis([xmin xmax -0.1 1.1] )
xlabel('x','fontsize', 14)
ylabel('c(x,t)','fontsize',14)
title(sprintf('Courant no = %1.3f r =%1.3f',Cr,r));
hold on

% set coefficients for numerical scheme
coef1 = r*(1 - Cr) - Cr * (Cr^2 - 3*Cr + 2)/6
coef2 = 1 - r*(2 - 3*Cr) + 0.5 * Cr * (Cr^2 - 2*Cr - 1)
coef3 = r*(1 - 3*Cr) - 0.5* Cr * (Cr^2 - Cr - 2)
coef4 = Cr * ( r + (Cr^2 -1)/6 )

for n = 1: nsteps
        
    % for loop for first order upwind scheme
    for j = 3 : N + 2
        cnew(j) = coef1*c(j+1) + coef2*c(j) + coef3*c(j-1) + coef4*c(j-2);
    end
    cnew(1) = cnew(3);
    cnew(2) = cnew(4);
    
    % update t and c
    t = t + dt;
    c = cnew;

    % plot initial distribution
    plot(x,c0,'g-');
    hold on
    % plot solution
    plot(x,c,'bo-','markerfacecolor','b','markersize',1);
    title(sprintf('Time = %1.3f  Courant no = %1.3f r = %1.3f',t,Cr,r));
    axis([xmin xmax -0.1 1.1] )
    xlabel('x','fontsize', 14)
    ylabel('c(x,t)','fontsize',14)
    hold off
    
    shg
    pause(dt);

end
% plot initial distribution
plot(x,c0,'g-');     
axis([xmin xmax -0.1 1.1] )
xlabel('x','fontsize', 14)
ylabel('c(x,t)','fontsize',14)
hold on
plot(x,c,'bo-','markerfacecolor','b','markersize',1);
title(sprintf('Time = %1.3f  Courant no = %1.3f r =%1.3f',t,Cr,r));
hold off
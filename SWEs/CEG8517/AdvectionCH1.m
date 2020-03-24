% first order upwind finite difference solution to 1d advection equation
% with exponential initial distribution
% clear workspace
clear
clc

%define variables
xmin = 0;
xmax = 1;
N = 100;            % no. nodes - 1
dt = 0.008;         % time step
t = 0;              % initial time
tmax = 0.5;         % max value of time
u = 1;              % flow velocity (celerity)

% discretization of domain
dx = (xmax - xmin)/N;
x = xmin - dx : dx : xmax + dx;  
% x values go from xmin to xmax in steps of dx with two 'ghost nodes'

% set initial conditions
c0 = exp(-200*(x-0.4).^2);
% plot(x,c0);         % plot initial distribution
c = c0;             % set initial conditions for c
cnew = c0;          % set initial conditions for c at next time step

% set up time loop
nsteps = tmax/dt;
Cr = u*dt/dx;           % Courant number

for n = 1: nsteps
        
    % for loop for first order upwind scheme
    for i = 2 : N + 2
        cnew(i) = c(i) - Cr*(c(i) - c(i-1));
        if (cnew(i)<1E-5)
            cnew(i) = 0;
        end
    end
    cnew(1) = cnew(3);
    cnew(N+3) = cnew(N+1);
    
    % update t and c
    t = t + dt;
    c = cnew;
    
    % calculate exact solution
    exact = exp(-200*(x - 0.4 - u*t).^2);

    % plot solution
    plot(x, exact, 'r-');
    hold on
    plot(x,c0,'g-');
    plot(x,c,'bo-','markerfacecolor','b','markersize',2);
    hold off
    axis([xmin xmax -0.1 1.1] )
    xlabel('x','fontsize', 14)
    ylabel('c(x,t)','fontsize',14)
    title(sprintf('Time = %1.3f  Courant no = %1.3f',t,Cr ))
    shg
    pause(dt);

end

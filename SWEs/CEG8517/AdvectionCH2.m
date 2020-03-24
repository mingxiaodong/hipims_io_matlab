% first order explicit upwind finite difference solution to 1d advection equation
% with exponential initial distribution
% clear workspace
clear
clc

%define variables
xmin = 0;
xmax = 1;
N = 100;            % no. nodes - 1
dt = 0.01;         % time step
t = 0;              % initial time
tmax = 5;         % max value of time
u = 50;              % flow velocity (celerity)

% discretization of domain
dx = (xmax - xmin)/N;
x = xmin - dx : dx : xmax + dx;  
% x values go from xmin to xmax in steps of dx with two 'ghost nodes'

% set initial conditions
% c0 = x*0;
% ind = x<0.2&x>=0;
% c0(ind) = 0;
% ind = x<0.3&x>=0.2;
% c0(ind) = 10*(x(ind)-0.2);
% ind = x<0.3&x>=0.2;
% c0(ind) = 10*(0.4-x(ind));
c0 = exp(-200*(x-20).^2);

plot(x,c0,'g-');     % plot initial distribution
format1 = 'Time = %1.3f  Courant no = %1.3f';
hold on

c = c0;             % set initial conditions for c
cnew = c0;          % set initial conditions for c at next time step

% set up time loop
nsteps = tmax/dt;
Cr = u*dt/dx;           % Courant number
tplotdiff = 10*dt;
tplot = tplotdiff;

    % calculate boundary conditions
    c(1) = c(3);
    c(N+3) = c(N+1);

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

    if (t>tplot)
        
        % plot exact solution
         plot(x, exact, 'r-');
        % plot numerical solution
        plot(x,c,'bo-','markerfacecolor','b','markersize',2);
        tplot = tplot+tplotdiff;
    end
    shg
    pause(dt);

end
    hold off
    axis([xmin xmax -0.1 1.1] )
    xlabel('x','fontsize', 14)
    ylabel('c(x,t)','fontsize',14)
    title(sprintf(format1,t,Cr ))
     
    % open output file
    outfile = 'AdvecExplicit.dat';
    file1 = fopen(outfile,'w');
    % write data to file
    format2='dx=%1.3f dt=%1.3f u=%1.3f Cr=%1.3f \r\n';
    fprintf(file1,format2,dx,dt,u,Cr);
    fprintf(file1,'t=%1.3f \r\n',t);
    fprintf(file1,'x c_init c_exact c_final \r\n');
    format2 = '%1.3f %1.3f %1.3f %1.3f \r\n'; 
    % for loop for output
    for i = 2 : N + 2
        fprintf(file1,format2,x(i), c0(i), exact(i),c(i));
    end
    fclose(file1);

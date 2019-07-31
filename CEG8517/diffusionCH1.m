% classic explicit finite difference solution to 1d diffusion equation
% with exponential initial distribution

% clear workspace
clear
clc

%define variables
xmin = 0;
xmax = 50;          % 50m domain
N = 100;            % no. of cells = no. nodes - 1
D = 10;             % diffusion coefficient
dt = 0.01;        % time step 
t = 0;              % initial time
tmax = 5;         % max value of time


% discretization of domain
dx = (xmax - xmin)/N;
x = xmin - dx : dx : xmax + dx;  
% x values go from xmin to xmax in steps of dx with two external nodes

% set initial conditions
c0 = exp(-0.1*(x-25).^2);

% plot initial distribution
plot(x,c0,'g-');     
axis([xmin xmax -0.1 1.1] )
xlabel('x','fontsize', 14)
ylabel('c(x,t)','fontsize',14)
hold on

c = c0;             % set initial conditions for c
cnew = c0;          % set initial conditions for c at next time step

% set up time loop
nsteps = tmax/dt;
r = D*dt/(dx*dx);   % Fourier number

% set number of time steps between plots
tplotnum = 50*dt;
tplot = tplotnum;

for n = 1: nsteps
        
    % for loop for classic explicit scheme
    for j = 2 : N + 2
        cnew(j) = r*c(j+1) + (1 - 2*r)*c(j) + r*c(j-1);
    end
    
    % boundary conditions
    cnew(1) = 2*r*c(2)+(1-2*r)*c(1);%c(1);
    cnew(2) = (1-r)*c(2)+r*c(1);
    cnew(N+3) = c(N+3);
    
    % update t and c
    t = t + dt;
    c = cnew;
    
    % calculate exact solution
    % exact = exp(-200*(x - 0.25 - u*t).^2);

    if (t>=tplot)
        
        % plot exact solution
        % plot(x, exact, 'r-');
        % plot numerical solution
        plot(x,c,'bo-','markerfacecolor','b','markersize',1);
        title(sprintf('Time = %1.3f r = %1.3f',t,r ))
        tplot = tplot + tplotnum;
    end
    shg
    pause(dt);

end
plot(x,c,'r-');
title(sprintf('Time = %1.3f r = %1.3f',t,r ))
hold off

    % open output file
    outfile = 'DiffusionExplicit.dat';
    file1 = fopen(outfile,'w');
    % write data to file
    format2='dx=%1.3f dt=%1.3f D=%1.3f r=%1.3f \r\n';
    fprintf(file1,format2,dx,dt,D,r);
    fprintf(file1,'t=%1.3f \r\n',t);
    fprintf(file1,'x c_init c_final \r\n');
    format2 = '%1.3f %1.3f %1.3f \r\n'; 
    % for loop for output
    for i = 2 : N + 2
        fprintf(file1,format2,x(i), c0(i), c(i));
    end
    fclose(file1);

%% 1D Shallow Water Equation with the HLL Riemann Solver
% Purpose: to simulate 1D free surface flow in an open channel
% * Godunove-type scheme is used** 1st order, HLL Riemann Solver
clear,clc
% * Constants
lx = 1400; %size of domain
M = 100; %number of cells
g = 9.81; %gravity acceleration

% * Flow variables
zb = zeros(1,M,'double'); % bed elevation above datum (m)
z = zeros(1,M,'double'); % zeta: water surface level (m), h = z-zb
qx = zeros(1,M,'double'); %flow flux qx=uh
zn = zeros(1,M,'double'); % convey new value of z(i)
qxn = zeros(1,M,'double'); % convey new value of qx(i)
x = zeros(1,M); % location of the domain

% * to initialize some varialbes
dx = lx/M; %cellsize
dt = 0.2; %time step
tout = 50; %terminating time
t = 0.0; %real time
for i = 1:M 
    %**set the domain attributes
    x(i) = (i-0.5)*dx;    
        %if x(i)>8 && x(i)<12 zbloc = 0.2-0.05*(x(i)-10)^2;
        %else zbloc = 0;
        %end
    zbloc = 0;%10+40*x(i)/lx+10*sin(pi*(4*x(i)/lx-1/2)); %set the bed level  
    zb(i) = zbloc;
    if i>=80 && i<=100
        z(i) = 80;  % water level through the channel
    else
        z(i) = 65;
    end
    qx(i) = 0; % flux through the channel    
end

% * plot innitial conditions
figure(1)
    % **bed profile
plot(x,zb,'k-','LineWidth',3)
fill([x(1) x x(M)],[0 z 0],[0 0.5 1])
hold on

plot(x,z,'b-')
fill([x x(M)],[zb 0],[0.3 0.3 0.3])
hold on


title({num2str(t)})
xlabel('location in channel (m)')
ylabel('water level (m)')

axis([0 lx 0 100])
hold off
disp('Now the figure shows the initial conditions. Please press Enter to continue')
pause

% * main programe
while t<tout
    
    t = t+dt;
    %**calculate new values
    for i = 1:M               
        %***set the bed elevation(zb),water surface level(z), and flux(qx) in east and west face          
            %**** east face based the two interfacing cells i and i+1    
        zb_L = zb(i); z_L = z(i); qx_L = qx(i);            
        if i == M %cell in the east boundary,need a ghost right cell
            [zb_R,z_R,qx_R] = bound(1,zb_L,z_L,qx_L); %bound function location: eastern = 1, western = 2;
        else
            zb_R = zb(i+1); z_R = z(i+1); qx_R = qx(i+1); %true right cell
        end       
        zbfe = (zb_L+zb_R)/2.0; % bed level in east face                            
            %****HLL Riemann Solver Function, Fe is a vector consisting of two elements: Fe1 and Fe2
        [Fe1,Fe2] = F_cal1(zbfe,z_L,z_R,qx_L,qx_R); 
        
            %**** west face based the two interfacing cells i-1 and i
        zb_R = zb(i); z_R = z(i); qx_R = qx(i);
        if i == 1 %cell in the west boundary,need a ghost left cell
           [zb_L,z_L,qx_L] = bound(2,zb_R,z_R,qx_R); %bound function location: eastern = 1, western = 2;
        else
           zb_L = zb(i-1); z_L = z(i-1); qx_L = qx(i-1); %true left cell
        end
        zbfw = (zb_L+zb_R)/2.0; % bed level in west face
            %****HLL Riemann Solver Function
        [Fw1,Fw2] = F_cal1(zbfw,z_L,z_R,qx_L,qx_R); 
                
        
        %***source term****
        s_dx = (zbfe-zbfw)/dx; % slope in cell i 
        s1 = 0;
        s2 = -g*z(i)*s_dx;        
        
        %***changing gap between U^n and U^n+1
        %consist of two element: dz and dqx
        dU1 = -dt*((Fe1-Fw1)/dx-s1);
        dU2 = -dt*((Fe2-Fw2)/dx-s2);
        
        %***calculate U^n+1 based on U^n and the gap
        zn(i) = z(i)+dU1; % z and qx are the elements of U^n
        qxn(i) = qx(i)+dU2;
        
    end
    
    if isnan(zn(i))
        disp('unstable')
        break
    end
    
    %if max(abs(z-zn))<=0.01
    %    disp(t)
    %    disp('no water level change any more')
    %    break
    %end...%
    
    %**update
    z = zn; 
    qx = qxn;
    
    %**plot results
    figure(1)
        %***bed profile
plot(x,zb,'k-','LineWidth',3)
fill([x(1) x x(M)],[0 z 0],[0 0.5 1])
hold on

plot(x,z,'b-')
fill([x x(M)],[zb 0],[0.3 0.3 0.3])
hold on
    
    title(['t = ' num2str(t)])
    xlabel('location in channel (m)')
    ylabel('water level (m)')
    axis([0 lx 0 100])
    pause(0.001)
    hold off
    
end


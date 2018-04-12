% 2nd order Godunov-Type Scheme for 1D SWEs
% zm qxm for U_star, the predicted one
% Purpsoe: to store flow data used by different program units.
clear,clc
global M dx g zb z qx zm qxm

% * Constants
lx = 14000; %length of domain
M = 50; %number of cells
g = 9.81; %gravity acceleration

% * Flow variables
zb = zeros(1,M,'double'); % bed elevation above datum (m)
z = zeros(1,M,'double'); % zeta: water surface level (m), h = z-zb
qx = zeros(1,M,'double'); %flow flux qx=uh
zn = zeros(1,M,'double'); % convey new value of z(i)
qxn = zeros(1,M,'double'); % convey new value of qx(i)
x = zeros(1,M); %? location of the domain

% * to initialize some varialbes
dx = lx/M; %cellsize
dt = 1; %time step
tout = 5000; %terminating time
t = 0.0; %real time

for i = 1:M 
    %**set the domain attributes
    x(i) = (i-0.5)*dx; 
    
    zb(i) = 10+40*x(i)/lx+10*sin(pi*(4*x(i)/lx-1/2)); %set the bed level
    z(i) = 65; % water level is 2 through out the channel
    qx(i) = 0; % flux is 4.42 through out the channel    
end

%% * plot innitial conditions
figure(1)
    % **bed profile
plot(x,zb,'b-','LineWidth',3)
hold on

plot(x,z,'ro','MarkerSize',4)
hold on

title({num2str(t)})
xlabel('itx\rm (m)')
ylabel('it\eta\rm (m)')
axis([0 14000 0 100])
hold off

%% * main programe
% allocate matrices for intermeidatevariables
zm=zeros(1,M,'double');
qxm=zeros(1,M,'double');
kmz=zeros(1,M,'double');
kmqx=zeros(1,M,'double');

while t<tout
    
t = t+dt;
%** 1st step of RK2 -- calculate U* and K(U^n)
for i=1:M
% east face
     [zb_L,z_L,qx_L] = Face(1,1,i);
     if i==M %boundary conditions
         [zb_R,z_R,qx_R] = bound_face(1,zb_L,z_L,qx_L);
     else
         [zb_R,z_R,qx_R] = Face(1,2,i+1);
     end
     zbfe = (zb_L+zb_R)/2;
     [Fe1,Fe2] = F_cal2(zbfe,z_L,z_R,qx_L,qx_R); %calculate flux vector east
% west face
     [zb_R,z_R,qx_R] = Face(1,2,i);
     if i==1 %boundary conditions
         [zb_L,z_L,qx_L] = bound_face(1,zb_R,z_R,qx_R);
     else
         [zb_L,z_L,qx_L] = Face(1,1,i-1);
     end
     zbfw = (zb_L+zb_R)/2;
     [Fw1,Fw2] = F_cal2(zbfw,z_L,z_R,qx_L,qx_R); %calculate flux vector west
%***source term
        sox = (zbfe-zbfw)/dx; % slope
        s1 = 0;
        s2 = -g*z(i)*sox;        
%***new value
        dzdt = (Fe1-Fw1)/dx-s1;
        dqxdt = (Fe2-Fw2)/dx-s2;
        
        zn(i) = z(i)-dt*dzdt; %where is z(i)come from
        qxn(i) = qx(i)-dt*dqxdt;
        kmz(i) = dzdt;
        kmqx(i) = dqxdt;
        zm(i) = z(i)-dt*dzdt;
        qxm(i) = qx(i)-dt*dqxdt;
end

%** 2nd step of RK2 -- calculate U^n+1 and K(U*)
for i=1:M
% east face
     [zb_L,z_L,qx_L] = Face(2,1,i);
     if i==M %boundary conditions
         [zb_R,z_R,qx_R] = bound_face(1,zb_L,z_L,qx_L);
     else
         [zb_R,z_R,qx_R] = Face(2,2,i+1);
     end
     zbfe = (zb_L+zb_R)/2;
     [Fe1,Fe2] = F_cal1(zbfe,z_L,z_R,qx_L,qx_R); %calculate flux vector east
% west face
     [zb_R,z_R,qx_R] = Face(2,2,i);
     if i==1 %boundary conditions
         [zb_L,z_L,qx_L] = bound_face(2,zb_R,z_R,qx_R);
     else
         [zb_L,z_L,qx_L] = Face(2,1,i-1);
     end
     zbfw = (zb_L+zb_R)/2;
     [Fw1,Fw2] = F_cal1(zbfw,z_L,z_R,qx_L,qx_R); %calculate flux vector west
%***source term
        sox = (zbfe-zbfw)/dx; % slope
        s1 = 0;
        s2 = -g*zm(i)*sox;        
%***new value
        dzdt = (Fe1-Fw1)/dx-s1;
        dqxdt = (Fe2-Fw2)/dx-s2;
        
        zn(i) = z(i)-dt*(kmz(i)+dzdt)/2.0; %where is z(i)come from
        qxn(i) = qx(i)-dt*(kmqx(i)+dqxdt)/2.0;
    
end    
%**update
    z = zn; 
    qx = qxn;
    
    %**plot results
    figure(1)
        %***bed profile
    plot(x,zb,'b-','LineWidth',3)
    hold on

    plot(x,z,'ro','MarkerSize',4)
    hold on
    
    title({num2str(t)})
    xlabel('itx\rm (m)')
    ylabel('it\eta\rm (m)')
    axis([0 14000 0 100])
    pause(0.0001)
    hold off
end

    


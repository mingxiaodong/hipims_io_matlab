%               Code for Introduction Lecture                  %
%
%
%--------------------------------------------------------------------%
%                1D shallow water equation solver                    %
% Purpose: tosimulate 1D                an open channel.% %% %% %% %%
%
%--------------------------------------------------------------------
clc,clear
%** constants
lx = 25.0;
M = 100;
g = 9.81;
%** flow variables
zb = zeros(1,M,'double'); % bed elevation (DEM)
z = zb; zn = zb; % water level (eta)
qx = zb; qxn = zb;% velocity
x = zeros(1,M); % cell location along channel
% Purpose: initialize variables
dx = lx/M; % cellsize
dt = 0.001; % timestep
tout = 50;
t = 0;
x = dx/2:dx:dx*(length(x)-1)+dx/2;
% initial condition
z = z+2;
qx = qx+4.42;
zb = 0.2-0.05*(x-10).^2;
zb(x<=8|x>=12) = 0;
figure(1); plot(x,[zb;z]);
ylim([0 2.2]) 
legend('river bed','water level')
%% -------------------------------------------------------------------
% main programe
figure(2)
while t<tout
    t = t+dt;
    %* calculate new values
    for i=1:M
        %** east face
        zb_L = zb(i); z_L = z(i); qx_L = qx(i);
            % whether in the right bound
        if i==M
            [zb_R,z_R,qx_R] = bound(1,zb_L,z_L,qx_L);
        else
            zb_R = zb(i+1); z_R = z(i+1); qx_R = qx(i+1);
        end
        [Fe1,Fe2] = F_cal(zb_L,z_L,qx_L,zb_R,z_R,qx_R);
        zbfe = (zb_L+zb_R)/2;
        
        %** west face
        zb_R = zb(i); z_R = z(i); qx_R = qx(i);
            % whether in the left bound
        if i==1
            [zb_L,z_L,qx_L] = bound(2,zb_R,z_R,qx_R);
        else
            zb_L = zb(i-1); z_L = z(i-1); qx_L = qx(i-1);
        end
        [Fw1,Fw2] = F_cal(zb_L,z_L,qx_L,zb_R,z_R,qx_R);
        zbfw = (zb_L+zb_R)/2;
        
        %** source term
        s1 = 0;
        sox = (zbfe-zbfw)/dx;
        s2 = -g*z(i)*sox;
        
        %** new value
        dzdt = (Fe1-Fw1)/dx - s1;
        dqxdt = (Fe2-Fw2)/dx - s2;
        
        zn(i) = z(i) - dt*dzdt;
        qxn(i) = qx(i) - dt*dqxdt;
    end
    %* update
    z = zn; qx = qxn;
    %* plot
    plot(x,[zb;z]);
    ylim([0 2.2])
    title(['t = ' num2str(t) 's'])
    pause(0.1)
%     disp(t)
end
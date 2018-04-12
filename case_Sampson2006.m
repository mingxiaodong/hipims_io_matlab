%% generate a one dimension DEM and initial condition
clear,clc
domainL = 10000;
cellsize = 100;
domainW = cellsize*3;
x11 = -domainL*0.5;
y11 = domainW/2-0.5*cellsize;
dx = cellsize;
dy = -cellsize;
r_dem =makerefmat(x11, y11, dx, dy);

h0 = 10;
a = 3000;
x_vec = -domainL/2:cellsize:domainL/2;
b_vec = h0*((x_vec/a).^2);
t = 0;
B = 5;
g = 9.8;
% s = sqrt(8*g*h0/a^2/2);
% eta_theo_fun =@(x_f,t_f) h0 - (a^2*B^2*s^2/(8*g^2*h0))*cos(2*s*t_f) - B^2/(4*g)- B*s*x_f/g;
tau = 0.00;
s = sqrt(8*g*h0/a^2-tau^2)/2;
eta_theo_fun =@(x_f,t_f) h0 + ...
    (a^2*B^2*exp(-tau*t_f)/(8*g^2*h0)).*...
    (-s*tau*sin(2*s*t_f) + (tau^2/4-s^2)*cos(2*s*t_f))...
    - B^2*exp(-tau*t_f)/(4*g)...
    - (exp(-tau*t_f/2)/g)*(B*s*cos(s*t_f)+0.5*tau*B*sin(s*t_f)).*x_f;


eta_theo = eta_theo_fun(x_vec,t);
h_vec = eta_theo-b_vec; h_vec(h_vec<0) = 0;
z_dem = repmat(b_vec,[3,1]);
h_grid = repmat(h_vec,[3,1]);
gauges_pos = [x_vec',x_vec'*0+y11-cellsize];
% figure
% plot(x_vec,b_vec,x_vec,eta_theo)
% ylim([min(b_vec), max(b_vec)])
% legend('bed','\eta value')
%% write file
fileName = 'G:\cudaSWEsSolver\WRRTest1\dem.txt';
arcgridwrite(fileName,z_dem,r_dem)
fileName = 'G:\cudaSWEsSolver\WRRTest1\h_initial.asc';
arcgridwrite(fileName,h_grid,r_dem)
fileName = 'G:\cudaSWEsSolver\WRRTest1\gauges.txt';
dlmwrite(fileName,gauges_pos)
% figure;
% hold on
% mapshow(z_dem,r_dem,'DisplayType','surface');
% scatter3(gauges_pos(:,1),gauges_pos(:,2),gauges_pos(:,2)*0+20,'r')
% axis square
% hold off
%% prepare function parameters
addpath C:\Users\b4042552\Dropbox\Matlab\GeoClasses
% addpath /Users/b4042552/Dropbox/Matlab/GeoClasses
CaseFolder = 'G:\cudaSWEsSolver\WRRTest1';
cd(CaseFolder)
h_Eta = 'h'; % write h file or eta file
DEMName = 'dem.txt';
[Z,R] = arcgridread(DEMName);
copyfile(DEMName,[CaseFolder '/input/mesh/DEM.txt'])
h_initial = arcgridread('h_initial.asc');
times_setup = [0 1200 250 1500];
dlmwrite([CaseFolder '/input/times_setup.dat'],times_setup,' ')
BoundCode = [2 0 0;...
             2 2 0];
gauges_position = dlmread('gauges.txt');
%%field set
FieldSetupV2(CaseFolder, Z, R,'h_Eta',h_Eta,...
        'BoundCode',BoundCode,... % boundary code
        'initial_hE',h_initial,... % boundary source
        'manning',0,... % manning value
        'GaugeCoor',gauges_position,... %gauges_position
        'WriteAllFiles','true');
%% read output data
z_h1 = arcgridread('h_250.asc');
z_h2 = arcgridread('h_500.asc');
eta_gauges = dlmread('eta_gauges.dat');
%% plot one time
figure;
subplot(2,2,1)
t = 250;
eta_theo = eta_theo_fun(x_vec,t);
h_line = plot(x_vec,[b_vec;b_vec+z_h1(2,:);eta_theo]);
h_line(1).LineWidth = 2; h_line(1).Color = 'k';
h_line(3).Color = 'b';
legend('bed','\eta simulated','\eta theory')
title('t = 250s')
xlim([-3000,0])
ylim([9.6 13])

subplot(2,2,2)
t = 250;
eta_theo = eta_theo_fun(x_vec,t);
h_line = plot(x_vec,[b_vec;b_vec+z_h1(2,:);eta_theo]);
h_line(1).LineWidth = 2; h_line(1).Color = 'k';
h_line(3).Color = 'b';
title('t = 250s')
xlim([-2000,-1500])
ylim([11 11.8])

subplot(2,2,3)
t = 500;
eta_theo = eta_theo_fun(x_vec,t);
h_line = plot(x_vec,[b_vec;b_vec+z_h2(2,:);eta_theo]);
h_line(1).LineWidth = 2; h_line(1).Color = 'k';
h_line(3).Color = 'b';
title('t = 500s')
xlim([0,3000])
ylim([9 15])

subplot(2,2,4)
t = 500;
eta_theo = eta_theo_fun(x_vec,t);
h_line = plot(x_vec,[b_vec;b_vec+z_h2(2,:);eta_theo]);
h_line(1).LineWidth = 2; h_line(1).Color = 'k';
h_line(3).Color = 'b';
title('t = 500s')
xlim([1500,2000])
ylim([11.5 13])
%% plot whole time
timeSeries = eta_gauges(:,1);
eta_all = eta_gauges(:,2:end);
figure
for i=1:length(timeSeries)
    t = timeSeries(i);
    eta_theo = eta_theo_fun(x_vec,t);
    eta_theo(eta_theo<b_vec) = nan;
    h_line = plot(x_vec,[b_vec;eta_all(i,:);eta_theo]);
    h_line(1).LineWidth = 2; h_line(1).Color = 'k';
    h_line(2).LineStyle = '--';
    h_line(3).Color = 'b';
    legend('bed','\eta simulated','\eta theory')
    title(['t = ' num2str(t) 's'])
    ylim([0 30])
    pause(0.5)
end
%%
figure
t_long = 1000;
for i = 1:length(t_long)
    t = t_long(i);
    eta_theo = eta_theo_fun(x_vec,t);
    eta_theo(eta_theo<b_vec) = nan;
    h_line = plot(x_vec,[b_vec;eta_theo]);
    h_line(1).LineWidth = 2; h_line(1).Color = 'k';
    legend('bed','\eta theory')
    title(['t = ' num2str(t) 's'])
    ylim([0 30])
    pause(0.5)
end
axis([-3500 3500 0 14])
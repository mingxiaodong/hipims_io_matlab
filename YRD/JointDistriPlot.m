%% 21-08-2015 plot the joint probability distribution of multi-hazard for each county
% clayton copula
clf,clc,clear
load CountyCopula0921
Codes = [310230;321182;330681;330481]; %county code
code = Codes(3);
%Chongming 310230
%Yangzhong 321182
%Zhuji 330681
%Haining 330481 innormal value
fid = find([CountyCopula.code]==code); % find the county's serial number in var 'CountyCopula'
CName = CountyCopula(fid).cname;
CPara = CountyCopula(fid).al;
RainMDcoef = CountyCopula(fid).RainGEV;
WindMDcoef = CountyCopula(fid).WindGEV;
CountyName = CountyCopula(fid).ename;
Mt = 40/length(CountyCopula(fid).Rain); %number of events

NoP = 40; %number of points for X and Y axes
Xm = max(CountyCopula(fid).Rain);
Ym = max(CountyCopula(fid).Wind);
X = linspace(RainMin2,Xm,NoP);
Y = linspace(WindMin,Ym,NoP);

[XX,YY] = meshgrid(X,Y);
Fx = gevcdf(XX,RainMDcoef(1),RainMDcoef(2),RainMDcoef(3));
Fy = gevcdf(YY,WindMDcoef(1),WindMDcoef(2),WindMDcoef(3));
F = copulacdf(CName,[Fx(:) Fy(:)],CPara);

Fxy = reshape(F,size(XX));
Txy = Mt./(1-Fx-Fy+Fxy);

%plot and save

%%*******plot Fxy and its contour*****
subplot(1,2,1);
surf(XX,YY,Fxy)
zlabel('Fxy','fontsize',20)
xlabel('Rainfall (mm)','fontsize',20)
ylabel('Wind speed (m/s)','fontsize',20)
title(CountyName ,'fontsize',20)
ax = gca;
ax.FontSize = 20;
axis([RainMin2,Xm,WindMin,Ym,0,1])

subplot(1,2,2);
v1 = [0.1 0.3 0.5 0.7 0.9]';
c = contour(X,Y,Fxy,v1,'linewidth',2);
%clabel(c,h,'LabelSpacing',300)
clabel(c,'fontsize',20)
title(CountyName ,'fontsize',20)
xlabel('Rainfall (mm)','fontsize',20)
ylabel('Wind Speed (m/s)','fontsize',20)
ax = gca;
ax.FontSize = 20;
axis([RainMin2,Xm,WindMin,Ym])
%save figure


%%*******plot Txy and its contour*****
figure;
subplot(1,2,1);
surf(XX,YY,Txy)
zlabel('Txy','fontsize',20)
xlabel('Rainfall (mm)','fontsize',20)
ylabel('Wind speed (m/s)','FontSize',20)
title(CountyName ,'fontsize',20)
ax = gca;
ax.FontSize = 20;
axis([RainMin2,Xm,WindMin,Ym,0,max(Txy(:))])

subplot(1,2,2);
v2 = [5 10 20 50 100]';
c = contour(X,Y,Txy,v2,'linewidth',2);
clabel(c,'fontsize',20)
title(CountyName ,'fontsize',20)
xlabel('Rainfall (mm)','fontsize',20)
ylabel('Wind Speed (m/s)','fontsize',20)
ax = gca;
ax.FontSize = 20;
axis([RainMin2,Xm,WindMin,Ym])

%% plot marginal distributin of rainfall and wind speed for Copula function
load CountyCopula
i = 113;
X = CountyCopula(i).Rain;
Y = CountyCopula(i).Wind;
rcoef = CountyCopula(i).RainGEV;
wcoef = CountyCopula(i).WindGEV;
CountyName = CountyCopula(i).ename;
figure
subplot(2,1,1)
cdfplot(X)
xlabel('X - Rainfall (mm)','FontSize',30);
ylabel('Fx','FontSize',30);
hold on
xr=min(X):max(X);
yrc=gevcdf(xr,rcoef(1),rcoef(2),rcoef(3));   % GEV cumulative distribution
plot(xr,yrc,'r','LineWidth',4)
legend('Empirical','Theoretical','Location','SE')
title('Rainfall distribution')
ax = gca;
ax.FontSize = 30;
hold off
subplot(2,1,2)
cdfplot(Y)
hold on
xw=min(Y):0.1:max(Y);    %
ywc=gevcdf(xw,wcoef(1),wcoef(2),wcoef(3)); %
plot(xw,ywc,'r','LineWidth',4)
legend('Empirical','Theoretical','Location','SE')
xlabel('Y - Wind Speed (m/s)','FontSize',30);
ylabel('Fy','FontSize',30);
title('Wind speed distribution')
ax = gca;
ax.FontSize = 30;
hold off
%% return period for set hazard magnitude
clc,clear
load CountyCopula0921
Rrank = [50 100 250]; % I, II, III,
Wrank = [10.8 17.2 24.5]; % A, B, C
% to store return period of each county for different hazard combination
ReInterval = zeros(140,9);
[X,Y] = meshgrid(Rrank,Wrank);
for i=1:140
    rcoef = CountyCopula(i).RainGEV;
    wcoef = CountyCopula(i).WindGEV;
    Fx = gevcdf(X(:),rcoef(1),rcoef(2),rcoef(3));
    Fy = gevcdf(Y(:),wcoef(1),wcoef(2),wcoef(3));
    Fxy = copulacdf(CountyCopula(i).cname,[Fx, Fy],CountyCopula(i).al);
    Mt = 40/length(CountyCopula(i).Rain);
    Txy = Mt./(1-Fx-Fy+Fxy);
    ReInterval(i,:) = Txy';
end
ReInterval = [[CountyCopula.code]',ReInterval];

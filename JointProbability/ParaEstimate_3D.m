%%load data with three variables
clear,clc
load('CountyEvents0921.mat','CountyEvent')
i = 139; % any county
% Source data of X, Y, Z
Xs = [CountyEvent{i,1}.TotalRain]';
Ys = [CountyEvent{i,1}.MaxRain]';
Zs = [CountyEvent{i,1}.MaxWind]';
[fx,x] = ecdf(Xs);
[fy,y] = ecdf(Ys);
[fz,z] = ecdf(Zs);
XYZcdf = zeros(length(Xs),3);

for n=1:length(Xs)
    XYZcdf(n,1) = min(fx(x == Xs(n)));
    XYZcdf(n,2) = min(fy(y == Ys(n)));
    XYZcdf(n,3) = min(fz(z == Zs(n)));
end


%%
Tau = corr(Ys,Zs,'Type','Kendall');
%Rho_s = corr(Ys,Zs,'Type','Spearman');
Cfamily = 'Frank';
alpha = copulaparam(Cfamily,Tau);
n = 1000;
U = copularnd(Cfamily,alpha,n);
figure
plot(U(:,1),U(:,2),'.')
hold on
plot(XYZcdf(:,2),XYZcdf(:,3),'xr')
title(['Kendall''s {\it\tau} = ', num2str(Tau),'  Copula Type: ', Cfamily])
xlabel('U1')
ylabel('U2')
hold off
%%
CoefX = gevfit(Xs);
CoefY = gevfit(Ys);
CoefZ = gevfit(Zs);
u1s = gevcdf(Xs,CoefX(1),CoefX(2),CoefX(3));
u2s = gevcdf(Ys,CoefY(1),CoefY(2),CoefY(3));
u3s = gevcdf(Zs,CoefZ(1),CoefZ(2),CoefZ(3));
%%
syms f(si,u1,u2,u3)
f(si, u1, u2, u3) = exp(...
    -(...
    (-log(u1)).^si+(-log(u2)).^si+(-log(u3)).^si...
    ).^(1/si)...
    );
fs = 0;
for i=1:length(Xs)
    fs =fs + log(subs(f,{u1,u2,u3},{u1s(i),u2s(i),u3s(i)}));
end

fd = diff(fs,si)/length(Xs);

%%
tic
theta = solve(fd,si);
theta = double(theta);
toc
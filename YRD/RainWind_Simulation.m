%% Rain and Wind Events Simulation
% created on 15-01-2016 by Ming
clear,clc
RWrand = dlmread('RW_Simulations.txt');
load('CountyEvents153010.mat','CountyEvent')
load('FT_Result153010.mat','FT_Result')
% to store the simulated value of Rain and Wind
RainWind_Rand = cell(140,2);
n = 1000;
for i = 1:140
    X = [CountyEvent{i}.MaxRain];
    Y = [CountyEvent{i}.MaxWind];
    Fx = RWrand(n*(i-1)+1:n*i,2);
    Fy = RWrand(n*(i-1)+1:n*i,3);
    
    [Rain,xi] = ksdensity(X,Fx,'function','icdf');
    [Wind,yi] = ksdensity(Y,Fy,'function','icdf');   
    RainWind_F = [Rain,Wind,Fx,Fy]; %based on ksdensity
    RainWind_Rand{i,1} = RainWind_F;
    
    Rcof = FT_Result(i).RainGEV;
    Wcof = FT_Result(i).WindGEV;
    Rain = icdf('gev',Fx,Rcof(1),Rcof(2),Rcof(3));
    Wind = icdf('gev',Fy,Wcof(1),Wcof(2),Wcof(3));
    RainWind_F = [Rain,Wind,Fx,Fy]; %based on Gev
    RainWind_Rand{i,2} = RainWind_F;
end
%% compare simulations and observations
load('CountyInfor')
i = 140;
% simulated rain and wind based on ksdensity
Xvk = RainWind_Rand{i,1}(:,1);
Yvk = RainWind_Rand{i,1}(:,2);
% simulated rain and wind based on gev distribution
Xvg = RainWind_Rand{i,2}(:,1);
Yvg = RainWind_Rand{i,2}(:,2);
% orgional rain and wind
Xo = [CountyEvent{i}.MaxRain];
Yo = [CountyEvent{i}.MaxWind];
subplot(1,2,1)
scatter(Xvk,Yvk,'*')
hold on 
scatter(Xo,Yo,'r*')
hold off
title([CountyInfor(i).ename, '--Kernel'])
subplot(1,2,2)
scatter(Xvg,Yvg,'*')
hold on 
scatter(Xo,Yo,'r*')
hold off
title([CountyInfor(i).ename, '--Gev'])

%% 2015-8-13 copula fitting and selection for each county, based on county events data
clear,clc
load CountyInfor
load FT_Result0921
CountyCode = [CountyInfor.code];
CountyCopula = CountyInfor;

for i=1:140    
    CountyCopula(i).Rain = FT_Result(i).Rain;
    CountyCopula(i).Wind = FT_Result(i).Wind;
    CountyCopula(i).RainGEV = FT_Result(i).RainGEV;
    CountyCopula(i).WindGEV = FT_Result(i).WindGEV;
    X = CountyCopula(i).Rain;
    Y = CountyCopula(i).Wind;
    rcoef = CountyCopula(i).RainGEV;  
    wcoef = CountyCopula(i).WindGEV;
        
    U = gevcdf(X,rcoef(1),rcoef(2),rcoef(3));    %CDF value of rainfall
    V = gevcdf(Y,wcoef(1),wcoef(2),wcoef(3));    %CDF value of wind speed
    % parameter of copula
    al_gum = copulafit('Gumbel',[U(:), V(:)]);         %
    al_cla = copulafit('Clayton',[U(:), V(:)]);        %
    al_fra = copulafit('Frank',[U(:), V(:)]);          %

%%Comparison between different types of Copula
    % Empirical CDF of X and Y
    [fx, Xsort] = ecdf(X);
    [fy, Ysort] = ecdf(Y);
    RealL = min(length(Xsort),length(Ysort));
    Fx = fx(2:RealL);
    Fy = fy(2:RealL);
    
    % Define Empirical Copula
    EC = @(u,v)mean((Fx <= u).*(Fy <= v));    
    
    % Calculate the Empirical Copula function value of original sample's CDF
    CUV = zeros(size(Fx));
    for j=1:numel(Fx)            % 
        CUV(j) = EC(Fx(j),Fy(j));
    end
    
    % Calculate the Fitted Copula function value of original sample's CDF
    Cgum = copulacdf('Gumbel',[Fx(:), Fy(:)],al_gum);
    Ccla = copulacdf('Clayton',[Fx(:), Fy(:)],al_gum);
    Cfra = copulacdf('Frank',[Fx(:), Fy(:)],al_fra);

    % Calculate Euclidean Distance
    d2Gumbel = (CUV-Cgum)'*(CUV-Cgum);
    d2Clayton = (CUV-Ccla)'*(CUV-Ccla);
    d2Frank = (CUV-Cfra)'*(CUV-Cfra);

    % Compare the distance and select the optimum Copula
    distance = [d2Gumbel d2Clayton d2Frank];
    [mindistan,b] = min(distance); %#ok<*ASGLU> 
    if b==1
        cname='Gumbel';
        al=al_gum;
    elseif b==2
        cname='Clayton';
        al=al_cla;
    elseif b==3
        cname='Frank';
        al=al_fra;
    else
        cname='erro';
    end
    CountyCopula(i).cname = cname;
    CountyCopula(i).al = al;
    
    %disp('The optimum copula is ')
    %disp([cname ' copula'])
end
save CountyCopula0921 CountyCopula RainMin1 RainMin2 WindMin
%% figure Empirical copula
clc
i = 28;
X = CountyCopula(i).Rain;
Y = CountyCopula(i).Wind;
[fx, Xsort] = ecdf(X);                      
[fy, Ysort] = ecdf(Y);
Fx = fx(2:end);
Fy = fy(2:end);

%% Define Empirical Copula
EC = @(u,v)mean((Fx <= u).*(Fy <= v)); 
[Udata,Vdata] = meshgrid(linspace(0,1,31)); % create meshgrid ranging in [0,1]*[0,1]
CopulaEmpirical = zeros(size(Udata));
for j = 6:numel(Udata)
    CopulaEmpirical(j) = EC(Udata(j),Vdata(j));
end
surf(Udata,Vdata,CopulaEmpirical)
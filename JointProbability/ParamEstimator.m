function [theta,thetaCI,temp] = ParamEstimator(X,Y,typeStr)
%
% Fisher and Switzer-2001-Graphical Assessment of Dependence Is a Picture Worth 100 Tests?
if size(X)~=size(Y)
    error('X and Y should have the same size')
end
ConL = 0.95; % confidence level
Z_score = norminv(1-(1-ConL)/2);
W = X*0;
W_td = X*0; % w tilde
n = numel(X);
% calculate H
for i=1:n
    cntW = 0;
    cntWt = 0;
    for j=1:n
        if X(j)<=X(i)&& Y(j)<=Y(i)
            cntW = cntW+1;
        end
        if X(i)<=X(j)&& Y(i)<=Y(j)
            cntWt = cntWt+1;
        end
    end
    W(i) = cntW/n;
    W_td(i) = cntWt/n;
end
W_bar = mean(W);
%% estimate based on Kendall's Tau
if strcmp(typeStr,'Tau')
    KenTau = 4*(n/(n-1))*W_bar-(n+3)/(n-1);
    syms t
    gfun = 9/2*t;
    gfun_D = diff(gfun); gfun_D = double(gfun_D);
    S_sq = sum((W+W_td-2*W_bar).^2)/n;
    S_sqrt = sqrt(S_sq);
    N_sigma = 4*S_sqrt*abs(gfun_D)/(sqrt(n));
    theta = subs(gfun,KenTau); theta= double(theta);
    %% estimate based on Spearmans' Rho
    % theta = 3*rho {h(rho)}
else
    Rho = SpearmansRho([X Y]);
    syms rho
    hfun = 3*rho;
    hfun_D = diff(hfun); hfun_D = double(hfun_D);
    % calculate sigma_n
    R = tiedrank(X);
    S = tiedrank(Y);
    %     An = mean((R.*S)/(n+1)^2);
    %     Bn = mean((R.*S).^2/(n+1)^4);
    %     Ci = X*0; Di = Ci; Ei = Ci;
    %     for i=1:n
    %         R_i1n = R(i)/(n+1);
    %         S_i1n = S(i)/(n+1);
    %         Cj = X*0;
    %         Dj = X*0;
    %         Ej = X*0;
    %         for j = 1:n
    %             Dj(j) = S_i1n*S(j)/(n+1)*max(R_i1n, R(j)/(n+1));
    %             Ej(j) = R_i1n * R(j)/(n+1) * max(S_i1n, S(j)/(n+1) );
    %             for k=1:n
    %                 Ck = X*0;
    %                 if R(k)<=R(i)&&S(k)<=S(j)
    %                     Ck(k) = 1;
    %                 end
    %             end
    %             Cj(j) = mean(Ck);
    %         end
    %         Ci(i) = mean(Cj)*R_i1n*S_i1n;
    %         Di(i) = mean(Dj);
    %         Ei(i) = mean(Ej);
    %     end
    %     Cn = mean(Ci)+1/4-An;
    %     Dn = mean(Di); En = mean(Ei);
    An = mean((R.*S)/(n+1)^2);
    Bn = mean((R.*S).^2/(n+1)^4);
    Ci = R*0;
    Di = Ci; Ei = Ci;
    for i=1:n
        Cj = R*0;
        Dj = R*0;
        Ej = R*0;
        for j=1:n
            Ck = R*0;
            for k=1:n
                if (R(k)<=R(i))&&(S(k)<=S(j))
                    indF = 1;
                else
                    indF = 0;
                end
                Ck(k) = indF*R(i)*S(i)/(n+1)/(n+1);
            end
            Cj(j) = sum(Ck)/n;
            Dj(j) = S(i)*S(j)/(n+1)/(n+1)*max(R(i)/(n+1),R(j)/(n+1));
            Ej(j) = R(i)*R(j)/(n+1)/(n+1)*max(S(i)/(n+1),S(j)/(n+1));
        end
        Ci(i) = sum(Cj)/n;
        Di(i) = sum(Dj)/n;
        Ei(i) = sum(Ej)/n;
    end
    Cn = sum(Ci)/n+0.25-An;
    Dn = sum(Di)/n;
    En = sum(Ei)/n;
    S_sq = 144*(-9*An^2+Bn+2*Cn+2*Dn+2*En);
    S_sqrt = sqrt(S_sq);
    N_sigma = S_sqrt*abs(hfun_D)/(sqrt(n));
    theta = subs(hfun,Rho); theta= double(theta);
end
thetaCI = [theta-Z_score*N_sigma,theta+Z_score*N_sigma];
% temp = [An Bn Cn Dn En];
temp = Z_score*N_sigma;%S_sqrt;%
end
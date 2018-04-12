function [chi_u,chiB_u,A_CI,B_CI,A_Sig,B_Sig]= QuantileAsymptoticalDependence(u,X,Y,N)
% The function is to calculate the chi dependence of bivariable sample(X,Y)
% t is the threshold, can be a scarlar or vector in (0 1)
% U and V are the Uniform Margins of X and Y respectively
% chi_u_l: lower bound of chi
% Reference: Coles,Heffernan,and Tawn, 2000. Dependence Measures for Extreme Value Analyses
% N = 200; % resample times
if nargin==3
    N = 200;
end
ind = isnan(X)|isnan(Y);
X(ind)=[];
Y(ind)=[];
[chi_u,chiB_u] = chi_calculation(X,Y,u);
% for confidence interval
if nargout>2
    n = numel(X); % sample size
    al = 0.95; % confident level
    A_re = nan(length(u),N); %chi_u
    B_re = A_re; %chiB_u
    for i=1:N
        Ind = floor(rand(n,1)*n+1);
        %     Ind = datasample((1:n)',n);
        X_re = X(Ind);%datasample(X,n,'Replace',true);
        Y_re = Y(Ind);%datasample(Y,n,'Replace',true);
        [A,B] = chi_calculation(X_re,Y_re,u);
        A_re(:,i) = A;
        B_re(:,i) = B;
    end
    A_re = sort(A_re,2,'ascend');
    B_re = sort(B_re,2,'ascend');
    A_CI = [A_re(:,round(N*(1-al)/2)) A_re(:,round(N*(1-(1-al)/2)))];
    B_CI = [B_re(:,round(N*(1-al)/2)) B_re(:,round(N*(1-(1-al)/2)))];
    
    sigV = 0.05;
    A_re = nan(length(u),N); %chi_u
    B_re = A_re; %chiB_u
    X_re = X;
    for i=1:N
        %     Ind = floor(rand(n,1)*n+1);
        Ind = datasample((1:n)',n,'Replace',false);
        Y_re = Y(Ind);%datasample(Y,n,'Replace',true);
        [A,B] = chi_calculation(X_re,Y_re,u);
        A_re(:,i) = A;
        B_re(:,i) = B;
    end
    A_re = sort(A_re,2,'descend');
    B_re = sort(B_re,2,'descend');
    A_Sig = A_re(:,round(N*sigV));%
    B_Sig = B_re(:,round(N*sigV));
end
end

function [A,B] = chi_calculation(X,Y,u)
% A: chi_u
% B: chiBar_u
td_X = quantile(X,u);
td_Y = quantile(Y,u);
A = zeros(size(u));
B = zeros(size(u));

for i = 1:numel(u)
    Fxy = sum(X<=td_X(i)&Y<=td_Y(i))/length(X);
    Fx  = sum(X<=td_X(i))/length(X);
    Fy  = sum(Y<=td_Y(i))/length(Y);
    F_xy = sum(X>td_X(i)&Y>td_Y(i))/length(X);
    A(i) = 2 - log(Fxy)/log((Fy+Fx)/2);%2 - log(Fxy)/log(Fx); %    
%     F_xy = 1-Fx-Fy+Fxy;
%     disp(F_xy)
%     if F_xy == 0
%         F_xy = 0.5/length(X);
%     end 
    B(i) = 2*log(1-Fx)/(log(F_xy))-1;
end


end
function [Chi,Lambda] = Chi_K_Plot_Fisher(X,Y)
%
% Fisher and Switzer-2001-Graphical Assessment of Dependence Is a Picture Worth 100 Tests?
if size(X)~=size(Y)
   error('X and Y should have the same size') 
end
H = X*0;
F = X*0;
G = X*0;
% W = X*0;
n = numel(X);
% calculate H
for i=1:n
    cntH = 0;
    cntF = 0;
    cntG = 0;
    for j=1:n
        if j~=i
            if X(j)<=X(i)&& Y(j)<=Y(i)
                cntH = cntH+1;
            end
            if X(j)<=X(i)
                cntF = cntF+1;
            end
            if Y(j)<=Y(i)
                cntG = cntG+1;
            end
        end
    end
    H(i) = cntH/(n-1);
    F(i) = cntF/(n-1);
    G(i) = cntG/(n-1);
end
Chi = (H-F.*G)./(F.*(1-F).*G.*(1-G)).^0.5;
S = sign((F-0.5).*(G-0.5));
Lambda = 4*S.*max([(F-0.5).^2,(G-0.5).^2],[],2);

% for i=1:n
%     funIntegral = @(w) w.*(-log(w)).*(w-w.*log(w)).^(i-1).*(1-w+w.*log(w)).^(n-i);
%     W(i) = n*nchoosek(n-1,i-1)*integral(funIntegral,0,1);
% end
% 
% K = [W H];
end
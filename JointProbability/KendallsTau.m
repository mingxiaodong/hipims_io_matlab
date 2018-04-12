function KenTau = KendallsTau(X,Y)
% Compute the Kendall's tau coefficient of concordance of the matrix X.
% formula of calculation is from Genest-2007-Everything you ...
% E.g. KenTau = KendallsTau(RankMatrix)
%
% Input:
%           X must be a N-by-K matrix, N is the number of
%           "candidate" and the K is the number of "judge"
% Outputs:
%           W = Kendall's coefficient of concordance
%
% Edited by Xiaodong Ming, 2016/9/14
%==========================================================================
if numel(X)~=numel(Y)
    error('input vectors should have the same length')
end
n = numel(X);
W = X*0;
for i=1:n
    cntW = 0;
    for j=1:n
            if X(j)<=X(i)&& Y(j)<=Y(i)
                cntW = cntW+1;
            end
    end
    W(i) = cntW/n;
end
W_bar = mean(W);
KenTau = 4*(n/(n-1))*W_bar-(n+3)/(n-1);
end
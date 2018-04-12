X = Xsample;
U = [(0.01:0.01:0.99)',(0.01:0.01:0.99)'];
Cn = copulaEcdf(X,U);

function Cn = copulaEcdf(X,U)
% calculate empirical copula of multiple variables
% X is a matrix whose column number represents the number of varibales
% U is a matrix representing cumulative probabilities for each variable
% Cn is the empirical cdf corresponding to p
n = size(X,1); % number of observations for each variable
d = size(X,2); % number of variables (dimensions)
if size(U,2)~=d
    error('dimensions of input parameter do not match ')
end
[~,R] = sort(X);
U_hat = R./(n+1);
indfun = @(A) sum( sum(U_hat<=repmat(A,[n,1]), 2) == d )/n; % A is the input para of indfun
U_table = table(U);
Cn_table = rowfun(indfun,U_table);
Cn = table2array(Cn_table);
end



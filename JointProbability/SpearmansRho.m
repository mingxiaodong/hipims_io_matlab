function SpmRho = SpearmansRho(X)
% Compute the Spearman's Rho coefficient of the matrix X.
% E.g. SpmRho = SpearmansRho(X)
%
% Input:
%           X must be a N-by-K matrix, N is the number of
%           "candidate" and the K is the number of "judge"
% Outputs:
%           W = Spearman's Rho coefficient
%
% Edited by Xiaodong Ming, 2016/9/14
%==========================================================================
n = size(X,1);
R = tiedrank(X); % get the ranks for data set of X 
R_bar = mean(R); % mean rank
RandkDeviation = R-repmat(R_bar,[n,1]); %Rank deviation from mean rand value
SST_rank = sum(prod(RandkDeviation,2)); %total sum of square
RDSquareSum = sum(RandkDeviation.^2); 
RDSS_prod = prod(RDSquareSum); 
SpmRho = SST_rank/sqrt(RDSS_prod);
end
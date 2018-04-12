function [x3points, y3points,l3points,w3points] = Nearest3BankPoints(bankPoints_xy,bankPoints_lw,bedPoints_xy)
%% find the nearest three bank points and return there coordinates
% 1. find the nearest points
% 2. find the two adjacent points
% bankPoints_xy is a two column vector represents the coordinates of the
% points on the bank
% bedPoints_xy is a two column vector, coordinates of the targeting points 
numBankPoints = length(bankPoints_xy);
SquareDiff_vec = @(A,B) (repmat(A,[1 length(B)])-repmat(B',[length(A) 1])).^2;
% Euclidean Distance between bed Points and bank Points
% row: poly points column: bank points
disEuc_BankCross = (SquareDiff_vec(bedPoints_xy(:,1),bankPoints_xy(:,1))...
                    +SquareDiff_vec(bedPoints_xy(:,2),bankPoints_xy(:,2)))...
                    .^0.5;
[~, min_ind] = min(disEuc_BankCross,[],2); 
min_ind(min_ind==1)=2; min_ind(min_ind==numBankPoints)=numBankPoints-1;
GetThreeCoords = @(inds,vecCoords) [vecCoords(inds-1) vecCoords(inds) vecCoords(inds+1)];
x3points = GetThreeCoords(min_ind,bankPoints_xy(:,1));
y3points = GetThreeCoords(min_ind,bankPoints_xy(:,2));
l3points = GetThreeCoords(min_ind,bankPoints_lw(:,1));
w3points = GetThreeCoords(min_ind,bankPoints_lw(:,2));

end
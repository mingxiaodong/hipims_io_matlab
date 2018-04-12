function bank1Points_lw = RelativeCoordsBankPoints(bank1Points,w_Value)
%% bank1Points_lw = RelativeCoordsBankPoints(bank1Points,w_Value)
% build the ralative coordinates for river bank points
% L coordinate: length along the river bank (ratio)
% W coordinate: width across the river bank (ratio)
% H coordinate: height of the river bank    (absolute)
    % function to calculate the distance between two neigbour points in one river bank
Distance_InterPoints = @(A)((A(1:end-1,1)-A(2:end,1)).^2+(A(1:end-1,2)-A(2:end,2)).^2).^0.5;
%%*****for bank 1
interDis_Bank1Points = Distance_InterPoints(bank1Points);
cumuDist_Bank1Points = cumsum(interDis_Bank1Points);
% L, W, and H coordinate of river bank 1
coor_L_Bank1Points = [0; cumuDist_Bank1Points/sum(interDis_Bank1Points)];
coor_W_Bank1Points = coor_L_Bank1Points*0+w_Value;
bank1Points_lw = [coor_L_Bank1Points,coor_W_Bank1Points];
end
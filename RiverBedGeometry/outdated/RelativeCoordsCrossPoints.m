function cross1Points_lw = RelativeCoordsCrossPoints(cross1Points,l_Value)
%%build the ralative coordinates for river bank points
% L coordinate: length along the river bank (ratio)
% W coordinate: width across the river bank (ratio)
    % function to calculate the distance between two neigbour points in one river bank
cross1Points = cross1Points(:,1:2);
Distance_InterPoints = @(A)((A(1:end-1,1)-A(2:end,1)).^2+(A(1:end-1,2)-A(2:end,2)).^2).^0.5;
%%*****for bank 1
interDis_Cross1Points = Distance_InterPoints(cross1Points);
cumuDist_Cross1Points = cumsum(interDis_Cross1Points);
% L, W, and H coordinate of river bank 1
coor_W_Cross1Points = [0; cumuDist_Cross1Points/sum(interDis_Cross1Points)];
coor_L_Cross1Points = coor_W_Cross1Points*0+l_Value;
cross1Points_lw = [coor_L_Cross1Points,coor_W_Cross1Points];
end
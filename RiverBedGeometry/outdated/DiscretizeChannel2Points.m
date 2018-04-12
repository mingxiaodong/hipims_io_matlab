function [X_grid,Y_grid,Z_grid] = DiscretizeChannel2Points(crossLine0,crossLine1,resolution)
%DiscretizeChannel2Points Discretize quadrangular/triangular channel
%section to feature points based on the elevation of the points in the 
%start and end cross lines of the river channel section. The banks line of
%this channel section must be straight!
%   crossLine1,crossLine2: x,y,z coordinates and z value of cross line points
%   resolution: resolution of the output points
%   Created by Xiaodong Ming on 2017-10-30.
%   See also DiscretizeChannel2Quadrangle
    DistanceTwoPoints = @(xy1,xy2) (sum((xy1-xy2).^2)).^0.5;
    crossLine0_D = DistanceTwoPoints(crossLine0(1,1:2),crossLine0(end,1:2));
    crossLine1_D = DistanceTwoPoints(crossLine1(1,1:2),crossLine1(end,1:2));
    bankLine0_D = DistanceTwoPoints(crossLine0(1,1:2),crossLine1(1,1:2));
    bankLine1_D = DistanceTwoPoints(crossLine0(end,1:2),crossLine1(end,1:2));
    n_bank  = round(max( bankLine0_D, bankLine1_D)/resolution)+1;
    n_cross = round(max(crossLine0_D,crossLine1_D)/resolution)+1;
    W = linspace(0,1, n_cross)';
    L = linspace(0,1,  n_bank)';
    X_grid = nan(n_cross,n_bank);
    Y_grid = X_grid;
    Z_grid = X_grid;
    crossLines0_W = RelativeCoordsCrossPoints(crossLine0,0);
    crossLines1_W = RelativeCoordsCrossPoints(crossLine0,1);
    crossLine0_mW = interp1(crossLines0_W(:,2),crossLine0(:,3),W);
    crossLine1_mW = interp1(crossLines1_W(:,2),crossLine1(:,3),W);
    lineraInterp = @(xyz0,xyz1,r) xyz0+r*(xyz1-xyz0);
    crossLine0_XYZ = lineraInterp(crossLine0(  1,:),crossLine0(end,:),W);
    crossLine1_XYZ = lineraInterp(crossLine1(  1,:),crossLine1(end,:),W);
    crossLine0_XYZ(:,3) = crossLine0_mW;
    crossLine1_XYZ(:,3) = crossLine1_mW;

    for i=1:n_cross
        bankLineXYZ  = lineraInterp(crossLine0_XYZ(i,:),crossLine1_XYZ(i,:),L);
        X_grid(i,:) = bankLineXYZ(:,1)';
        Y_grid(i,:) = bankLineXYZ(:,2)';
        Z_grid(i,:) = bankLineXYZ(:,3)';
    end

end


function z_new = LevelBound(z_dem,r_dem,boundFrame,levelValue,boundExpWidth)
%z_new = LEVELBOUND(z_dem,r_dem,boundInd,levelValue,boundExpWidth)    Level a part of 
%   outline bound defined by boundFrame to a same elevation value (levelValue) 
%   boundFrame is a n*4 numeric matrix, when it is empty[], the whole
%   outline bound will be levelled.
%   boundExpWidth is the number of grids expand outside of the boundary,
%   the default value is 0.
if nargin==3
    boundExpWidth = 1;
end
[~,Bound_Cell_ID] = gen_VB_ID(z_dem,r_dem,boundFrame);
z_new = z_dem;
z_ind = zeros(size(z_dem)); %index of grids to be levelled
if isempty(boundFrame)
    z_ind(Bound_Cell_ID==1) = 1; % outline bound
else
    z_ind(Bound_Cell_ID>1) = 1; % other bounds apart from outline bound
end
if boundExpWidth>1
    z_ind = rasterCal(z_ind,boundExpWidth-1,'max');
end
z_ind(Bound_Cell_ID==0) = 0;
z_new(z_ind>=1) = levelValue;
end
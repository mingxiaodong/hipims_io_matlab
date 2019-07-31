function h = MapVelocity(Ugrid,Vgrid,R,varargin)
% draw velocity map, return the Quiver object
% Syntax
%   MapVelocity(Ugrid,Vgrid,R)
%   MapVelocity(Ugrid,Vgrid,R,resolution)
%   MapVelocity(Ugrid,Vgrid,R,ax)
%   MapVelocity(Ugrid,Vgrid,R,resolution,ax)
% Ugrid: velocity grid in X direction
% Vgrid: velocity grid in Y direction
% R: reference matrix of the velocity grids 
% resolution: spatial resolution to resample velocity grids
% ax: axis to plot map
% See Also: ResampleRaster,Raster2FeaturePoints
% Created by Xiaodong Ming
resolution=R(2);
ax = gca;        
disp(varargin)
if nargin==4
    if isnumeric(varargin{1})
        resolution = varargin{1};
    elseif isobject(varargin{1})
        ax = varargin{1};
    end
elseif nargin==5
    ind1 = cellfun(@isnumeric,varargin);
    resolution = varargin{ind1};
    ind2 = cellfun(@isobject,varargin);
    ax = varargin{ind2};
end
if resolution~=R(2)
    [Ugrid,~] = ResampleRaster(Ugrid,R,resolution);
    [Vgrid,R] = ResampleRaster(Vgrid,R,resolution);
end
[X,Y] = Raster2FeaturePoints(Ugrid,R);
Y = Y-R(2);
Ugrid(abs(Ugrid)<0.0001&abs(Ugrid)>4)=nan;
Vgrid(abs(Vgrid)<0.0001&abs(Vgrid)>4)=nan;

h = quiver(ax,X,Y,Ugrid,Vgrid,3,'k');

end
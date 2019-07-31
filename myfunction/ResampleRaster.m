function [z_New,r_new] = ResampleRaster(Z,R,varargin)
% Resample a raster grid to a new spatial resolution
% [z_New,r_new] = ResampleRaster(Z,R,newResolution) resample a raster matrix
% [z_New,r_new] = ResampleRaster(Z,R,interpMethod)
% [z_New,r_new] = ResampleRaster(Z,R,newResolution,interpMethod)
% newResolution: 
% interpMethod: default linear
% function Raster2FeaturePoints is required in this function
interpMethod = 'linear';
newResolution = [];
z_outsize = [];
r_output = [];
if length(varargin)==1
    if isnumeric(varargin{1})
        newResolution = varargin{1};
    else
        interpMethod = varargin{1};
    end
elseif length(varargin)==2
    if isnumeric(varargin{1})&& isnumeric(varargin{2})
        z_outsize = varargin{1};
        r_output = varargin{2};
        newResolution = r_output(2);
    elseif isnumeric(varargin{1})
        newResolution = varargin{1}; interpMethod = varargin{2};
    else
        newResolution = varargin{2}; interpMethod = varargin{1};
    end
elseif length(varargin)==3
    z_outsize = varargin{1};
    r_output = varargin{2};
    interpMethod = varargin{3};
else
    error('wrong number of arguments')
end
dx = R(2,1); dy = R(1,2);
dx_new = newResolution; dy_new = -newResolution;
[m,n] = size(Z);
x11 = R(3,1)+dx;
y11 = R(3,2)+dy;
x = x11:dx:x11+dx*(n-1);
y = y11:dy:y11+dy*(m-1);
x_new = x11:dx_new:x(end);
y_new = y11:dy_new:y(end);
y = fliplr(y);
[X,Y] = ndgrid(x,y);
% X = X'; Y = Y';
Z = Z';
F = griddedInterpolant(X,Y,Z,interpMethod);%default: 'linear',%'ExtrapolationMethod','nearest');
if isempty(r_output)
    r_new = makerefmat(x11,y11,dx_new,dy_new);
    [X_new,Y_new] = ndgrid(x_new,y_new);
else
    r_new = r_output;
    [X_new,Y_new] = Raster2FeaturePoints(z_outsize,r_output);
    X_new = X_new'; Y_new = Y_new';
end
z_New = F(X_new,Y_new);
z_New = z_New';
z_New = flipud(z_New);
end

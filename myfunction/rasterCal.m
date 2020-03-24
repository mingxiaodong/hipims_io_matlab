function z_out = rasterCal(Z,varargin)
% z_output = rasterCal(Z,calRange,method) raster calculation
% z_output = rasterCal(Z,calRange,method,nanflag) raster calculation
% calRange: number of neighbouring cells to be calculated
% method: min, max, mode, median, sum
% nanFlag: 'omitnan','includenan'
nanFlag = 'omitnan';
if isempty(varargin)
   calRange = 1;
   method = 'mean';
elseif length(varargin)==1
    if ischar(varargin{1})
        calRange = 1;
        method = varargin{1};
    else
        calRange = varargin{1};
        method = 'mean';
    end
elseif length(varargin)==2
    calRange = varargin{1};
    method = varargin{2};
else
    calRange = varargin{1};
    method = varargin{2};
    nanFlag = varargin{3};
end
Z_expand_X = [repmat(Z(:,1),1,calRange), Z , repmat(Z(:,end),1,calRange)];
Z_expand = [repmat(Z_expand_X(1,:),calRange,1);...
              Z_expand_X ;...
              repmat(Z_expand_X(end,:),calRange,1)];
[firstInd1,firstInd2] = meshgrid(1:(2+calRange),1:(2+calRange));
firstInd = [firstInd1(:),firstInd2(:)];
[m,n] = size(Z);
z_adjacent = zeros(m,n,length(firstInd));
for i = 1:length(firstInd)
    m0 = firstInd(i,1);
    n0 = firstInd(i,2);
    z0 = Z_expand(m0:m0+m-1,n0:n0+n-1);
    z_adjacent(:,:,i) = z0;
end
if strcmp(method,'min')||strcmp(method,'max')
    eval(['z0 = ' method '(z_adjacent,[],3,''' nanFlag ''');'])
elseif strcmp(method,'mode')
    eval(['z0 = ' method '(z_adjacent,3);'])
else
    eval(['z0 = ' method '(z_adjacent,3,''' nanFlag ''');'])
end
z_out = z0;
end
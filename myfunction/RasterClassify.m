function Z_classified = RasterClassify(Z,breakValues)
%RasterClassify Classify the raster value into categories based on
%break values
%   Detailed explanation goes here
Z_classified = nan(size(Z));
categoryValue = 1:length(breakValues)+1;
for i=1:length(breakValues)+1
    switch i
        case 1
            ind = Z<breakValues(i);
        case length(breakValues)+1
            ind = Z>=breakValues(i-1);
        otherwise
            ind = Z>=breakValues(i-1)&Z<breakValues(i);
    end
    Z_classified(ind) = categoryValue(i);
end

end


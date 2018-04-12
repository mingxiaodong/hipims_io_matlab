classdef DomainRaster
    properties
        Z % 
        R
        PixelFlag_Valiad
        PixelFlag_OutlineBound
        PixelID_Valiad %size: same as Z; values: 0,1,...,n-1 (n: number of valid pixels)
        PixelID_OutlineBound %size: same as Z; values: 0,1 (1 means outline bound pixels)
        PixelID_AllBound %pixel ID of all Boundaries(outline+Input/Output bound) 2,3,...n+1 (n means the number of IO boundary)
        Frame_IOBound % cell of two-column matrixes: X and Y coordiantes of frame points
    end
    
    methods
        function PixelID_Valiad = getPixelID_Valiad(obj)
            if isempty(obj.PixelFlag_Valiad)
                obj.PixelFlag_Valiad = ~isnan(obj.Z);
            end
            ValiadInd_Change = flipud(obj.PixelFlag_Valiad');
            PixelID_Valiad = nan(size(ValiadInd_Change));
            ValiadInd_Change = find(ValiadInd_Change==true);
            [orgID,newID] = sort(ValiadInd_Change,'ascend');
            PixelID_Valiad(orgID) = flipud(newID)-1;
            PixelID_Valiad = flipud(PixelID_Valiad)';
        end
        function PixelID_OutlineBound = getPixelID_OutlineBound(obj)
            Z_expand = [obj.Z(:,1), obj.Z, obj.Z(:,end)];
            Z_expand = [Z_expand(1,:); Z_expand; Z_expand(end,:)];
            Z_cal = Z_expand(1:end-2,2:end-1) + Z_expand(3:end,2:end-1)+...
                    Z_expand(2:end-1,1:end-2) + Z_expand(2:end-1,3:end);
            PixelFlag = isnan(Z_cal)&~isnan(obj.Z);
            PixelID_OutlineBound = zeros(size(obj.Z));
            PixelID_OutlineBound(PixelFlag)=1;
        end
        function PixelID_AllBound = getPixelID_AllBound(obj)
            [XX_Z,YY_Z] = Raster2FeaturePoints(obj.Z,obj.R);
            if isempty(obj.PixelID_OutlineBound)
               obj.PixelID_OutlineBound = getPixelID_OutlineBound(obj);
            end
            PixelID_AllBound = obj.PixelID_OutlineBound;
            if ~iscell(obj.Frame_IOBound)
                masks = {obj.Frame_IOBound};
            else
                masks = obj.Frame_IOBound;
            end
            for i=1:length(masks)
                mask = masks{i};
                X = mask(:,1); Y = mask(:,2);
                if size(mask,1) == 2
                    rectangle_xy = [min(X),min(Y);max(X),min(Y);...
                        max(X),max(Y);min(X),max(Y);min(X),min(Y)];
                    X = rectangle_xy(:,1);
                    Y = rectangle_xy(:,2);
                end
                in = inpolygon(XX_Z(:),YY_Z(:),X,Y);
                ind_outlineBound = obj.PixelID_OutlineBound==1;
                ind_IO_Bound = in&ind_outlineBound(:);
                PixelID_AllBound(ind_IO_Bound)=i+1;
            end
        end
    end
end   
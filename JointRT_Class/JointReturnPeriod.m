classdef JointReturnPeriod
    properties
        Xtable
        Ytable
        Ztable
        note
        th1 = [0.98 0.98 0.98]; %extreme threshold for X,Y,Z
        th2 = [0.9 0.9 0.9]; %extreme threshold for XY,XZ,YZ
        th3 = 0.85; %%extreme threshold for XYZ
        cNames2D = {'clayton','clayton','frank'};
        cParms2D = [0.0311 0.0835 0.0937];
        cName3D = 'frank';
        cParm3D = 1.12;
    end
    methods
        function RT1struct = univariateRT(XdataTable,thP)
            Xdata = XdataTable.Value;
            Ddata = XdataTable.Date;
            yearnum = numel(Ddata)/365.25;
            RT1struct.yearnum = yearnum;
            thV = quantile(Xdata,thP);
            ind = Xdata>=thV;
            Xe = Xdata(ind);
            De = Ddata(ind);
            PureEvent = ExtractDuplicateEvents(De,Xe,3);
            RT1struct.PureEvent = PureEvent;
            RT1struct.thV = thV;
            Xe_pure = PureEvent.Event;
            RT1struct.Mt = yearnum/height(PureEvent);            
            distName = 'Generalized Extreme Value';%'Gamma';%
            pd = fitdist(Xe_pure,distName);
            h = kstest(Xe_pure,'CDF',pd);
            RT1struct.pd = pd;
            RT1struct.kstest = h;
            h = chi2gof(Xe_pure,'CDF',pd);
            RT1struct.chi2test = h;
        end
    end
    methods(Static)
        function F = ArchimCopulaCDF(cname,u,alpha)
            claytoncdf = @(u1,u2,alpha) (1 +(u1.^(-alpha)-1+u2.^(-alpha)- 1)).^(-1./alpha);
            frankcdf = @(u1,u2,alpha) ...
                -1/alpha .* log(1 + exp(-(-log((exp(-alpha .*...
                u1) - 1)./(exp(-alpha) - 1)) + -log((exp(-alpha .* u2) - 1)./(exp(-alpha) -...
                1)))) .* (exp(-alpha) - 1));
            gumbelcdf = @(u1,u2,alpha) ...
                exp(-((-log(u1)).^alpha + (-log(u2)).^alpha).^(1/alpha));
            amhcdf = @(u1,u2,alpha)...
                (1-alpha)./(exp(log((1-alpha.*(1-u1))./u1)+ log((1-alpha.*(1-u2))./u2))-alpha);
            frankcdf3 = @(u1,u2,u3,alpha)...
                -1/alpha * log(...
                1 + exp(-(...
                -log((exp(-alpha*u1)-1)./(exp(-alpha)-1))...
                + -log((exp(-alpha * u2) - 1)./(exp(-alpha)-1))...
                + -log((exp(-alpha * u3) - 1)./(exp(-alpha) - 1))...
                ))...
                *(exp(-alpha) - 1)...
                );
            claytoncdf3 = @(u1,u2,u3,alpha)...
                (1 + ...
                ( u1.^(-alpha) - 1 + u2.^(-alpha) - 1 + u3.^(-alpha) - 1)...
                )...
                .^(-1./alpha);
            gumbelcdf3 = @(u1,u2,u3,alpha)...
                exp(-((-log(u1)).^alpha + (-log(u2)).^alpha +...
                (-log(u3)).^alpha).^(1/alpha));
            
            u1=u(:,1);u2=u(:,2);
            if size(u,2)==3
                u3 = u(:,3);
                if strcmpi(cname,'clayton')
                    F = claytoncdf3(u1,u2,u3,alpha);
                elseif strcmpi(cname,'frank')
                    F = frankcdf3(u1,u2,u3,alpha);
                elseif strcmpi(cname,'gumbel')
                    F = gumbelcdf3(u1,u2,u3,alpha);
                else
                    error(['do not surpport 3D copula ' cname])
                end
            else
                if strcmpi(cname,'clayton')
                    F = claytoncdf(u1,u2,alpha);
                elseif strcmpi(cname,'frank')
                    F = frankcdf(u1,u2,alpha);
                elseif strcmpi(cname,'gumbel')
                    F = gumbelcdf(u1,u2,alpha);
                elseif strcmpi(cname,'amh')
                    F = amhcdf(u1,u2,alpha);
                else
                    error(['do not surpport 2D copula ' cname])
                end
            end
        end
        
        
    end
end
            
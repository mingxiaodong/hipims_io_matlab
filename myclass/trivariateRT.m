classdef trivariateRT
    %bivariateRT Summary of this class goes here
    %   Detailed explanation goes here    
    properties
        Xdata % observed value
        XDT   % observing date and time
        thV   % value threshold
        thP   % propotion threshold
        unit  % the unit of Xdata
        name  % the name of Xdata
        yearnum % number of years of the observation
        PureEvent % identified events exluding neigbouring records
        Mt
        pd % probability distribution
        testKS
        testChi
        DependStruct
        copulaName
        copulaParam
        copulaCDF
    end 
    methods
        function obj = trivariateRT(DataTable1,DataTable2,DataTable3,varargin)
            %univariateRT Construct an instance of this class
            % obj = univariateRT(Xtable,Ytable)          
            %   pass xdata and datetime
            matchedTable = MatchRecordsOnDate(DataTable1,DataTable2,DataTable3);
            obj.XDT = table2array(matchedTable(:,1));
            obj.Xdata = table2array(matchedTable(:,2:end));
            [I,~] = ind2sub(size(obj.Xdata),find(isnan(obj.Xdata)));
            ind = unique(I);
            obj.Xdata(ind,:)=[];
            obj.XDT(ind)=[];
            obj.yearnum = numel(obj.XDT)/365.25;
            if isempty(varargin)
                obj.thP = [0.95, 0.95, 0.95];
            elseif numel(varargin)==1
                obj.thP = varargin{1};
            else
                if strcmpi(varargin{1},'P')
                    obj.thP = varargin{2};
                elseif strcmpi(varargin{1},'V')
                    obj.thV = varargin{2};
                end
            end
            if numel(varargin)==3
                dayGap = varargin{3};
            else
                dayGap = 0;
            end
            if isempty(obj.thP)
                EP = nan(1,size(obj.Xdata,2));
                for i=1:numel(EP)
                    EP(i) = (sum(obj.Xdata(:,i)>=obj.thV(i))+1-0.44)/(numel(obj.Xdata(:,i))+0.12);
                end
                obj.thP = 1-EP;
            end
            if isempty(obj.thV)
                obj.thV = nan(1,size(obj.Xdata,2));
                for i=1:numel(obj.thV)
                    obj.thV(i) = quantile(obj.Xdata(:,i),obj.thP(i));
                end
            end
            ind = obj.Xdata(:,1)>=obj.thV(1)&...
                  obj.Xdata(:,2)>=obj.thV(2)&...
                  obj.Xdata(:,3)>=obj.thV(3);
            obj.PureEvent = ExtractDuplicateEvents(obj.XDT(ind),obj.Xdata(ind,:),dayGap);
            obj.Mt = obj.yearnum/height(obj.PureEvent);
            obj.DependStruct = TriCorr(obj);
%             obj = distFit(obj);
        end
        
        function obj = distFit(obj,varargin)
            % varargin: distribution names(one or three) 
            defaultDistName='GeneralizedExtremeValue';
            expectedDistName = {'GeneralizedExtremeValue','GeneralizedPareto'};
            
            validClassPosObj = @(x) isa(x,'trivariateRT');
            validStringParam = @(x) any(validatestring(x,expectedDistName));
            p_arg = inputParser;
            addRequired(p_arg,'object',validClassPosObj);
            addOptional(p_arg,'distname1',defaultDistName,validStringParam);
            addOptional(p_arg,'distname2',defaultDistName,validStringParam);
            addOptional(p_arg,'distname3',defaultDistName,validStringParam);
            parse(p_arg,obj,varargin{:});
            

            distname1 = p_arg.Results.distname1;
            distname2 = p_arg.Results.distname2;
            distname3 = p_arg.Results.distname3;
            distnames = {distname1,distname2,distname3};
            pdCell = cell(1,3);
            for i=1:length(distnames)
                if strcmp(distnames{i},'GeneralizedPareto')
                    pdCell{i} = fitdist(obj.PureEvent.Event(:,i),distnames{i},'Theta',obj.thV(i)-0.000001);
                else
                    pdCell{i} = fitdist(obj.PureEvent.Event(:,i),distnames{i});
                end
            end
            pd1 = pdCell{1};
            pd2 = pdCell{2};
            pd3 = pdCell{3};
            Xe_pure1 = obj.PureEvent.Event(:,1);
            Xe_pure2 = obj.PureEvent.Event(:,2);
            Xe_pure3 = obj.PureEvent.Event(:,3);
            h1 = kstest(Xe_pure1,'CDF',[Xe_pure1, cdf(pd1,Xe_pure1)]);
            h2 = kstest(Xe_pure2,'CDF',[Xe_pure2, cdf(pd2,Xe_pure2)]);
            h3 = kstest(Xe_pure3,'CDF',[Xe_pure3, cdf(pd3,Xe_pure3)]);
            obj.testKS = [h1 h2 h3];
            h1 = chi2gof(Xe_pure1,'CDF',pd1);
            h2 = chi2gof(Xe_pure2,'CDF',pd2);
            h3 = chi2gof(Xe_pure3,'CDF',pd3);
            obj.testChi = [h1 h2 h3];
            obj.pd = [pd1,pd2,pd3];
            
            disp('GOT fit test results:')
            disp([pd1.DistributionName,' for ' obj.name{1}])
            disp([pd2.DistributionName,' for ' obj.name{2}])
            disp([pd3.DistributionName,' for ' obj.name{3}])
            disp(['KS Test: h=' num2str(obj.testKS)])
            disp(['Chi-squre Test: h=' num2str(obj.testChi)])
        end
        
        function Output1 = TriCorr(obj)
            Output1 = struct('type','Spearman','rho',[],'pvalue',[]);
            Output1(2).type = 'Kendall';
            [Rho_S,P_S] = corr(obj.PureEvent.Event,'Type','Spearman');% Spearman’s rho
            Output1(1).rho = Rho_S;
            Output1(1).pvalue = P_S;
            [Rho_K,P_K] = corr(obj.PureEvent.Event,'Type','Kendall');% Kendall's Tau
            Output1(2).rho = Rho_K;
            Output1(2).pvalue = P_K;            
        end
        
        function obj = CopulaFit(obj,name,alpha)
            obj.copulaName = name;
            obj.copulaParam = alpha;
            if strcmpi(name,'clayton')
                cdffun = @(u1,u2,u3)...
                    (1 + ...
                    ( u1.^(-alpha) - 1 + u2.^(-alpha) - 1 + u3.^(-alpha) - 1)...
                    )...
                    .^(-1./alpha);
            elseif strcmpi(name,'frank')
                cdffun = @(u1,u2,u3)...
                    -1/alpha * log(...
                    1 + exp(-(...
                    -log((exp(-alpha*u1)-1)./(exp(-alpha)-1))...
                    + -log((exp(-alpha * u2) - 1)./(exp(-alpha)-1))...
                    + -log((exp(-alpha * u3) - 1)./(exp(-alpha) - 1))...
                    ))...
                    *(exp(-alpha) - 1)...
                    );
            elseif strcmpi(name,'gumbel')
                cdffun = @(u1,u2,u3)...
                    exp(-((-log(u1)).^alpha + (-log(u2)).^alpha +...
                    (-log(u3)).^alpha).^(1/alpha));
            end
            obj.copulaCDF = cdffun;
        end
        
        function Fxyz = JointCDF(obj,x1,x2,x3,varargin)
            %jointF = JointCDF(obj,x1,x2,x3,objUni1,,objUni2,,objUni3)
            if isempty(varargin)
                u1 = cdf(obj.pd(1),x1);
                u2 = cdf(obj.pd(2),x2);
                u3 = cdf(obj.pd(3),x3);
                Fxyz = obj.copulaCDF(u1,u2,u3);
            elseif numel(varargin)==3&&sum(cellfun(@isobject,varargin))==3
                obj1 = varargin{1};
                obj2 = varargin{2};
                obj3 = varargin{3};
                u1 = cdf(obj1.pd,x1);
                u2 = cdf(obj2.pd,x2);
                u3 = cdf(obj3.pd,x3);
                Fxyz = obj.copulaCDF(u1,u2,u3);
            end
        end
        
        function [RTV,obj] = JointRT(obj,x1,x2,x3,varargin)
            %RTV = JointRT(obj,x1,x2,x3,objUni1,,objUni2,,objUni3)
            if isempty(varargin)
                u1 = cdf(obj.pd(1),x1);
                u2 = cdf(obj.pd(2),x2);
                u3 = cdf(obj.pd(3),x3);
                Fxy = obj.copulaCDF(u1,u2,u3*0+1);
                Fxz = obj.copulaCDF(u1,u2*0+1,u3);
                Fyz = obj.copulaCDF(u1*0+1,u2,u3);
                jointF3 = JointCDF(obj,x1,x2,x3);
                RTV = obj.Mt./(1-u1-u2-u3+Fxy+Fxz+Fyz-jointF3);
            else
                obj1 = varargin{1};
                obj2 = varargin{2};
                obj3 = varargin{3};
                u1 = cdf(obj1.pd,x1);
                u2 = cdf(obj2.pd,x2);
                u3 = cdf(obj3.pd,x3);
%                 ind1 = x1<=obj1.thV;
%                 if sum(ind1)>0
%                     u1(ind1) = cdf(obj.pd(1),x1(ind1));
%                 end
%                 ind2 = x2<=obj2.thV;
%                 if sum(ind2)>0
%                     u2(ind2) = cdf(obj.pd(2),x2(ind2));
%                 end
%                 ind3 = x3<=obj3.thV;
%                 if sum(ind3)>0
%                     u3(ind3) = cdf(obj.pd(3),x3(ind3));
%                 end
%                 ind = ind1|ind2|ind3;
                Fxy = obj.copulaCDF(u1,u2,u3*0+1);
                Fxz = obj.copulaCDF(u1,u2*0+1,u3);
                Fyz = obj.copulaCDF(u1*0+1,u2,u3);
                jointF3 = obj.copulaCDF(u1,u2,u3);
                
                meanMt = (obj1.Mt+obj2.Mt+obj3.Mt)/3;
                RTV = meanMt./(1-u1-u2-u3+Fxy+Fxz+Fyz-jointF3);
%                 RTV0 = obj.Mt./(1-u1-u2-u3+Fxy+Fxz+Fyz-jointF3);
%                 RTV(ind) = RTV0(ind);
            end
        end
        
        function XYZ = iJointRT(obj,RT,obj1,obj2,obj3)
            xlow = obj1.thV;
            ylow = obj2.thV;
            zlow = obj3.thV;
            x1 = [xlow,xlow+range(obj1.PureEvent.Event)*10]; x2 = ylow; x3 = zlow;
            Sol1 = solveRTwithOneVar(obj,x1,x2,x3,obj1,obj2,obj3,RT);
            x1 = xlow; x2 = [ylow,ylow+range(obj2.PureEvent.Event)*10]; x3 = zlow;
            Sol2 = solveRTwithOneVar(obj,x1,x2,x3,obj1,obj2,obj3,RT);            
            x1 = xlow; x2 = ylow; x3 = [zlow,zlow+range(obj3.PureEvent.Event)*10];
            Sol3 = solveRTwithOneVar(obj,x1,x2,x3,obj1,obj2,obj3,RT);
            xsol = [Sol1(1),Sol2(1),Sol3(1)];
            ysol = [Sol1(2),Sol2(2),Sol3(2)];
            zsol = [Sol1(3),Sol2(3),Sol3(3)];
%             x1 = linspace(min(xsol),max(xsol),150);
            x1 = min(xsol):0.5:max(xsol);
%             x2 = linspace(min(ysol),max(ysol),150);
            x2 = min(ysol):5:max(ysol);
%             x3 = linspace(min(zsol),max(zsol),150);
            x3 = min(zsol):0.02:max(zsol);

%             XYV = [x1',x2',x3'];
            [XX,YY,ZZ] = meshgrid(x1,x2,x3);
            TT = JointRT(obj,XX(:),YY(:),ZZ(:),obj1,obj2,obj3);
            TT = reshape(TT,size(XX));
            Tdiffs = abs(TT-RT);
            [~,Ir] = min(Tdiffs,[],3);
            [a,b] = ind2sub(size(TT(:,:,1)),(1:numel(TT(:,:,1)))');            
            I = sub2ind(size(TT),a,b,Ir(:));
            Ir = Ir(:); I(Ir==1)=[];
%             ZZ(Tdiffs>0.5)=nan;
            XYZ = [XX(I),YY(I),ZZ(I)];
%             XYZ(isnan(XYZ(:,3)),:)=[];
        end

        function Xsol = solveRTwithOneVar(obj,x1,x2,x3,obj1,obj2,obj3,T)
            % x1, x2, x3: two are scalar(given value), the other is two 
            % elements vector (to be solved) showing the range of solve
            if numel(x2)==2
                x2 = linspace(x2(1),x2(2),range(x2)*100);
                x1 = x2*0+x1;
                x3 = x2*0+x3;
            elseif numel(x1)==2
                x1 = linspace(x1(1),x1(2),range(x1)*100);
                x2 = x1*0+x2;
                x3 = x1*0+x3;
            else
                x3 = linspace(x3(1),x3(2),range(x3)*100);
                x2 = x3*0+x2;
                x1 = x3*0+x1;
            end
            RTV = JointRT(obj,x1,x2,x3,obj1,obj2,obj3);
            Tdiffs = abs(RTV-T);
            [~,I] = min(Tdiffs);
            if numel(I)>1
                I=I(1);
            end
            Xsol = [x1(I) x2(I) x3(I)];
        end        
    %% visualization
    function Output1 = ScatterPlot(obj)
        Output1 = obj.PureEvent.Event;
        depStruct = obj.DependStruct;
        names = obj.name; varNames = names;
        units = obj.unit;
        for i=1:3
            varNames{i} = [names{i}, '(' units{i},')'];
        end
        [~,AX,BigAx,H,~] = plotmatrix(Output1);
        % axis tight
        for i=1:3
            if H(i).NumBins<6
                H(i).NumBins=6;
            end
            for j=1:3
                ax = AX(i,j);
                ax.XLim = [min(Output1(:,j))*0.9,max(Output1(:,j)*1.1)];
                ax.YLim = [min(Output1(:,i))*0.9,max(Output1(:,i))*1.1];
                ax.YTickLabelRotation=90;
                if i~=j
                    x = ax.XLim(1)+range(ax.XLim)*0.05;
                    y = ax.YLim(1)+range(ax.YLim)*0.9;
                    str1 = ['\rho_S = ' num2str(depStruct(1).rho(i,j),'%.4f')];%,...
%                         ', \itp = ',num2str(depStruct(1).pvalue(i,j),'%.2f')];
                    str2 = ['\tau_K = ' num2str(depStruct(2).rho(i,j),'%.4f')];%,...
%                         ', \itp = ',num2str(depStruct(2).pvalue(i,j),'%.2f')];
                    if i>j
                        text(ax,x,y,str1)
                    else
                        text(ax,x,y,str2)
                    end
                end
            end
        end
        text(BigAx,[0.1 0.45 0.8], repmat(0.05,1,3), varNames);
        text(BigAx,repmat(-0.06,1,3), [1.7 0.85 0]+0.5, varNames,'Rotation',90);        
    end
    
    function Output1 = QQPlot(obj)
        for i=1:3
            subplot(1,3,i)
            qqplot(obj.PureEvent.Event(:,i),obj.pd(i))
            ylabel(['Observed ' obj.name{i} ': ' obj.unit{i} ])
            xlabel(['Theoretical ' obj.name{i} ': ' obj.unit{i} ])
            title(obj.pd(i).DistributionName)
            axis image square
        end
        Output1 = obj.PureEvent.Event;
        
    end
    
    function CDFPlot(obj)
        alpha = 0.05;
        xsample = obj.PureEvent.Event;
        [~,x1] = ecdf(xsample(:,1));
        [~,x2] = ecdf(xsample(:,2));
        [~,x3] = ecdf(xsample(:,3));
        x1CDF = obj.CDFwithBounds(obj.pd(1),x1,alpha);
        x2CDF = obj.CDFwithBounds(obj.pd(2),x2,alpha);
        x3CDF = obj.CDFwithBounds(obj.pd(3),x3,alpha);
        subplot(1,3,1)
        ecdf(x1,'bounds','on');
        hold on
        plot(x1CDF(:,1),x1CDF(:,2),'r-')
        plot(x1CDF(:,1),x1CDF(:,3:4),'r--')
        hold off
        axis image square
        
        subplot(1,3,2)
        ecdf(x2,'bounds','on');
        hold on
        plot(x2CDF(:,1),x2CDF(:,2),'r-')
        plot(x2CDF(:,1),x2CDF(:,3:4),'r--')
        hold off
        axis image square
        
        subplot(1,3,3)
        ecdf(x3,'bounds','on');
        hold on
        plot(x3CDF(:,1),x3CDF(:,2),'r-')
        plot(x3CDF(:,1),x3CDF(:,3:4),'r--')
        hold off
        axis image square
    end
    
    function Output1 = RTPlot(obj,varargin)
        Output1 = obj.PureEvent.Event;
        % xsample = RTPlot(obj,xrange)        
        xsample1 = obj.PureEvent.Event(:,1);
        xsample2 = obj.PureEvent.Event(:,2);
        xsample3 = obj.PureEvent.Event(:,3);
        if numel(varargin)<3
            xrange1 = linspace(min(xsample1),min(xsample1)+range(xsample1),100);
            xrange2 = linspace(min(xsample2),min(xsample2)+range(xsample2),100);
            xrange3 = linspace(min(xsample3),min(xsample3)+range(xsample3),100);
        else
            xrange1 = varargin{1};
            xrange2 = varargin{2};
            xrange3 = varargin{3};
        end
        pd1 = obj.pd(1); pd2 = obj.pd(2); pd3 = obj.pd(3);
        [RT_e1, RT_t1] = obj.CalculateReturnPeriod(xsample1,xrange1,obj.Mt,pd1);
        [RT_e2, RT_t2] = obj.CalculateReturnPeriod(xsample2,xrange2,obj.Mt,pd2);
        [RT_e3, RT_t3] = obj.CalculateReturnPeriod(xsample3,xrange3,obj.Mt,pd3);
        LogFlag = 'linear';
        if numel(varargin)>=1&&numel(varargin)<3&&strcmpi('bound',varargin{1})
            BoundFlag = true;
            if numel(varargin) == 2
                LogFlag = varargin{2};
            end
        else
            BoundFlag = false;
        end
        
        subplot(1,3,1)
        plot(RT_t1(:,1),RT_t1(:,2),'b')
        hold on
        plot(RT_e1(:,1),RT_e1(:,2),'r*')
        if BoundFlag
            plot(RT_t1(:,1),RT_t1(:,3:4),'b--')
        end
        hold off
        axis image square
        ax = gca; ax.YScale = LogFlag;
        title(['Return Period of ' obj.name{1}])
        xlabel([obj.name{1} ': ' obj.unit{1}])
        ylabel('Return period: year')
        
        subplot(1,3,2)
        plot(RT_t2(:,1),RT_t2(:,2),'b')
        hold on
        plot(RT_e2(:,1),RT_e2(:,2),'r*')
        if BoundFlag
            plot(RT_t2(:,1),RT_t2(:,3:4),'b--')
        end
        hold off
        axis image square
        ax = gca; ax.YScale = LogFlag;
        title(['Return Period of ' obj.name{2}])
        xlabel([obj.name{2} ': ' obj.unit{2}])
        ylabel('Return period: year')
        
        subplot(1,3,3)
        plot(RT_t3(:,1),RT_t3(:,2),'b')
        hold on
        plot(RT_e3(:,1),RT_e3(:,2),'r*')
        if BoundFlag
            plot(RT_t3(:,1),RT_t3(:,3:4),'b--')
        end
        hold off
        axis image square
        ax = gca; ax.YScale = LogFlag;
        title(['Return Period of ' obj.name{3}])
        xlabel([obj.name{3} ': ' obj.unit{3}])
        ylabel('Return period: year')
    end
    
    end
    
    methods(Static)
        function output = CDFwithBounds(pd,x,alpha)
            % output = [x,Fx,Fx_lo,Fx_up];
            if size(x,1)==1
                x = x';
            end
            ci = paramci(pd,'Alpha',alpha);
            Fx = cdf(pd,x);
            if size(ci,2)==3
                Fx_lo = cdf(pd.DistributionName,x,ci(1,1),ci(1,2),ci(1,3));
                Fx_up = cdf(pd.DistributionName,x,ci(2,1),ci(2,2),ci(2,3));
            elseif size(ci,2)==2
                Fx_lo = cdf(pd.DistributionName,x,ci(1,1),ci(1,2));
                Fx_up = cdf(pd.DistributionName,x,ci(2,1),ci(2,2));
            elseif size(ci,2)==1
                Fx_lo = cdf(pd.DistributionName,x,ci(1,1));
                Fx_up = cdf(pd.DistributionName,x,ci(2,1));
            end
            output = [x,Fx,Fx_lo,Fx_up];
        end
        
        function [RT_e, RT_t] = CalculateReturnPeriod(xempirical,xtheoretical,Mt,pd)
            [f,x,flo,fup]  = ecdf(xempirical);            
            xCDF = trivariateRT.CDFwithBounds(pd,xtheoretical,0.05);
            RT_e = Mt./(1-[f,flo,fup]);
            RT_e = [x,RT_e];
            RT_t = Mt./(1-xCDF(:,2:4));
            RT_t = [xCDF(:,1),RT_t];            
        end
        
        function PureEvent = ExtractDuplicateEvents(Dates,Observations,Daygap)
            % PureEvent = ExtractDuplicateEvents remove the multiple events in the
            %   short duration defined by Daygap. The maximum value for each variable
            %   will be selected as the represent values in that duration.
            % Dates: datetime array gives the time of observations
            % Observations: martix with columns representing the variables
            % Daygap: integer give the duration of extraction
            % PureEvent: return value, a table of variable Date and Event
            [Dates,I] = sort(Dates);
            Observations = Observations(I,:);
            gaps = [Daygap+1;Dates(2:end)-Dates(1:end-1)];
            indGaps = gaps<=Daygap;
            Date = Dates;
            Event = Observations;
            for i=1:length(Dates)
                if indGaps(i)
                    Date(i) = Date(i-1);
                    Event(i,:) = max(Event(i-1:i,:));
                end
            end
            PureEvent = table(Date,Event);
            [~,ia,~] = unique(Date,'last');
            PureEvent = PureEvent(ia,:);
        end
    end
end
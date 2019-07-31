classdef bivariateRT
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
        function obj = bivariateRT(DataTable1,DataTable2,varargin)
            % bivariateRT Construct an instance of this class
            % obj = univariateRT(Xtable,Ytable)
            % obj = univariateRT(Xtable,Ytable,P,probValues)
            % obj = univariateRT(Xtable,Ytable,V,TrueValues)
            % obj = univariateRT(Xtable,Ytable,V,TrueValues,dagGap) 
            %   pass xdata and datetime
            matchedTable = MatchRecordsOnDate(DataTable1,DataTable2);
            obj.XDT = table2array(matchedTable(:,1));
            obj.Xdata = table2array(matchedTable(:,2:end));
            [I,~] = ind2sub(size(obj.Xdata),find(isnan(obj.Xdata)));
            ind = unique(I);
            obj.Xdata(ind,:)=[];
            obj.XDT(ind)=[];
            obj.yearnum = numel(obj.XDT)/365.25;
            if isempty(varargin)
                obj.thP = [0.95 0.95];
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
                dayGap = 3;
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
            ind = obj.Xdata(:,1)>=obj.thV(1)&obj.Xdata(:,2)>=obj.thV(2);
            obj.PureEvent = ExtractDuplicateEvents(obj.XDT(ind),obj.Xdata(ind,:),dayGap);
            obj.Mt = obj.yearnum/height(obj.PureEvent);
            obj.DependStruct = BiCorr(obj);
%             obj = distFit(obj);
        end
        
        function obj = distFit(obj,varargin)
            % varargin: distribution names 
            defaultDistName='GeneralizedExtremeValue';
            expectedDistName = {'GeneralizedExtremeValue','GeneralizedPareto'};
            
            validClassPosObj = @(x) isa(x,'bivariateRT');
            validStringParam = @(x) any(validatestring(x,expectedDistName));
            p_arg = inputParser;
            addRequired(p_arg,'object',validClassPosObj);
            addOptional(p_arg,'distname1',defaultDistName,validStringParam);
            addOptional(p_arg,'distname2',defaultDistName,validStringParam);
            parse(p_arg,obj,varargin{:});

            Xe_pure1 = obj.PureEvent.Event(:,1);
            Xe_pure2 = obj.PureEvent.Event(:,2);
            distname1 = p_arg.Results.distname1;
            distname2 = p_arg.Results.distname2;
            if strcmp(distname1,'GeneralizedPareto')
                pd1 = fitdist(Xe_pure1,distname1,'Theta',obj.thV(1)-0.000001);
            else
                pd1 = fitdist(Xe_pure1,distname1);
            end
            
            if strcmp(distname2,'GeneralizedPareto')
                pd2 = fitdist(Xe_pure2,distname2,'Theta',obj.thV(2)-0.000001);
            else
                pd2 = fitdist(Xe_pure2,distname2);
            end
            % GOT Test
            % h = 1 indicates the rejection of the null hypothesis that
            %   data in vector Xe_pure comes from the defined distribution 
            h1 = kstest(Xe_pure1,'CDF',[Xe_pure1, cdf(pd1,Xe_pure1)]);
            h2 = kstest(Xe_pure2,'CDF',[Xe_pure2, cdf(pd2,Xe_pure2)]);
            obj.testKS = [h1 h2];
            h1 = chi2gof(Xe_pure1,'CDF',pd1);
            h2 = chi2gof(Xe_pure2,'CDF',pd2);
            obj.testChi = [h1 h2];
            obj.pd = [pd1,pd2];
            disp('GOT fit test results:')
            disp([pd1.DistributionName,' for var1 & ', pd2.DistributionName,' for var2'])
            disp(['KS Test: h=' num2str(obj.testKS)])
            disp(['Chi-squre Test: h=' num2str(obj.testChi)])
                       
        end
        
        function Output1 = BiCorr(obj)
            Output1 = struct('type','Spearman','rho',[],'pvalue',[]);
            Output1(2).type = 'Kendall';
            Xe = obj.PureEvent.Event(:,1);
            Ye = obj.PureEvent.Event(:,2);
            [Rho_S,P_S] = corr(Xe,Ye,'Type','Spearman');% Spearman’s rho
            Output1(1).rho = Rho_S;
            Output1(1).pvalue = P_S;
            [Rho_K,P_K] = corr(Xe,Ye,'Type','Kendall');% Kendall's Tau
            Output1(2).rho = Rho_K;
            Output1(2).pvalue = P_K;            
        end
        
        function obj = CopulaFit(obj,name,alpha)
            obj.copulaName = name;
            obj.copulaParam = alpha;            
            if strcmpi(name,'clayton')
                cdffun = @(u1,u2) (1 +(u1.^(-alpha)-1+u2.^(-alpha)- 1)).^(-1./alpha);
            elseif strcmpi(name,'frank')
                cdffun = @(u1,u2) ...
                -1/alpha .* log(1 + exp(-(-log((exp(-alpha .*...
                u1) - 1)./(exp(-alpha) - 1)) + -log((exp(-alpha .* u2) - 1)./(exp(-alpha) -...
                1)))) .* (exp(-alpha) - 1));
            elseif strcmpi(name,'gumbel')
                cdffun = @(u1,u2) ...
                exp(-((-log(u1)).^alpha + (-log(u2)).^alpha).^(1/alpha));
            elseif strcmpi(name,'amh')
                cdffun = @(u1,u2)...
                (1-alpha)./(exp(log((1-alpha.*(1-u1))./u1)+ log((1-alpha.*(1-u2))./u2))-alpha);
            elseif strcmpi(name,'indep')
                cdffun = @(u1,u2) u1.*u2;
            end
            obj.copulaCDF = cdffun;
        end
        
        function jointF = JointCDF(obj,x1,x2,varargin)
            % jointF = JointCDF(objDH,x1,x2,objSH1,objSH2)
            if isempty(varargin)
                u1 = cdf(obj.pd(1),x1);
                u2 = cdf(obj.pd(2),x2);
                jointF = obj.copulaCDF(u1,u2);
            elseif numel(varargin)==2&&sum(cellfun(@isobject,varargin))==2
                obj1 = varargin{1};
                obj2 = varargin{2};
                u1 = cdf(obj1.pd,x1);
                u2 = cdf(obj2.pd,x2);
                jointF = obj.copulaCDF(u1,u2);
            end
        end
        
        function RTV = JointRT(obj,x1,x2,varargin)
            % RTV = JointRT(objDH,x1,x2,objSH1,objSH2)
            if isempty(varargin)
                u1 = cdf(obj.pd(1),x1);
                u2 = cdf(obj.pd(2),x2);
                jointF = JointCDF(obj,x1,x2);
                RTV = obj.Mt./(1-u1-u2+jointF);
            else
                obj1 = varargin{1};
                obj2 = varargin{2};
                u1 = cdf(obj1.pd,x1);
                u2 = cdf(obj2.pd,x2);
%                 ind1 = x1<=obj1.thV;
%                 if sum(ind1)>0
%                 u1(ind1) = cdf(obj.pd(1),x1(ind1));
%                 end
%                 ind2 = x2<=obj2.thV;
%                 if sum(ind2)>0
%                 u2(ind2) = cdf(obj.pd(2),x2(ind2));
%                 end
%                 ind = ind1&ind2;
                jointF = obj.copulaCDF(u1,u2);
                meanMt = (obj1.Mt+obj2.Mt)/2;
                RTV = meanMt./(1-u1-u2+jointF);
%                 RTV0 = obj.Mt./(1-u1-u2+jointF);
%                 RTV(ind) = RTV0(ind);
            end
        end
        
        function XYV = iJointRT(obj,RT,obj1,obj2)
            % RT: scalar(year)
            xlow = obj1.thV;
            ylow = obj2.thV;
            x1 = xlow; x2 = [ylow,ylow+range(obj2.PureEvent.Event)*20];
            Sol1 = solveRTwithOneVar(obj,x1,x2,obj1,obj2,RT);            
            x1 = [xlow,xlow+range(obj1.PureEvent.Event)*20]; x2 = ylow;
            Sol2 = solveRTwithOneVar(obj,x1,x2,obj1,obj2,RT);  
            xsol = [Sol1(1),Sol2(1)];
            ysol = [Sol1(2),Sol2(2)];
            x1 = linspace(min(xsol),max(xsol),1000);
            x2 = linspace(min(ysol),max(ysol),1000);
            [XX,YY] = meshgrid(x1,x2);
            TT = JointRT(obj,XX(:),YY(:),obj1,obj2);
            TT = reshape(TT,size(XX));
            [~,Ir] = min(abs(TT-RT),[],2);
%             TT0 = TT*0;
%             Ir = unique(Ir);
%             TT0(:,Ir)=nan; I = isnan(TT0);
            I = sub2ind(size(TT),(1:size(TT,1))',Ir);
            XYV = [XX(I),YY(I)];
        end
        
        function Xsol = solveRTwithOneVar(obj,x1,x2,obj1,obj2,T)
            % x1 and x2: one is scalar(given value), the other is two 
            % elements vector (to be solved) showing the range of solve
            if numel(x1)==1&&numel(x2)==2
                x2 = linspace(x2(1),x2(2),1000);
                x1 = x2*0+x1;
            else
                x1 = linspace(x1(1),x1(2),1000);
                x2 = x1*0+x2;
            end
            RTV = JointRT(obj,x1,x2,obj1,obj2);
            Tdiffs = abs(RTV-T);
            [~,I] = min(Tdiffs);
            if numel(I)>1
                I=I(1);
            end
            Xsol = [x1(I) x2(I)];
        end
        
    %% visualizationP
    function Output1 = ScatterPlot(obj)
        % plot scatter plot for the two varaibles and calculate their
        % dependence value
        Output1 = obj.PureEvent.Event;
        scatter(obj.PureEvent.Event(:,1),obj.PureEvent.Event(:,2),'filled')
        str1 = ['\rho_S = ' num2str(obj.DependStruct(1).rho,'%.4f'),...
            ', \itp = ',num2str(obj.DependStruct(1).pvalue,'%.2f')];
        str2 = ['\tau_K = ' num2str(obj.DependStruct(2).rho,'%.4f'),...
            ', \itp = ',num2str(obj.DependStruct(2).pvalue,'%.2f')];
        ax = gca;
        x = ax.XLim(2)-range(ax.XLim)*0.3;
        y = ax.YLim(1)+range(ax.YLim)*0.9;
        text(x,y,[str1; str2])
        xlabel(obj.unit{1})
        ylabel(obj.unit{2})
    end
    
    function Output1 = QQPlot(obj,fig_h)
        % plot QQplot for each variable
        % QQPlot(obj)
        % QQPlot(obj,fig_h)        
        % fig_h: a handle of a figure to plot
        % return two axes
        if nargin==2
            ax1 = axes(fig_h,'OuterPosition',[0,0,0.5,0.5]);
            ax2 = axes(fig_h,'OuterPosition',[0.5,0.5,0.5,0.5]);
            
        else
            ax1 = subplot(1,2,1);
            ax2 = subplot(1,2,2);
        end
        axCell = {ax1,ax2};
        for i=1:2
            axes(axCell{i})
            qqplot(obj.PureEvent.Event(:,i),obj.pd(i))
            ylabel(['Observed ' obj.name{i} ': ' obj.unit{i} ])
            xlabel(['Theoretical ' obj.name{i} ': ' obj.unit{i} ])
            title(obj.pd(i).DistributionName)
            axis image square
%             axCell{i}=gca;
        end
        Output1 = axCell;
    end
    
    function axCell=CDFPlot(obj,fig_h)
        % Plot Marginal CDF curves for each variable
        % axCell=CDFPlot(obj);
        % axCell=CDFPlot(obj,fig_h);
        % fig_h: a handle of a figure to plot
        % return two axes
        
        if nargin==2
            ax1 = axes(fig_h,'OuterPosition',[0,0,0.5,0.5]);
            ax2 = axes(fig_h,'OuterPosition',[0.5,0.5,0.5,0.5]);            
        else
            ax1 = subplot(1,2,1);
            ax2 = subplot(1,2,2);
        end
        axCell = {ax1,ax2};
        xsample = obj.PureEvent.Event;
        for i=1:2
            axes(axCell{i})
            [~,x1] = ecdf(xsample(:,i));
            ci = paramci(obj.pd(i));
            Fx1_CI0 = cdf(obj.pd(i).DistributionName,x1,ci(1,1),ci(1,2),ci(1,3));
            Fx1_CI1 = cdf(obj.pd(i).DistributionName,x1,ci(2,1),ci(2,2),ci(2,3));
            ecdf(xsample(:,i));%,'bounds','on'
            hold on
            plot(x1,cdf(obj.pd(i),x1),'r-')
            plot(x1,[Fx1_CI0,Fx1_CI1],'r--')
            hold off
            xlabel(['x: ' obj.name{i} '(' obj.unit{i} ')'])
            ylabel('CDF: F(x)')
            title(obj.pd(i).DistributionName)
            axis image square
%             axCell{i}=gca;
        end
        
    end
    
    function Output1 = RTPlot(obj,varargin)
        Output1 = obj.PureEvent.Event;
        % xsample = RTPlot(obj,xrange)        
        xsample1 = obj.PureEvent.Event(:,1);
        xsample2 = obj.PureEvent.Event(:,2);
        if isempty(varargin)
            xrange1 = linspace(min(xsample1),max(xsample1)*1.5,100);
            xrange2 = linspace(min(xsample2),max(xsample2)*1.5,100);
        elseif numel(varargin)>=2
            xrange1 = varargin{1};
            xrange2 = varargin{2};
        end
        pd1 = obj.pd(1); pd2 = obj.pd(2);
        [RT_e1, RT_t1] = obj.CalculateReturnPeriod(xsample1,xrange1,obj.Mt,pd1);
        [RT_e2, RT_t2] = obj.CalculateReturnPeriod(xsample2,xrange2,obj.Mt,pd2);
        
        if numel(varargin)==3&&strcmpi('bound',varargin{3})
            BoundFlag = true;
        else
            BoundFlag = false;
        end
        
        subplot(1,2,1)
        plot(RT_t1(:,1),RT_t1(:,3),'b')
        hold on
        plot(RT_e1(:,1),RT_e1(:,3),'r*')
        if BoundFlag
            plot(RT_t1(:,1),RT_t1(:,[2,4]),'b--')
        end
        hold off
        axis image square
        title(['Return Period of ' obj.name{1}])
        xlabel([obj.name{1} ': ' obj.unit{1}])
        ylabel('Return period: year')
        
        subplot(1,2,2)
        plot(RT_t2(:,1),RT_t2(:,3),'b')
        hold on
        plot(RT_e2(:,1),RT_e2(:,3),'r*')
        if BoundFlag
            plot(RT_t2(:,1),RT_t2(:,[2,4]),'b--')
        end
        hold off
        axis image square
        title(['Return Period of ' obj.name{2}])
        xlabel([obj.name{2} ': ' obj.unit{2}])
        ylabel('Return period: year')
    end
    
    function JointDistPlot(obj,zvalue,varargin)
        % JointDistPlot(obj,'ReturnPeriod','Surface')
        % JointDistPlot(obj,'ReturnPeriod','Contour')
        % JointDistPlot(obj,'CDF','Contour')
        % JointDistPlot(obj,'CDF','Contour',0.1:0.2:0.9)
        % zvalue : CDF, ReturnPeriod
        
        x = obj.PureEvent.Event(:,1);
        y = obj.PureEvent.Event(:,2);
        xlow = min(x);
        ylow = min(y);
        showObv = false;
%         [~,~,jointFemp] = obj.copulaEcdf([x,y]);
        rangeRatio = [1.5 1.5];
        plottype = 'Contour';
        Levels = [20 50 100 200 500];
        if ~isempty(varargin)
            ind = cellfun(@isnumeric,varargin);
            if sum(ind)>=1
            numericVar = varargin(ind);
            for i=1:numel(numericVar)
                if numel(numericVar{i})<=2
                    rangeRatio = numericVar{i};
                else
                    Levels = numericVar{i};
                end
            end
            end
            ind = cellfun(@ischar,varargin);
            params = varargin(ind);
            if ismember('Surface',params)
                plottype = 'Surface';
            else
                plottype = 'Contour';
            end
            ind = cellfun(@isobject,varargin);
            objs = varargin(ind);% the univariateRT objects
            if ~isempty(objs)
                xlow = objs{1}.thV;
                ylow = objs{2}.thV;
                showObv = false;
            end
        end
        
        
        XX = linspace(xlow,xlow+range(x)*rangeRatio(1),100);
        YY = linspace(ylow,ylow+range(y)*rangeRatio(2),100);
        [XX,YY] = meshgrid(XX,YY);         
        if strcmp(zvalue,'CDF')
            if max(Levels)>=1
                Levels = 0.1:0.2:0.9;
            end
            ZZ = JointCDF(obj,XX(:),YY(:));
            zname = 'CDF';
        elseif strcmp(zvalue,'ReturnPeriod')
            if max(Levels)<1
                Levels = [10 20 50 100 200 500];
            end
            ZZ = JointRT(obj,XX(:),YY(:));
            zname = 'Return Period (year)';
        end
        
        if strcmp(zvalue,'CDF')&&numel(objs)==2
            ZZ = JointCDF(obj,XX(:),YY(:),objs{1},objs{2});
        elseif strcmp(zvalue,'ReturnPeriod')&&numel(objs)==2
            ZZ = JointRT(obj,XX(:),YY(:),objs{1},objs{2});
        end
        ZZ = reshape(ZZ,size(XX));                
        if strcmp(plottype,'Contour')
            contour(XX,YY,ZZ,Levels,'Show','on')
            if showObv
                hold on
                scatter(x,y,'r*')
                hold off
            end
            axis image square
        elseif strcmp(plottype,'Surface')
            surf(XX,YY,ZZ)
            zlabel(zname)
            box on
            axis image square
            view(-42,9)
        end        
        xlabel([ obj.name{1} ' (' obj.unit{1} ')'])
        ylabel([ obj.name{2} ' (' obj.unit{2} ')'])
        
%         Output1 = [jointFemp jointFe];        
    end
    
    function JointRTCompare(obj)
        X = obj.PureEvent.Event;
        [x1,x2,Cn] = obj.copulaEcdf(X,'2D');
        surfc(x1,x2,Cn)
    end
    
    function fig_h=MarginalFitCompare(obj)
        % plot CDF and QQplot for each variable in a figure
        fig_h = figure('Name','Marginal Distribution');
%         pos2extent = @(x) [x(1:2),x(1)+x(3),x(2)+x(4)];
        extent2pos = @(x) [x(1:2),x(3)-x(1),x(4)-x(2)];
        axCell1 = CDFPlot(obj,fig_h);
        axCell2 = QQPlot(obj,fig_h);
        
        extent = [0,0.5,0.5,1];
        axCell1{1}.OuterPosition=extent2pos(extent);
        extent = [0.5,0.5,1,1];
        axCell1{2}.OuterPosition=extent2pos(extent);
        extent = [0,0,0.5,0.5];
        axCell2{1}.OuterPosition=extent2pos(extent);
        extent = [0.5,0,1,0.5];
        axCell2{2}.OuterPosition=extent2pos(extent);
        fig_h.InnerPosition
        fig_h.InnerPosition(3)=fig_h.InnerPosition(4);
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
            xCDF = bivariateRT.CDFwithBounds(pd,xtheoretical,0.05);
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
        function [output1,output2,Cn] = copulaEcdf(X,varargin)
            % calculate empirical copula of multiple variables
            % X is a matrix whose column number represents the number of varibales
            % U is a matrix representing cumulative probabilities for each variable
            % Cn is the empirical cdf corresponding to p
            n = size(X,1); % number of observations for each variable
            d = size(X,2); % number of variables (dimensions)
            
            funecdf = @(x) sum(X(:,1)<=x)/(n+1);
            x1 = sort(X(:,1));
            f1 = arrayfun(funecdf,x1);
            
            funecdf = @(x) sum(X(:,2)<=x)/(n+1);
            x2 = sort(X(:,2));
            f2 = arrayfun(funecdf,x2);

            if ~isempty(varargin)&&strcmpi(varargin{1},'2D')
                [x1,x2] = meshgrid(x1,x2);
                [f1,f2] = meshgrid(f1,f2);
            end
            U = [f1(:),f2(:)];
            [~,R] = sort(X);
            U_hat = R./(n+1);
            indfun = @(A) sum( sum(U_hat<=repmat(A,[n,1]), 2) == d )/n; % A is the input para of indfun
            U_table = table(U);
            Cn_table = rowfun(indfun,U_table);
            Cn = table2array(Cn_table);
            Cn = reshape(Cn,size(x1));
            output1 = x1;
            output2 = x2;
            if ~isempty(varargin)
                if strcmpi(varargin{end},'F')                
                    output1 = f1;
                    output2 = f2;
                end
            end
        end
    end
end


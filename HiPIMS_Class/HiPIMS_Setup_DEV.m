classdef HiPIMS_Setup
    properties
        Z % DEM matrix
        R % DEM spatial reference
        CaseFolder % dir of input folder
        BoundStruct %structure: Position([X Y]), Type('rigid','open'), T_h (time, water depth), T_huv (time, velocity)
        h0 = 0; % initial water depth
        uv0 = 0; % initial water velocity = [u0, v0]
        RainMask = 0;
        RainSources ={[0,0;1,0]} % cell, consisting of time series of rainfall rate (m/s)
        manning = 0.035;
        sewer_sink = 0;
        hydraulic_conductivity = 0;
        cumulative_depth = 0;
        capillary_head = 0;
        water_content_diff = 0;
        LanduseMask
        GaugeCoor = [];
        DefenceMonitor %3 cols: X,Y,Water Elevation(failing threshold)
        DefenceFailureZ
        readme
    end
    
    properties (Access = private)
        ID_ValidCell % the ID of valid cell starting from 0 at lower left corner and increasing by 1 rightwards and upwards
        ID_BoundCell % the ID of boundary cells, 1 represents the outline bound
        BoundTypeCode % BoundTypeCode.h; BoundTypeCode.hU; code to define the bound type
        AllInputFileNames = {'all',...
            'DEM',... % DEM.txt at mesh folder
            'h',...%initial water depth
            'hU',...  %initial water velocity
            'precipitation',... %initial precipitation
            'z',...   %bed elevation
            'manning',... %coefficient of friction
            'sewer_sink',...
            'cumulative_depth',...
            'hydraulic_conductivity',...
            'capillary_head',...
            'water_content_diff',...
            'precipitation_mask',...
            'precipitation_source_all',...
            'gauges_pos',... %position of gauges
            'h_BC',... %Boundary condition water depth
            'hU_BC',... %Boundary condition flow per unit width
            'eta' %initial water elevation
            };
        FieldFolder
        MeshFolder
        History
    end
    
    methods
        %% initialization
        function obj = HiPIMS_Setup(Z_dem,R_dem,varargin)
            % obj = HiPIMS_Setup(Z_dem,R_dem)
            % obj = HiPIMS_Setup(Z_dem,R_dem,CaseFolder)
            obj.Z = Z_dem;
            obj.R = R_dem;
            if length(varargin)==1
                obj.CaseFolder = varargin{1};
            elseif length(varargin)>1
                error('The max number of arguments is 3: Z,R,CaseFolder')
            else
                obj.CaseFolder = pwd;
            end
            obj.FieldFolder = [obj.CaseFolder '/input/field/'];
            obj.MeshFolder = [obj.CaseFolder '/input/mesh/'];
            obj.ID_ValidCell = GetValidID(obj);
            % default bound type: single outline bound, open with no water outside
            obj.BoundStruct = struct('Type','open','Position','outline',...
                'T_h',[0 0; 1 0],'T_huv',[0 0 0; 1 0 0]);
            obj.ID_BoundCell= GetBoundID(obj);
            obj.History = struct('DateTime',datestr(datetime('now')),'Log','Created');
            %struct('DateTime',datetime('now'),'Log','Created');
        end
        
        function obj = setBoundary(obj,boundStruct)
            % obj = setBoundary(obj,boundStruct)
            % boundStruct: Position([X Y]), Type('rigid','open'), T_h (time, water depth), T_huv (time, velocity)
            obj.BoundStruct = boundStruct;
            obj.ID_BoundCell = GetBoundID(obj);
            
            for i=1:length(obj.BoundStruct)
                
                % decompose flow into velocity x and y
                T_huv = obj.BoundStruct(i).T_huv;
                if size(T_huv,2)==2 % uv is given as flow m^3/s
                    flow = T_huv(:,2);
                    BoundWidth = abs(obj.R(2))*length(obj.ID_BoundCell(i).BoundIndAtZ);
                    if BoundWidth > 0
                        Q_per_m = flow/BoundWidth; % m^2/s
                        [I,J] = ind2sub(size(obj.Z),obj.ID_BoundCell(i).BoundIndAtZ); %range(cols):delta_x, range(a):delta_y
                        theta = atan(range(J)/range(I)); %theta: the angle between bound line and the vertical line
                        hu = Q_per_m*cos(theta); 
                        hv = Q_per_m*sin(theta); %m^2/s
                    else
                        hu = flow*0; hv = flow*0; % no inflow bound, width is zero
                    end
                    obj.BoundStruct(i).T_huv = [T_huv(:,1),hu,hv];
                end
            end
            obj.History = writeHistory(obj,'Boundary changed');
        end
        
        function obj = setRainfall(obj,rainMask,rainSources)
            obj.RainMask = rainMask;
            obj.RainSources = rainSources;
            if ~iscell(rainSources)
                rainRate = rainSources(:,2:end);
            else
                rainRate = rainSources{1}(:,2:end);
            end
            if max(rainRate(:))>1.0000e-03
                warning('The unit of rainfall might be wrong')
            end
            obj.History = writeHistory(obj,'Rainfall changed');
        end
        
        function obj = Rotation(obj)
            row0 = size(obj.Z,1);
            col0 = size(obj.Z,2);
            [x0,y0] = pix2map(obj.R,row0,col0);
            obj.Z = fliplr(obj.Z');
            obj.R = makerefmat(0,col0*obj.R(2),obj.R(2),-obj.R(2));
            % rotate GaugeCoor
            x = obj.GaugeCoor(:,1);y = obj.GaugeCoor(:,2);
            obj.GaugeCoor = [y-y0,-(x-x0)];
            % rotate BoundStruct
            newBoundStruct = obj.BoundStruct;
            for i=1:length(newBoundStruct)
                if isnumeric(newBoundStruct(i).Position)
                    newPosition = newBoundStruct(i).Position;
                    newPosition = [newPosition(:,2)-y0,x0-newPosition(:,1)];
                    newBoundStruct(i).Position = newPosition;
                end
                if ~isempty(newBoundStruct(i).T_huv)
                    T_huv = newBoundStruct(i).T_huv;
                    T_huv = [T_huv(:,1),T_huv(:,3),T_huv(:,2)];
                    newBoundStruct(i).T_huv = T_huv;
                end
            end
            obj.BoundStruct = newBoundStruct;
            obj.ID_ValidCell = GetValidID(obj);
            obj.ID_BoundCell= GetBoundID(obj);
            % rotate defence failure variables
            if ~isempty(obj.DefenceMonitor)
                x = obj.DefenceMonitor(:,1);y = obj.DefenceMonitor(:,2);
                obj.DefenceMonitor = [y-y0,-(x-x0),obj.DefenceMonitor(:,3)];
            end
            newDefenceFailureZ = obj.DefenceFailureZ;
            if ~isempty((newDefenceFailureZ))
                for i=1:length(newDefenceFailureZ)
                    ind_pre = newDefenceFailureZ{i}(:,1);
                    newDEMValue = newDefenceFailureZ{i}(:,2);
                    ind_rot = HiPIMS_Setup.rotateInd(size(obj.Z,2),size(obj.Z,1),ind_pre);
                    newDefenceFailureZ{i} = [ind_rot newDEMValue];
                end
            end
            obj.DefenceFailureZ = newDefenceFailureZ;
            % rotate Z type variables
            obj.manning = fliplr(obj.manning');
            obj.sewer_sink = fliplr(obj.sewer_sink');
            obj.hydraulic_conductivity = fliplr(obj.hydraulic_conductivity');
            obj.cumulative_depth = fliplr(obj.cumulative_depth');
            obj.capillary_head = fliplr(obj.capillary_head');
            obj.water_content_diff = fliplr(obj.water_content_diff');
            obj.h0 = fliplr(obj.h0');
            if numel(obj.uv0)>2
                obj.uv0 = [fliplr(obj.uv0(:,1:end/2)') fliplr(obj.uv0(:,end/2+1:end)')];
            end
            %
            obj.History = writeHistory(obj,'DEM Roatated');
        end
        
        function obj = setDefenceFailureZ(obj,lineShp,changeValue,ChangeMethod)
            [~, changeIndV] = AmendDEM(obj.Z,obj.R,lineShp,changeValue,ChangeMethod);
            obj.DefenceFailureZ = changeIndV;
            obj.History = writeHistory(obj,'Defence Failure DEM set');
        end
        
        %% write files
        function output = writeInputFile(obj,varargin)
            if isempty(varargin)
                filenames = {'all'};
            else
                filenames = varargin; % file name without extension
            end
            [~, Locb] = ismember(filenames,obj.AllInputFileNames);
            Locb(Locb==0)=[];
            if isempty(Locb)
                disp(filenames)
                error('Files name not in the list')
            end
            WriteFlag = Locb; % flag of files to be written
            creatFolders(obj)
            obj.FieldFolder = [obj.CaseFolder '/input/field/'];
            obj.MeshFolder = [obj.CaseFolder '/input/mesh/'];
            obj.ID_ValidCell = GetValidID(obj);
            if WriteFlag(1)==1 % all
                WriteFlag = 1:17;
            end
            output = obj.AllInputFileNames(WriteFlag); output = output';
            for i = WriteFlag
                switch i
                    case 1 % write all input files
                        disp('Writing all input flies...')
                    case 2 % DEM
                        writeName = [obj.MeshFolder 'DEM.txt'];
                        HiPIMS_Setup.Arcgridwrite(writeName,obj.Z,obj.R)
                        disp(writeName)
                    case 3 % initial water depth
                        writeName = [obj.FieldFolder 'h.dat'];
                        zValue = obj.h0+zeros(size(obj.Z));
                        writeFile_zType(obj,writeName,zValue,'h')
                        disp(writeName)
                    case 4 % initial water velocity
                        writeName = [obj.FieldFolder 'hU.dat'];
                        if isscalar(obj.uv0)
                            zValue = obj.uv0+zeros(size([obj.Z obj.Z]));
                        else
                            zValue = obj.uv0;
                        end
                        writeFile_zType(obj,writeName,zValue,'hU');
                        disp(writeName)
                    case 5 % initial precipitation
                        writeName = [obj.FieldFolder 'precipitation.dat'];
                        zValue = 0+zeros(size(obj.Z));
                        writeFile_zType(obj,writeName,zValue,[])
                        disp(writeName)
                    case 6 % bed elevation
                        writeName = [obj.FieldFolder 'z.dat'];
                        zValue = obj.Z;
                        writeFile_zType(obj,writeName,zValue,[])
                        disp(writeName)
                    case 7 % manning
                        writeName = [obj.FieldFolder 'manning.dat'];
                        zValue = obj.manning+zeros(size(obj.Z));
                        writeFile_zType(obj,writeName,zValue,[])
                        disp(writeName)
                    case 8 % sewer_sink
                        writeName = [obj.FieldFolder 'sewer_sink.dat'];
                        zValue = obj.sewer_sink+zeros(size(obj.Z));
                        writeFile_zType(obj,writeName,zValue,[])
                        disp(writeName)
                    case 9 % cumulative_depth
                        writeName = [obj.FieldFolder 'cumulative_depth.dat'];
                        zValue = obj.cumulative_depth+zeros(size(obj.Z));
                        writeFile_zType(obj,writeName,zValue,[])
                        disp(writeName)
                    case 10 % hydraulic_conductivity
                        writeName = [obj.FieldFolder 'hydraulic_conductivity.dat'];
                        zValue = obj.hydraulic_conductivity+zeros(size(obj.Z));
                        writeFile_zType(obj,writeName,zValue,[])
                        disp(writeName)
                    case 11 % capillary_head
                        writeName = [obj.FieldFolder 'capillary_head.dat'];
                        zValue = obj.capillary_head+zeros(size(obj.Z));
                        writeFile_zType(obj,writeName,zValue,[])
                        disp(writeName)
                    case 12 % water_content_diff
                        writeName = [obj.FieldFolder 'water_content_diff.dat'];
                        zValue = obj.water_content_diff+zeros(size(obj.Z));
                        writeFile_zType(obj,writeName,zValue,[])
                        disp(writeName)
                    case 13 % precipitation_mask
                        writeName = [obj.FieldFolder 'precipitation_mask.dat'];
                        rainMask = obj.RainMask+zeros(size(obj.Z));
                        writePrecipitationMask(obj,writeName,rainMask)
                    case 14 % precipitation_source
                        writePrecipitationSource(obj)
                    case 15 % gauges position
                        fileID = fopen([obj.FieldFolder 'gauges_pos.dat'],'w');
                        if isempty(obj.GaugeCoor)
                            obj.GaugeCoor = [0 0];
                        end
                        fprintf(fileID,'%.3f %.3f\n',obj.GaugeCoor(1:end-1,:)');
                        fprintf(fileID,'%.3f %.3f',obj.GaugeCoor(end,:));
                        fclose(fileID);
                        disp([obj.FieldFolder 'gauges_pos.dat'])
                    case 16 % h_BC
                        n = 0;
                        for j=1:length(obj.BoundStruct)
                            if ~isempty(obj.BoundStruct(j).T_h)
                                dlmwrite([obj.FieldFolder 'h_BC_' num2str(n) '.dat'],...
                                    obj.BoundStruct(j).T_h,'precision',10,'delimiter',' ')
                                n=n+1;
                                disp([obj.FieldFolder 'h_BC_' num2str(n) '.dat for Bound' num2str(j)])
                            end
                        end
                    case 17 % hU_BC
                        n = 0;
                        for j=1:length(obj.BoundStruct)
                            if ~isempty(obj.BoundStruct(j).T_huv)
                                dlmwrite([obj.FieldFolder 'hU_BC_' num2str(n) '.dat'],...
                                    obj.BoundStruct(j).T_huv,'precision',10,'delimiter',' ')
                                n=n+1;
                                disp([obj.FieldFolder 'hU_BC_' num2str(n) '.dat for Bound' num2str(j)])
                            end
                        end
                    case 18 % initial water elevation
                        writeName = [obj.FieldFolder 'eta.dat'];
                        zValue = obj.h0+obj.Z; zValue(isnan(zValue))=0;
                        writeFile_zType(obj,writeName,zValue,'h')
                end
            end
        end
        
        function obj = writeFailureMonitor(obj)
            mat3c = obj.DefenceMonitor;
            if ~isempty(mat3c)
                fileID = fopen([obj.FieldFolder 'failure_monitor.dat'],'w');
                fprintf(fileID,'%.3f %.3f %.3f\n',mat3c(1:end-1,:)');
                fprintf(fileID,'%.3f %.3f %.3f',mat3c(end,:));
                fclose(fileID);
                disp([obj.FieldFolder 'failure_monitor.dat'])
            end
        end
        
        function obj = writeDefenceFailureZ(obj)
            defenceFailureZ = obj.DefenceFailureZ;
            if ~isempty(defenceFailureZ)
                for i=1:length(defenceFailureZ)
                    z_changed = obj.Z; ind = defenceFailureZ{i}(:,1);
                    z_changed(ind) = defenceFailureZ{i}(:,2);
                    writeName = [obj.FieldFolder 'z_' num2str(i) '.dat'];
                    zValue = z_changed;
                    writeFile_zType(obj,writeName,zValue,[])
                    disp(writeName)
                end
            end
        end
        
        %% visualization
        function Logs=GetHistory(obj)
            Logs = obj.History;
        end
        
        function GeneralMap(obj)
            % GeneralMap(obj)
            % plot DEM and boundary cells
            [~,lgd] = PlotBoundID(obj);
            hold on
            if ~isempty(obj.GaugeCoor)
                mapshow(obj.GaugeCoor(:,1),obj.GaugeCoor(:,2),'DisplayType','Point')
                lgd.String{end} = 'Gauge Position';
            end
            mapshow(obj.Z,obj.R,'DisplayType','surface')            
            zdatam('allline',max(obj.Z(:)))
            hold off
            lgd.String(end)=[];
        end        
        
        function [BoundID,lgd] = PlotBoundID(obj)
            % plot DEM cells
            BoundID = obj.ID_BoundCell;
            hold on
            colorStr = {'k','r','b','m','y','c','g'};
            for i=1:length(BoundID)
                [row,col] = ind2sub(size(obj.Z),BoundID(i).BoundIndAtZ);
                [x,y] = pix2map(obj.R,row,col);
                n = i;
                while n>length(colorStr)
                    n=n-length(colorStr)+1;
                end
                if ~isempty(x)
                   plot(x,y,'s','MarkerEdgeColor',colorStr{n},'MarkerFaceColor',colorStr{n})
                end
            end
            hold off
            axis image; grid on; box on
            lgd = legend;
            lgd.String{1} = 'Outline Boundary';
            if length(lgd.String)>1
                for i=2:length(lgd.String)
                    lgd.String{i} = ['IO Boundary ', num2str(i-1)];
                end
            end
            lgd.Box = 'off'; %lgd.Location = 'best';
            lgd.Color = 'none';
        end
        
        function ValidID = PlotValidID(obj)
            mapshow(obj.ID_ValidCell,obj.R,'DisplayType','surface')
            axis image;
            ValidID = obj.ID_ValidCell;
            %             title('Valid cells ID')
        end
    end
    methods (Access = private)
        function NewHistory = writeHistory(obj,str)
            NewHistory = obj.History;
            n = length(NewHistory);
            NewHistory(n+1).DateTime = datestr(datetime('now'));
            NewHistory(n+1).Log = str;
        end
        
        function ValidID = GetValidID(obj)
            % get the valid ID of DEM girds starting from 0,1,2,...
            Z_changed = flipud(obj.Z); Z_changed = Z_changed';
            validInd = find(~isnan(Z_changed));
            ValidID = nan(size(Z_changed));
            ValidID(validInd) = 0:numel(validInd)-1;
            ValidID = ValidID';ValidID = flipud(ValidID);
        end
        
        function [boundID, OutlineBound]= GetBoundID(obj)
            % get the cell ID of DEM girds
            % ID starts from 0,1,2,...
            % obj must have a struct named BoundStruct with fields: BoundPosition
            z_dem = obj.Z; r_dem = obj.R;
            boundStruct = obj.BoundStruct;
            boundID = struct('BoundIndAtZ',[],'Type',[],'CodeH',[3 0 0],'CodeHU',[3 0 0]);
            Z_u = [z_dem(2:end,:);z_dem(end,:)+nan];
            Z_d = [z_dem(1,:)+nan;z_dem(1:end-1,:)];
            Z_l = [z_dem(:,2:end),z_dem(:,end)+nan];
            Z_r = [z_dem(:,1)+nan,z_dem(:,1:end-1)];
            Z_NeigbourPlus = Z_u+Z_d+Z_l+Z_r+z_dem;
            OutlineBound = isnan(Z_NeigbourPlus)&~isnan(z_dem);
            % the index of outline cells at Z matrix
            OutlineBound = find(OutlineBound); 
            
            [row,col] = ind2sub(size(z_dem),OutlineBound);
            [xq,yq] = pix2map(r_dem,row,col); % XY coords of outline points
            numBound = length(boundStruct);
            if numBound==0 %outline boundary, boundStruct is empty
                boundID.BoundIndAtZ = OutlineBound;
                boundID.Type = 'open'; %outlint boundary, boundStruct is 1
            elseif numBound==1 %only outline boundary
                boundID.BoundIndAtZ = OutlineBound;
                if isempty(boundStruct(1).Type)
                    boundID.Type = 'open';
                else
                    boundID.Type = boundStruct.Type;
                end
            else % there are more boundaries apart from outline, starting from 2
                OutlineBoundRemain = OutlineBound;
                for i=2:numBound
                    %IO boundary
                    X = boundStruct(i).Position(:,1);
                    Y = boundStruct(i).Position(:,2);
                    if sum(~isnan(X)) == 2 %if only two points for bound region polygon
                        X = [min(X);max(X);max(X);min(X)];
                        Y = [min(Y);min(Y);max(Y);max(Y)];
                        boundStruct(i).Position = [X Y];
                    end
                    ind = inpolygon(xq,yq,X,Y);
                    boundID(i).BoundIndAtZ = OutlineBound(ind);
                    boundID(i).Type = boundStruct(i).Type;
                    OutlineBoundRemain(ind) = nan;
                end
                OutlineBoundRemain(isnan(OutlineBoundRemain)) = [];
                boundID(1).BoundIndAtZ = OutlineBoundRemain;
                if isempty(boundStruct(1).Type)
                    boundID(1).Type = 'open';
                else
                    boundID(1).Type = boundStruct(1).Type;
                end
            end
            % define bound type code
            h_cnt = 0; hU_cnt = 0;
            for i=1:length(boundStruct)
                if strcmp(boundStruct(i).Type,'rigid')
                    boundID(i).CodeH  = [2 0 0];
                    boundID(i).CodeHU = [2 2 0];
                else
                    if isempty(boundStruct(i).T_h)
                        boundID(i).CodeH  = [2 0 0];
                    else
                        boundID(i).CodeH  = [3 0 h_cnt]; h_cnt = h_cnt+1;
                    end
                    if isempty(boundStruct(i).T_huv)
                        boundID(i).CodeHU  = [2 1 0];
                    else
                        boundID(i).CodeHU  = [3 0 hU_cnt]; hU_cnt = hU_cnt+1;
                    end
                end
            end
            
        end
        
        function creatFolders(obj)
            caseFolder = obj.CaseFolder;
            if ~exist([caseFolder '/input'],'dir')
                mkdir([caseFolder '/input']);
                mkdir([caseFolder '/input'],'field');
                mkdir([caseFolder '/input'],'mesh')
            elseif ~exist([caseFolder '/input/field'],'dir')
                mkdir([caseFolder '/input'],'field');
            elseif ~exist([caseFolder '/input/mesh'],'dir')
                mkdir([caseFolder '/input'],'mesh')
            end
            if ~exist([caseFolder '/output'],'dir')
                mkdir([caseFolder '/output']);
            end
        end
        
        function writePrecipitationSource(obj)
            rainSourceValue = obj.RainSources;
            if iscell(rainSourceValue)||size(rainSourceValue,2)==2
                numOfRainSource = length(rainSourceValue);
                for i = 1:numOfRainSource
                    Time_Rain = rainSourceValue{i};
                    Time_Rain(isnan(Time_Rain(:,2)),:)=[];
                    filename = [obj.FieldFolder 'precipitation_source_' num2str(i-1) '.dat'];
                    fileID = fopen(filename,'w');
                    fprintf(fileID,'%09.2f %.15e\n',Time_Rain');
                    fclose(fileID);
                end
                disp([obj.FieldFolder 'precipitation_source_0-' num2str(i)])
            else
                numOfRainSource = size(rainSourceValue,2)-1;
                formatSpec = repmat('%.15e ',[1 numOfRainSource] );
                formatSpec = ['%09.2f ' formatSpec(1:end-1) '\n'];
                fileID = fopen([obj.FieldFolder 'precipitation_source_all.dat'],'w');
                fprintf(fileID,'%u\n',numOfRainSource);
                fprintf(fileID,formatSpec,rainSourceValue');
                fclose(fileID);
                disp([obj.FieldFolder 'precipitation_source_all.dat'])
            end
        end
        
        function writePrecipitationMask(obj,filename,rainMask)
            % precipitation_mask.data
            id = obj.ID_ValidCell; id = id(:);
            id_RM = [id rainMask(:)]; % zeros(numel(Valid_ID),1) rainfall mask value is 0
            id_RM(isnan(id),:)=[]; % delete nodata value
            id_RM = sortrows(id_RM,1);
            fileID = fopen(filename,'w');
            fprintf(fileID,'$Element Number\n'); fprintf(fileID,'%d\n',length(id_RM));
            fprintf(fileID,'$Element_id  Value\n'); fprintf(fileID,'%-12d %d\n',id_RM');
            fclose(fileID);
            disp(filename)
        end
        %%
        function writeFile_zType(obj,filename,zValue,TypeFlag)
            % this type of input file has boundary type code in the end of
            % the file
            % TypeFlag = 'h','hU',[]?others?
            id_BoundCell = obj.ID_BoundCell; % struct
            id_bC = {id_BoundCell.BoundIndAtZ}'; % id cell2mat(
            % id & bound code
            for i=1:length(id_BoundCell)
                if strcmp(TypeFlag,'h')
                    if i>1 % add an value of initial water depth for IO boundary
                        zValue(id_BoundCell(i).BoundIndAtZ) = 0.0001;
                    end
                    repcodes = id_BoundCell(i).CodeH;
                elseif strcmp(TypeFlag,'hU')
                    if i>1
                        zValue(id_BoundCell(i).BoundIndAtZ) = 0.0001;
                    end
                    repcodes = id_BoundCell(i).CodeHU;
                else
                    repcodes = [2 0 0];
                end
                codes = repmat(repcodes,length(id_BoundCell(i).BoundIndAtZ),1);
                id_b0 = id_BoundCell(i).BoundIndAtZ;
                id_b0 = obj.ID_ValidCell(id_b0);
                id_bC{i} = [id_b0 codes];
            end
            id_bC = cell2mat(id_bC);
            id_bC = sortrows(id_bC);
            %           id_bC = createID_bCode_Matrix(obj,TypeFlag);
            id = obj.ID_ValidCell; id = id(:);
            if numel(zValue)==numel(obj.Z)
                id_zV = [id,zValue(:)];
            elseif numel(zValue)==2*numel(obj.Z)
                id_zV = [id,...
                    zValue(1:numel(obj.Z))',...
                    zValue(1+numel(obj.Z):end)'];
            end
            id_zV(isnan(id),:)=[];
            id_zV = sortrows(id_zV);
            %           id_zV = createID_zValue_Matrix(obj,zValue);
            if size(id_zV,2)==3
                formatSpec1 = '%-12d %.8f %.8f\n';
            else
                formatSpec1 = '%-12d %.8f\n';
            end
            fileID = fopen(filename,'w');
            % print valid cell and their initial value
            fprintf(fileID,'$Element Number\n');
            fprintf(fileID,'%d\n',length(id_zV));
            fprintf(fileID,'$Element_id  Value\n');
            fprintf(fileID,formatSpec1,id_zV');
            % print boundary cell and their code
            fprintf(fileID,'$Boundary Numbers\n');
            fprintf(fileID,'%d\n',length(id_bC));
            fprintf(fileID,'$Element_id Boundary_type\n');
            fprintf(fileID,'%-12d %d %d %d\n',id_bC');
            fclose(fileID);
        end
        %%
        function id_zV = createID_zValue_Matrix(obj,zValue)
            id = obj.ID_ValidCell; id = id(:);
            if numel(zValue)==numel(obj.Z)
                id_zV = [id,zValue(:)];
            elseif numel(zValue)==2*numel(obj.Z)
                id_zV = [id,...
                    zValue(1:numel(obj.Z))',...
                    zValue(1+numel(obj.Z):end)'];
            end
            id_zV(isnan(id),:)=[];
            id_zV = sortrows(id_zV);
        end
        %%
        function id_bC = createID_bCode_Matrix(obj,TypeFlag)
            % TypeFlag = 'h','hU',[]?others?
            id_BoundCell = obj.ID_BoundCell; % struct
            id_bC = {id_BoundCell.BoundIndAtZ}'; % id cell2mat(
            % id & bound code
            for i=1:length(id_BoundCell)
                if strcmp(TypeFlag,'h')
                    repcodes = id_BoundCell(i).CodeH;
                elseif strcmp(TypeFlag,'hU')
                    repcodes = id_BoundCell(i).CodeHU;
                else
                    repcodes = [2 0 0];
                end
                codes = repmat(repcodes,length(id_BoundCell(i).BoundIndAtZ),1);
                id_b0 = id_BoundCell(i).BoundIndAtZ; id_b0 = obj.ID_ValidCell(id_b0);
                id_bC{i} = [id_b0 codes];
            end
            id_bC = cell2mat(id_bC);
            id_bC = sortrows(id_bC);
        end        
    end
    methods (Static)
        function ind_rot = rotateInd(nrow_pre,ncol_pre,ind_pre)
            [row,col] = ind2sub([nrow_pre,ncol_pre],ind_pre);
            row_rot = col;
            col_rot = nrow_pre-row+1;
            ind_rot = sub2ind([ncol_pre,nrow_pre],row_rot,col_rot);
        end
        %%
        function Arcgridwrite(filename,Z,R)
            %ARCGRIDWRITE Write gridded data set to Arc ASCII Grid Format
            x00 = R(3,1);
            y00 = R(3,2);
            gridSize = abs(R(2));
            nrows = size(Z,1);
            ncols = size(Z,2);
            
            % XLLCORNER and YLLCORNER are the coordinates of the
            % lower left corner of the lower left cell of Z.
            xllcorner = x00 + gridSize/2;
            yllcorner = y00 - gridSize*nrows - gridSize/2;
            NODATA_value = -9999;
            Z(isnan(Z)) = NODATA_value;
            
            fileID = fopen(filename,'w');
            fprintf(fileID,'ncols    %d\n', ncols);
            fprintf(fileID,'nrows    %d\n', nrows);
            fprintf(fileID,'xllcorner    %f\n', xllcorner);
            fprintf(fileID,'yllcorner    %f\n', yllcorner);
            fprintf(fileID,'cellsize    %f\n', gridSize);
            fprintf(fileID,'NODATA_value    -9999\n');
            dlmwrite(filename,Z,'-append','delimiter','\t')
            fclose(fileID);
        end
        %%
        function [z_dem_changed,pix_IndValue] = AmendDEM(Z,R,lineShp,changeValue,varargin)
            %Z_New = AmendDEM(Z,R,lineShp,changeValue) change DEM value based on polyline
            %   Detailed explanation goes here
            % Z is DEM matrix, R is the reference matrix
            % lineShp is plolyline struct containing X and Y coords for line points
            % changeValue is an array containing the change values for each line
            % z_dem_changed = AmendDEM(Z,R,lineShp,changeValue,'Replace') exactly
            % replace original DEM with the changeValue.  It is default calculation choice
            % z_dem_changed = AmendDEM(Z,R,lineShp,changeValue,'Add') add a positive or
            % negtive value to the original DEM value.
            % z_dem_changed = AmendDEM(Z,R,lineShp,changeValue,'Positive') only
            % increase DEM values, if changeValue is lower than original, ignore the change
            
            % load defenceFailureData
            % defenceFailureData_ind = cell(length(defenceFailureData),2);
            
            % Created by Xiaodong Ming at 8/3/2018
            z_dem = Z; r_dem = R;
            changeLine_shp = lineShp;
            if isempty(varargin)
                calMethod = 'Replace';
            else
                calMethod = varargin{1};
            end
            %%find the points to be changed
            gridSize = abs(r_dem(2));
            z_dem_changed = z_dem;
            pix_IndValue = cell(length(changeLine_shp),1);
            for i = 1:length(changeLine_shp) % number of lines
                us_cre_lev = min(changeValue(i,:)); % changeValue upstream
                ds_cre_lev = max(changeValue(i,:)); % changeValue downstream
                x_vertex = [changeLine_shp(i).X]'; x_vertex(isnan(x_vertex)) = [];
                y_vertex = [changeLine_shp(i).Y]'; y_vertex(isnan(y_vertex)) = [];
                if size(x_vertex,1) == 1
                    x_vertex = x_vertex';
                end
                if size(y_vertex,1) == 1
                    y_vertex = y_vertex';
                end
                x_all = x_vertex;
                y_all = y_vertex;
                n = numel(x_all);
                % obtain xy coors from each line segment based on a unit of gridSize
                for j=1:numel(x_vertex)-1 %number of vertex in one polyline
                    x1 = x_vertex(j); x2 = x_vertex(j+1);
                    y1 = y_vertex(j); y2 = y_vertex(j+1);
                    if x1 ~= x2
                        xy_slope = (y2-y1)/(x2-x1);
                        xy_cosin = 1/sqrt(1+xy_slope^2);
                        gridGap = xy_cosin*gridSize*(x2-x1)/abs(x1-x2);
                        x_segment = x1:gridGap:x2;
                        y_segment = y1 + (x_segment-x1)*xy_slope;
                    else %x1==x2 means a vertical line segment
                        if y1==y2 %means they are duplicate points, not a line segment
                            continue
                        else
                            gridGap = gridSize*(y2-y1)/abs(y1-y2);
                            y_segment = y1:gridGap:y2;
                            x_segment = x1 + y_segment*0;
                        end
                    end
                    x_all(n+1:n+numel(x_segment)) = x_segment;
                    y_all(n+1:n+numel(x_segment)) = y_segment;
                    n = numel(x_all);
                end
                pix_RowCol = map2pix(r_dem,x_all, y_all); % cols and rows
                pix_RowCol = round(pix_RowCol);
                pix_RowCol(pix_RowCol<=0)=nan; % check the rows and cols out of Z range
                pix_RowCol(pix_RowCol(:,1)>size(z_dem,1),1)= nan;
                pix_RowCol(pix_RowCol(:,2)>size(z_dem,2),2)= nan;
                pix_Ind = sub2ind(size(z_dem),pix_RowCol(:,1),pix_RowCol(:,2));
                pix_Ind(isnan(pix_Ind)) = [];
                pix_Ind = unique(pix_Ind,'stable');
                %     pix_Ind = unique(pix_Ind,'stable'); % indice converted from cols and rows and sorted
                %**************change the value for each polyline**************************
                % from upstream to downstream
                if ~isempty(pix_Ind)
                    if z_dem_changed(pix_Ind(1))>=z_dem_changed(pix_Ind(end))
                        newValues = (linspace(us_cre_lev,ds_cre_lev,numel(pix_Ind)))';
                    else
                        newValues = (linspace(ds_cre_lev,us_cre_lev,numel(pix_Ind)))';
                    end
                    if strcmpi(calMethod,'Positive')
                        newValues = max(newValues,z_dem_changed(pix_Ind));
                    elseif strcmpi(calMethod,'Add')
                        newValues = newValues+z_dem_changed(pix_Ind);
                        %else Replace
                    end
                    z_dem_changed(pix_Ind) = newValues;
                    pix_IndValue{i} = [pix_Ind,newValues];
                end
            end
        end
    end
end
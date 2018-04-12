classdef HiPIMS_Case
    properties
        CaseFolder % dir of input folder
        BoundStruct %structure: Position([X Y]), Type('rigid','open'), T_h (time, water depth), T_uv (time, velocity)
        h0 = 0; % initial water depth
        uv0 = 0; % initial water velocity = [u0, v0]
        RainMask = 0;
        RainSources ={[0,0;1,0]} % cell, consisting of time series of rainfall rate (m/s)        
        precipitation_mask = 0;
        manning = 0.035; 
        sewer_sink = 0;
        hydraulic_conductivity = 0;
        cumul_depth = 0;
        capillary_head = 0;
        water_content_diff = 0;
        LanduseMask
        GaugeCoor = [0 0];        
    end
    properties (Access = private)
        Z % DEM matrix
        R % DEM spatial reference
        ID_ValidCell % the ID of valid cell starting from 0 at lower left corner and increasing by 1 rightwards and upwards
        ID_BoundCell % the ID of boundary cells, 1 for outline bound
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
        function obj = HiPIMS_Case(Z_dem,R_dem,varargin)
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
                'T_h',[0 0; 1 0],'T_uv',[0 0 0; 1 0 0]);
            
            obj.ID_BoundCell= GetBoundID(obj);
            obj.History = struct('DateTime',datestr(datetime('now')),'Log','Created');
            %struct('DateTime',datetime('now'),'Log','Created');
        end
        
        function obj = setBoundary(obj,boundStruct)
            obj.BoundStruct = boundStruct;
            obj.ID_BoundCell = GetBoundID(obj);
            % decompose flow into velocity x and y
            for i=1:length(obj.BoundStruct)
                T_uv = obj.BoundStruct(i).T_uv;
                if size(T_uv,2)==2 % uv is given as flow
                    flow = T_uv(:,2);
                    BoundWidth = abs(obj.R(2))*length(obj.ID_BoundCell(i).BoundIndAtZ);
                    if BoundWidth > 0
                        Q_per_m = flow/BoundWidth; % m^2/s
                        [I,J] = ind2sub(size(obj.Z),obj.ID_BoundCell(i).BoundIndAtZ); %range(cols):delta_x, range(a):delta_y
                        theta = atan(range(J)/range(I)); %theta: the angle between bound line and the vertical line
                        hu = Q_per_m*cos(theta); hv = Q_per_m*sin(theta); %m^2/s
                    else
                        hu = flow*0; hv = flow*0; % no inflow bound, width is zero
                    end
                    obj.BoundStruct(i).T_uv = [T_uv(:,1),hu,hv];
                end
            end
            
            NewHistory = obj.History;
            n=length(NewHistory);
            NewHistory(n+1).DateTime = datestr(datetime('now'));
            NewHistory(n+1).Log = 'Boundary changed';
            obj.History = NewHistory;
        end
        
        function obj = setRainfall(obj,rainMask,rainSources)
            obj.RainMask = rainMask;
            obj.RainSources = rainSources;
            NewHistory = obj.History;
            n=length(NewHistory);
            NewHistory(n+1).DateTime = datestr(datetime('now'));
            NewHistory(n+1).Log = 'Rainfall changed';
            obj.History = NewHistory;
        end
        
 
        %% write files        
        function output1 = writeInputFile(obj,varargin)
            if isempty(varargin)
                filenames = {'all'};
            else
                filenames = varargin; % file name without extension 
            end             
            [~, Locb] = ismember(filenames,obj.AllInputFileNames);
            Locb(Locb==0)=[]; 
            WriteFlag = Locb; % flag of files to be written
            creatFolders(obj)
            obj.FieldFolder = [obj.CaseFolder '/input/field/'];
            obj.MeshFolder = [obj.CaseFolder '/input/mesh/'];
            obj.ID_ValidCell = GetValidID(obj);
            if WriteFlag(1)==1 % all
                WriteFlag = 1:17;
            end
            output1 = obj.AllInputFileNames(WriteFlag); output1 = output1';
            for i = WriteFlag
                switch i
                    case 1 % write all input files
                        disp('Writing all input flies...')
                    case 2 % DEM
                        writeName = [obj.MeshFolder 'DEM.txt'];
                        Arcgridwrite(writeName,obj.Z,obj.R)
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
                    case 9 % cumul_depth
                        writeName = [obj.FieldFolder 'cumul_depth.dat'];
                        zValue = obj.cumul_depth+zeros(size(obj.Z));
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
                            if ~isempty(obj.BoundStruct(j).T_uv)
                            dlmwrite([obj.FieldFolder 'hU_BC_' num2str(n) '.dat'],...
                                obj.BoundStruct(j).T_uv,'precision',10,'delimiter',' ')
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
        
        %% visualization
        function Logs=GetHistory(obj)
            Logs = obj.History;
        end
        function GeneralMap(obj)
            [~] = PlotBoundID(obj);
            hold on
            mapshow(zeros(size(obj.Z)),obj.R,'Cdata',obj.Z,'DisplayType','surface')
            demcmap(obj.Z)
            lgd = legend;
            lgd.String(end)=[];
            hold off
        end
        function ValidID = PlotValidID(obj)
            mapshow(obj.ID_ValidCell,obj.R,'DisplayType','surface')
            axis image;
            ValidID = obj.ID_ValidCell;
%             title('Valid cells ID')
        end        
        function BoundID = PlotBoundID(obj)
            BoundID = obj.ID_BoundCell;
            hold on
            for i=1:length(BoundID)
                [row,col] = ind2sub(size(obj.Z),BoundID(i).BoundIndAtZ);
                [x,y] = pix2map(obj.R,row,col);
                if ~isempty(x)
                    if i==1
                        plot(x,y,'k.')
                    else
                        plot(x,y,'*')
                    end
                end
            end
            hold off
            axis image; grid on; box on
            lgd = legend;
            lgd.String{1} = 'Outline boundary';
            if length(lgd.String)>1
                for i=2:length(lgd.String)
                    lgd.String{i} = ['Boundary ', num2str(i-1)];
                end
            end
            lgd.Box = 'off'; %lgd.Location = 'best'; 
            lgd.Color = 'none';
        end
        
    end
    methods (Access = private)
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
            OutlineBound = find(OutlineBound);
            OutlineBoundLeft = OutlineBound;
            [row,col] = ind2sub(size(z_dem),OutlineBound);
            [xq,yq] = pix2map(r_dem,row,col); % XY coords of outline points
            numBound = length(boundStruct);
            if numBound==0 %outlint boundary, boundStruct is empty
                boundID.BoundIndAtZ = OutlineBound;
                boundID.Type = 'open'; %outlint boundary, boundStruct is 1
            elseif numBound==1 %outlint boundary
                boundID.BoundIndAtZ = OutlineBound;
                if isempty(boundStruct(1).Type)
                    boundID.Type = 'open';
                else
                    boundID.Type = boundStruct.Type;
                end
            else % there are IO boundaries, starting from 2
                for i=2:numBound
                    %IO boundary
                    X = boundStruct(i).Position(:,1);
                    Y = boundStruct(i).Position(:,2);
                    if sum(~isnan(X)) == 2 %only two points for bound region polygon
                        X = [min(X);max(X);max(X);min(X)];
                        Y = [min(Y);min(Y);max(Y);max(Y)];
                        boundStruct(i).Position = [X Y];
                    end
                    ind = inpolygon(xq,yq,X,Y);
                    boundID(i).BoundIndAtZ = OutlineBound(ind);
                    boundID(i).Type = boundStruct(i).Type;
                    OutlineBoundLeft(ind) = nan;
                end
                OutlineBoundLeft(isnan(OutlineBoundLeft)) = [];
                boundID(1).BoundIndAtZ = OutlineBoundLeft;
                if isempty(boundStruct(1).Type)
                    boundID(1).Type = 'open';
                else
                    boundID(1).Type = boundStruct.Type;
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
                    if isempty(boundStruct(i).T_uv)
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
        
        function writeFile_zType(obj,filename,zValue,TypeFlag)
            % TypeFlag = 'h','hU',[]?others?
            id_BoundCell = obj.ID_BoundCell; % struct
            id_bC = {id_BoundCell.BoundIndAtZ}'; % id cell2mat(
            % id & bound code
            for i=1:length(id_BoundCell)
                if strcmp(TypeFlag,'h')
                    if i>1 % add a minimum water depth for IO boundary
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
    end
end
function [outputArg1,outputArg2] = CombineMultiGPUResults(caseFolder,resultFileName)
%CombineMultiGPUResults Summary of this function goes here
%   Detailed explanation goes here
list = dir(caseFolder);
list(~[list.isdir])=[];
folderNum = nan(length(list),1);
for i=1:length(list)
    folderNum(i) = str2double(list(i).name);
end
folderNum_max = max(folderNum);
Z = []; z1 = Z;
gauge_Ind1 = [];
gauge_Value1 = [];
gauge_time_Value = [];
gauge_Ind = [];
for i=folderNum_max:-1:0
    if strcmp(resultFileName(end-2:end),'asc')
        oneFileName = [caseFolder '\'  num2str(i) '\output\' resultFileName];
        if i==folderNum_max
            [z0,R] = arcgridread(oneFileName);
        else
            [z1,~] = arcgridread(oneFileName);
            z1(1:2,:) = [];
        end
        if isempty(z1)
            Z = z0;
        else
            Z = [z0;z1];
        end
        z0 = Z;
    elseif strcmp(resultFileName(end-2:end),'txt')
        oneFileName = [caseFolder '\'  num2str(i) '\input\mesh\' resultFileName];
        if i==folderNum_max
            [z0,R] = arcgridread(oneFileName);
        else
            [z1,~] = arcgridread(oneFileName);
            z1(1:2,:) = [];
        end
        if isempty(z1)
            Z = z0;
        else
            Z = [z0;z1];
        end
        z0 = Z;
    elseif length(resultFileName)>10&&strcmp(resultFileName(end-9:end),'gauges.dat')
        % read DEM XY range & gauges position
        fileID = fopen([caseFolder '/' num2str(i) '/input/mesh/DEM.txt']);
        M = textscan(fileID,'%s %f\n',5); M = M{2};
        domainRectangle = [M(3),M(4);...
            M(3),M(4)+M(5)*M(2);...
            M(3)+M(5)*M(1),M(4)+M(5)*M(2);...
            M(3)+M(5)*M(1),M(4)];
        fclose(fileID);
        gauges_XY = dlmread([caseFolder '/' num2str(i) '/input/field/gauges_pos.dat']);
        gauges_Index = dlmread([caseFolder '/' num2str(i) '/input/field/gauges_index.dat']);
        in = inpolygon(gauges_XY(:,1),gauges_XY(:,2),...
            domainRectangle(:,1),domainRectangle(:,2));
        gauge_Ind0 = gauges_Index(in);
        % read output value in gauges
        gauges_output = dlmread([caseFolder '/'  num2str(i) '/output/' resultFileName]);
        timeS = gauges_output(:,1);
        if strcmp(resultFileName(1:2),'hU')
            valueS = gauges_output(:,2:end)+nan;
            if ~isempty(valueS)
            valueSx = gauges_output(:,2:2:end);
            valueSx = valueSx(:,in');
            valueSy = gauges_output(:,3:2:end);
            valueSy = valueSy(:,in');
            if isempty(valueSx)
                valueS = [];
            else
                valueS(:,1:2:end) = valueSx;
                valueS(:,2:2:end) = valueSy;
            end
            end
        else
            valueS = gauges_output(:,2:end);
            valueS = valueS(:,in');
        end
        if isempty(gauge_Ind1)
            gauge_Ind = gauge_Ind0;
        else
            gauge_Ind = [gauge_Ind1;gauge_Ind0];
        end
        gauge_Ind1 = gauge_Ind;
        if isempty(gauge_Value1)
            gauge_Value1 = [timeS valueS];
            gauge_time_Value = gauge_Value1;
        else
            gauge_time_Value = [gauge_Value1 valueS];
        end
        gauge_Value1 = gauge_time_Value;
        fprintf('%u ',i)
    end
end
if ~isempty(gauge_time_Value)
    [gauge_Ind,Isort,~] = unique(gauge_Ind);
    if strcmp(resultFileName(1:2),'hU')
        gauge_Value = gauge_time_Value(:,2:end)+nan;
        gauge_ValueX = gauge_time_Value(:,2:2:end);
        gauge_ValueY = gauge_time_Value(:,3:2:end);
        gauge_ValueX = gauge_ValueX(:,Isort');
        gauge_ValueY = gauge_ValueY(:,Isort');
        gauge_Value(:,1:2:end) = gauge_ValueX;
        gauge_Value(:,2:2:end) = gauge_ValueY;
    else
        gauge_Value = gauge_time_Value(:,2:end);
        gauge_Value = gauge_Value(:,Isort');
    end
    outputArg1 = [gauge_time_Value(:,1) gauge_Value];
    outputArg2 = gauge_Ind;
elseif ~isempty(Z)
    outputArg1 = Z;
    outputArg2 = R;
end
end


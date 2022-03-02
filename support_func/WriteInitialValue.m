function WriteInitialValue(writePath,fileName,Cell_ID,Bound_ID,InitialValue,BoundType)
% default bound type : 1: normal bound; 2: input bound; 3: output bound
% Bound_ID: -1: shared bound; 1: outline bound; 2,3,4,...: IO bound
% Cell_ID: nan, 0 ~ number_of_valid_cell-1
% fileNames: 'h','eta','hU' are three special files, others including
% 'z','precipitation','manning'...
% BoundType_Vec is cell containing 2*3 matrix  [h;hU]*{bound1 2 3}
adjustValue = 0.001; % avoid all the value to be zero in 'h', 'eta', 'hU';
Outline_Bound_Value = 2; %[2 0 0]
Shared_Bound_Value = 4; %[4 0 0]
BoundType_Vec = cell(1,size(BoundType,2)/3);
for i = 1:length(BoundType_Vec)
    BoundType_Vec{i} = BoundType(:,(1:3)+(i-1)*3);
end

%% ****ValidID_InitialValue: valid cells' ID and their initial value
if strcmp(fileName,'h')||strcmp(fileName,'eta') % avoid the depth in the flow-input-bound cells to be nil
    InitialValue(Bound_ID==2) = InitialValue(Bound_ID==2)+adjustValue;
end 

ValidCellValue = InitialValue(:); %Valid cells' value
if strcmp(fileName,'hU') % file name is hU, hU consist of [hu hv]
    ValidID_InitialValue = [Cell_ID(:), ValidCellValue(1:end/2), ValidCellValue(end/2+1:end)];
    formatSpec1 = '%-12d %.8f %.8f\n';    % format for valid cells' ID and their initial values 
else  %if filename is z, h, precepitation, or manning   
    ValidID_InitialValue = [Cell_ID(:), ValidCellValue];
    formatSpec1 = '%-12d %.8f\n'; % format for valid cells' ID and their initial values
end
ValidID_InitialValue(isnan(ValidID_InitialValue(:,1))|ValidID_InitialValue(:,1)<0,:) = []; % remove the invalid cells and their value
ValidID_InitialValue(isnan(ValidID_InitialValue)) = 0; % replace the NaN initial value to zero
ValidID_InitialValue = sortrows(ValidID_InitialValue,1);

%% ****BondType_Vector_Print: boundary cells' valid ID & their boundary type vec
formatSpec2 = '%-12d %d %d %d\n'; % format for bound cells' ID and their boundary vector
Valid_Bound_ID = [Cell_ID(:),Bound_ID(:)]; 
Valid_Bound_ID(isnan(Valid_Bound_ID(:,1)),:) = []; % remove the invalid cells
Valid_Bound_ID(Valid_Bound_ID(:,2)==0,:)=[]; % remove the non-boundary cells
Valid_Bound_ID = sortrows(Valid_Bound_ID,1); % only boundary cells left

bound_cnt = length(Valid_Bound_ID); % number of bound cells
BondCell_TypeVector = zeros(bound_cnt,3); % default Value [0 0 0]
% default outline bound type vector [2 0 0]
BondCell_TypeVector(Valid_Bound_ID(:,2)== 1,1) = Outline_Bound_Value;
% default shared bound type vector [4 0 0]
BondCell_TypeVector(Valid_Bound_ID(:,2)==-1,1) = Shared_Bound_Value;

if strcmp(fileName,'h')|| strcmp(fileName,'eta')||strcmp(fileName,'hU')
    if strcmp(fileName,'h')|| strcmp(fileName,'eta'); i = 1;
    elseif strcmp(fileName,'hU'); i = 2;
    else error('Please check the fileName of h, eta and hU')
    end
    % boundary type vector for [h; hU]* [bound1, bound2, bound3]
    for m = 1:length(BoundType_Vec) %m: number of bound types
        ind = (Valid_Bound_ID(:,2) == m); % find the the corresponding bound ID 
        BondCell_TypeVector(ind,:) = repmat(BoundType_Vec{m}(i,:),sum(ind),1);
    end    
end
Bond_ValidID_TypeVector = [Valid_Bound_ID(:,1), BondCell_TypeVector];

%% print files
% if the value in 'h.dat' or 'hU.dat' is pure nil,add a small positive value to the first element
if strcmp(fileName(1),'h') && sum(ValidID_InitialValue(:,2))==0 
    ValidID_InitialValue(1,2) = adjustValue;
end
fileID = fopen([writePath fileName '.dat'],'w');
% print valid cell and their initial value
fprintf(fileID,'$Element Number\n'); fprintf(fileID,'%d\n',length(ValidID_InitialValue));
fprintf(fileID,'$Element_id  Value\n'); fprintf(fileID,formatSpec1,ValidID_InitialValue');
% print boundary cell and their vector
fprintf(fileID,'$Boundary Numbers\n'); fprintf(fileID,'%d\n',bound_cnt);
fprintf(fileID,'$Element_id Boundary_type\n'); fprintf(fileID,formatSpec2,Bond_ValidID_TypeVector');
fclose(fileID);
end

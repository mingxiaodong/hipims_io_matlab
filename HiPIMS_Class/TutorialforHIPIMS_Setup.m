clc
load('TestData.mat', 'BoundStruct')
[Z,R] = arcgridread('TestDEM.txt');
% create case
caseObj = HiPIMS_Setup(Z,R); % casefolder is the current folder
caseObj = HiPIMS_Setup(Z,R,'C:\Users\b4042552\Dropbox\Matlab\HiPIMS_Class');
caseObj = setBoundary(caseObj,BoundStruct);
%%
fig = figure; fig.Color = 'w';
GeneralMap(caseObj)
% A= writeInputFile(caseObj,'hU');
%%
BoundID = PlotBoundID(caseObj);
IND = {BoundID.BoundIndAtZ}';
IND = cell2mat(IND);
[row,col] = ind2sub(size(Z),IND);
[X,Y] = pix2map(R,row,col);
%% sort the boundary points
clc
x0 = X(1); y0 = Y(1);
IND = nan(size(X)); IND(1) = 1;
XY0 = [X Y]; XY0(1,:)=nan;
for i=2:numel(X)
    ed = (XY0(:,1)-x0).^2+(XY0(:,2)-y0).^2;
    M = min(ed);
    if ~isnan(M)
        I = find(ed==M);
        IND(i) = I(end);
        x0 = X(I(end)); y0 = Y(I(end));
        XY0(I,:)=nan;
    end
end
IND(isnan(IND)) = [];
reorderedx = X(IND);
reorderedy = Y(IND);
plot(reorderedx, reorderedy, 'bo-');
axis equal
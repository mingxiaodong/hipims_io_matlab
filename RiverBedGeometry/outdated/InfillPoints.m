function outPutPoints = InfillPoints(inputPoints,infillTimes)
% linear infill points(centre)
if nargin == 1
    infillTimes=1;
end
for i = 1:infillTimes
    inputPoints_half = (inputPoints(1:end-1,:)+inputPoints(2:end,:))/2;
    newX = [inputPoints(1:end-1,1),inputPoints_half(:,1)];
    newY = [inputPoints(1:end-1,2),inputPoints_half(:,2)];
    newX = newX'; newX = newX(:);
    newY = newY'; newY = newY(:);
    if size(inputPoints,2)==3
        newZ = [inputPoints(1:end-1,3),inputPoints_half(:,3)];
        newZ = newZ'; newZ = newZ(:);
        outPutPoints = [[newX,newY,newZ];inputPoints(end,:)];
    else
        outPutPoints = [[newX,newY];inputPoints(end,:)];
    end
    inputPoints = outPutPoints;
end
end
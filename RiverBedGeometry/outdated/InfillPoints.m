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
    outPutPoints = [[newX,newY];inputPoints(end,:)];
    inputPoints = outPutPoints;
end
end
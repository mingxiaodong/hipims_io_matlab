function ValueSum = IndicationVector(varargin)
% get the sum of indication function value between vector and scalar
% ValueSum = IndicationVector(Vec1,s1);
if mod(length(varargin),2)~=0
    error('wrong number of inputs')
end
nCon = length(varargin)/2;
Ind0 = 0;
for i=1:nCon
    Vec1 = varargin{i*2-1};
    s1 = varargin{i*2};
    s1 = repmat(s1,size(Vec1));
    Ind1 = Vec1<=s1;
    if i==1
        Ind0 = Ind1;
    else
        Ind0 = and(Ind0,Ind1);
    end
end
ValueSum = sum(Ind0);
end
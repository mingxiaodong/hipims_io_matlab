function [points_l,points_w] = RelativeCoordsInterp(bank1Points_xy,bank1Points_lw,...
    bank2Points_xy,bank2Points_lw,points_xy)
% x & y: matrix, each row convey the x & y coordinate of the points to be
% interpolated for point [xq(i),yq(i)]
% v_L & v_W: vector with the same length of xq, convey L and W values of
% the points to be interpolated
% xq,yq: the x and y coordinates of the targeting points

[x_interp1,y_interp1,l_interp1,w_interp1] = Nearest3BankPoints...
    (bank1Points_xy,bank1Points_lw,points_xy);% column: num of points to be interp
% bank2Points_xy = bank2Points_xy(end:-1:1,:);
[x_interp2,y_interp2,l_interp2,w_interp2] = Nearest3BankPoints...
    (bank2Points_xy,bank2Points_lw,points_xy);% row: num ofcrossline points
x_interp = [x_interp1,x_interp2];
y_interp = [y_interp1,y_interp2];
l_interp = [l_interp1,l_interp2];
w_interp = [w_interp1,w_interp2];

xq = points_xy(:,1);
yq = points_xy(:,2);
coor_L = nan(length(xq),1);
coor_W = coor_L;
for i=1:length(xq)
    x_ToBeInterp = x_interp(i,:)';
    y_ToBeInterp = y_interp(i,:)';
    v1_ToBeInterp = l_interp(i,:)';
    v2_ToBeInterp = w_interp(i,:)';
    weight = -2;
    coor_L(i) = IDW(x_ToBeInterp,y_ToBeInterp,v1_ToBeInterp,...
        xq(i),yq(i),weight,'ng',length(x_ToBeInterp));
    coor_W(i) = IDW(x_ToBeInterp,y_ToBeInterp,v2_ToBeInterp,...
        xq(i),yq(i),weight,'ng',length(x_ToBeInterp));
end
points_l = coor_L;
points_w = coor_W;

end

function[Vint]=IDW(xc,yc,vc,x,y,e,r1,r2)
%%%%%%%%%%%%%%%%%
%%% INPUTS
%xc = stations x coordinates (columns) [vector]
%yc = stations y coordinates (rows) [vector]
%vc = variable values on the point [xc yc]
%x = interpolation points  x coordinates [vector]
%y = interpolation points y coordinates [vector]
%e = distance weight
%r1 --- 'fr' = fixed radius ;  'ng' = neighbours
%r2 --- radius lenght if r1 == 'fr' / number of neighbours if  r1 =='ng'
%%% OUTPUTS
%Vint --- Matrix [length(y),length(x)] with interpolated  variable values
%%% EXAMPLES
%%% --> V_spa=IDW(x1,y1,v1,x,y,-2,'ng',length(x1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simone Fatichi -- simonef@dicea.unifi.it
%   Copyright 2009
%   $Date: 2009/06/19 $
%   $Updated 2012/02/24 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vint=zeros(length(y),length(x));
xc=reshape(xc,1,length(xc));
yc=reshape(yc,1,length(yc));
vc=reshape(vc,1,length(vc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  strcmp(r1,'fr')
    if  (r2<=0)
        disp('Error: Radius must be positive')
        return
    end
    for i=1:length(x)
        for j=1:length(y)
            D=[]; V=[]; wV =[]; vcc=[];
            D= sqrt((x(i)-xc).^2 +(y(j)-yc).^2);
            if min(D)==0
                disp('Error: One or more stations have the coordinates of an interpolation point')
                return
            end
            vcc=vc(D<r2); D=D(D<r2);
            V = vcc.*(D.^e);
            wV = D.^e;
            if isempty(D)
                V=NaN;
            else
                V=sum(V)/sum(wV);
            end
            Vint(j,i)=V;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    if (r2 > length(vc)) || (r2<1)
        disp('Error: Number of neighbours not congruent with data')
        return
    end
    for i=1:length(x)
        for j=1:length(y)
            D=[]; V=[]; wV =[];vcc=[];
            D= sqrt((x(i)-xc).^2 +(y(j)-yc).^2);
            if min(D)==0
                disp('Error: One or more stations have the coordinates of an interpolation point')
                return
            end
            [D,I]=sort(D);
            vcc=vc(I);
            V = vcc(1:r2).*(D(1:r2).^e);
            wV = D(1:r2).^e;
            V=sum(V)/sum(wV);
            Vint(j,i)=V;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
end
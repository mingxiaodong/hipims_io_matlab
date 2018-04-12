%% three dimensional Gumbel?Hougaard copula
%u1,u2,u3: the cdf of three variables
%si: theta, the parameter of Gumbel?Hougaard copula
%C3_Gum: the 3D Gumbel-Hougaard copula
clear,clc

%%load data with three variables
load('CountyEvents0921.mat')
i = 114; % any county
% Source data of X, Y, Z
Xs = [CountyEvent{i,1}.TotalRain];
Ys = [CountyEvent{i,1}.MaxRain];
Zs = [CountyEvent{i,1}.MaxWind];
CoefX = gevfit(Xs);
CoefY = gevfit(Ys);
CoefZ = gevfit(Zs);
step = 50; 
p1=0.000001;
p2=0.9;
Xrange = gevinv([p1,p2],CoefX(1),CoefX(2),CoefX(3));
Yrange = gevinv([p1,p2],CoefY(1),CoefY(2),CoefY(3));
Zrange = gevinv([p1,p2],CoefZ(1),CoefZ(2),CoefZ(3));

x = linspace(Xrange(1),Xrange(2),step);%rand(100,1));
y = linspace(Yrange(1),Yrange(2),step);%rand(100,1);
z = linspace(Zrange(1),Zrange(2),step);%rand(100,1);
%pz = linspace(0.01,0.9,step)';
%z = gevinv(pz,CoefZ(1),CoefZ(2),CoefZ(3));

[X,Y,Z] = meshgrid(x,y,z);
u1 = gevcdf(X(:),CoefX(1),CoefX(2),CoefX(3));
u2 = gevcdf(Y(:),CoefY(1),CoefY(2),CoefY(3));
u3 = gevcdf(Z(:),CoefZ(1),CoefZ(2),CoefZ(3));

%copula
si = 2.69;
C3_Gum = exp(...
    -(...
    (-log(u1)).^si+(-log(u2)).^si+(-log(u3)).^si...
    ).^(1/si)...
    );
C3_Gum = reshape(C3_Gum,size(X));

%%dynamaic plot of 3D joint CDF
curFolder = cd;
figure
for l = 1:step
    XYx = X(:,:,l);
    XYy = Y(:,:,l);
    surf(XYx,XYy,C3_Gum(:,:,l))
    xlabel('X')
    ylabel('Y')
    zlabel('F_{xyz} - Gumbel Joint CDF')
    Pmax = C3_Gum(length(z),length(z),l);
    axis([Xrange,Yrange,0,0.9]) 
    title(['When Z = num2str F_z = ',num2str(gevcdf(z(l),CoefZ(1),CoefZ(2),CoefZ(3)),'%.2f')])
    pause(0.4)
    PictureName = [curFolder '\pictures\' num2str(l,'%03u')];
    print(PictureName,'-djpeg','-r100')
end
%% make animation
cd([curFolder '/pictures'])
AnimationName = [curFolder '/Copula3D.gif'];
list = dir('*.jpg');
N = size(list);
for i=1:N
    list(i).datenum = str2double(list(i).name(1:end-4));
end
[~,sortID] = sortrows([list.datenum].'); list = list(sortID); clear sortID
for i = 1:N
    img = imread(list(i).name);
    imshow(img);
    frame=getframe(gcf);
    im=frame2im(frame);
    [A,map]=rgb2ind(img,256);
    if i == 1
		imwrite(A,map,AnimationName,'gif','LoopCount',Inf,'DelayTime',0.2);
	else
		imwrite(A,map,AnimationName,'gif','WriteMode','append','DelayTime',0.2);
    end
end


clc
clear
n = 10;
for ii = 1 : n
    [best_fuzzy,best_topsis] = imopsofunc();
    best_f{ii} = best_fuzzy;
    best_t{ii} = best_topsis;
end

best_f
best_t

for i = 1:10
    plot(best_f(1,i),best_f(2,i),'bs');
    plot(best_t(1,i),best_t(2,i),'ko');
    hold on
end
x0 = rep_costs(1,:)';
y0 = rep_costs(2,:)';
z0 = rep_costs(3,:)';
uij = 30:5:45; 
vij=0:0.05:0.2;  %这是定义自变量u，v的范围
[xi,yj]=meshgrid(uij,vij); %这是画出网格点
z2=griddata(x0,y0,z0,xi,yj,'cubic');%这是作立方插值运算，你的坐标为x0,y0,z0,
mesh(xi,yj,z2); %这就开始作曲面图了
hold on;%保持图形不变
plot3(x0,y0,z0,'mo');%这作的是你坐标的散点图
hold off;


n = 100000;
delta1 = rand(n, 1)*2*pi;
delta2 = rand(n, 1)*2*pi;
delta3 = rand(n, 1)*2*pi;
delta4 = rand(n, 1)*2*pi;

x = cos(delta4) - cos(delta2);
y = cos(delta1) - cos(delta3);
z = sin(delta1) + sin(delta2) + sin(delta3) + sin(delta4);
scatter3(x, y, z, '.y')    
hold on
%%
[x, y] = deal(linspace(-2, 2, 120));
zMax = 2*sin(acos(x/2));
[X, Y] = meshgrid(x, y);
Z = bsxfun(@plus, zMax, zMax');
surf(X, Y, Z)
colormap hot
surf(X, Y, -Z)

%%

x = rep_costs(1,:);
y = rep_costs(2,:);
z = rep_costs(3,:);
scatter3(x, y, z, '.y')    
hold on
%%
% [x, y] = deal(linspace(-2, 2, 120));
x = linspace(30, 45, 100);
y = linspace(0, 0.02, 100);
zMax = 2*sin(acos(x/2));
[X, Y] = meshgrid(x, y);
Z = bsxfun(@plus, zMax, zMax');
surf(X, Y, Z)
colormap hot
surf(X, Y, -Z)


%%
x0 = rep_costs(1,:)';
y0 = rep_costs(2,:)';
z0 = rep_costs(3,:)';

X = rep_costs';
C=convhulln(X);
C=C(:);   %将C变换成一列
C=unique(C);
XI=X(C,:);%
figure
scatter3(XI(:,1),XI(:,2),XI(:,3),'.B');
hold off
figure
scatter3(XI(:,1),XI(:,2),XI(:,3),'.B');
tu_x=XI(:,1);%包络面上点的X坐标
tu_y=XI(:,2);%包络面上点的Y坐标
tu_z=XI(:,3);%包络面上点的Z坐标
index=find(tu_z>=0);
tu_xb=tu_x(index);
tu_yb=tu_y(index);
tu_zb=tu_z(index);
[fitobject,gof,output]=fit([tu_xb,tu_yb],tu_zb,'poly23');
plot(fitobject,[tu_xb,tu_yb],tu_zb)
[fitobject1,gof1,output1]=fit([tu_xb,tu_yb],tu_zb,'poly33');
plot(fitobject1,[tu_xb,tu_yb],tu_zb)
hold off
figure
[fitobject2,gof2,output2]=fit([tu_xb,tu_yb],tu_zb,'smoothingspline');
plot(fitobject2,[tu_xb,tu_yb],tu_zb)
saveas(gcf,'C:\Users\Feng Da\Documents\MATLAB\AGCAVC-IMOPSO\results\39bus-100-500-200-3.fig');

tri = delaunay(X(:,1),X(:,2));
hold off
figure
trisurf(tri,X(:,1),X(:,2),X(:,3));
hidden off
shading interp
hold on 
plot3(rep_costs(1,idx),rep_costs(2,idx),rep_costs(3,idx),'ro');
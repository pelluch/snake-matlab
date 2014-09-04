clc;
close all;
% clear all;

P = [20,20;20,40;25,45;30,45;50,35;50,30;40,30;40,35;60,50;60,20];

[P2, P3] = convhull(P(:,1)',P(:,2)');
P4 = MonotoneConvHull(P);

xx = -1:.05:1;
yy = abs(sqrt(xx));
[x, y] = pol2cart(xx,yy);
k = convhull(x,y);
figure;
plot(x(k),y(k),'r-',x,y,'b+');

figure;
plot(P(P2,1),P(P2,2),'r-', P(:,1), P(:,2),'b+');


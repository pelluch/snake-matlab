% test Resample Polygon
clc;
clear all;
close all;

x = [0,0;0,10];
maxdist = 200;

y = ResamplePolygon(x, maxdist);
y = [y ; x(end,:)];
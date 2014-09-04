% run DoubleSnake for different datasets
%% Load data
clc;
close all;


% clear all;
% Read VTK 3D
for ind = 1:12
    % for ind = 1:24
    [V(:,:,:,ind), info] = ReadData3D(strcat('C:/DocsMaracuya/Datos/Snakes/SAguirre/_OriginalReductor_',int2str(ind-1),'.vtk'));
    % [V(:,:,:,ind), info] = ReadData3D(strcat('C:/epxImages/SAguirre/KarisEjeCorto/_OriginalReductor_',int2str(ind-1),'.vtk'));
end
%% Parameters

segLimits.startS=4;
segLimits.endS=13;
segLimits.startF = 11;
segLimits.endF = 11;
segLimits.dir = 1;        % 0 = base to apex

circle.center = [ 140, 160 ] ;
circle.radius = 16;
circle.numPoints = 100;

P1Options.Gamma = 1; % Time step, default 1
P1Options.Iterations = 500; %  Number of iterations, default 100
P1Options.Sigma1 = 3; %  Sigma used to calculate image derivatives, default 10
P1Options.Wline = 0; %  Attraction to lines, if negative to black lines otherwise white lines , default 0.04
P1Options.Wedge = 2; %  Attraction to edges, default 2.0
P1Options.Wterm = 0; % Attraction to terminations of lines (end points) and corners, default 0.01
P1Options.Sigma2 = 2; % Sigma used to calculate the gradient of the edge energy image (which gives the image force), default 20
% options (Gradient Vector Flow)
P1Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors, default 0.2. (Warning setting this to high >0.5 gives an instable Vector Flow)
P1Options.GIterations = 1; % Number of GVF iterations, default 0
P1Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
P1Options.Alpha = .1; % Membrame energy  (first order), default 0.2
P1Options.Beta = 1; %  Thin plate energy (second order), default 0.2
P1Options.Delta = +0.009; %  Baloon force, default 0.1
P1Options.Kappa = .001; %  Weight of external image force, default 2
P1Options.r = 1; %  Weight of spring
% P1Options.r = 0.; %  Weight of spring

% Second snake

P2Options.minWidth = 10; %  Distance between two snakes in pixels.
P2Options.maxWidth = 16; %  Distance between two snakes in pixels.
P2Options.medianKernel = 4; % Median filter before clustering kernel
P2Options.minTopCluster = 3; % First class to be joined, up to numClust
P2Options.numClust = 4; % Number of kmeans clusters


P2Options.nPoints = 100; %: Number of contour points, default 100
P2Options.Gamma = .25; % Time step, default 1
P2Options.Iterations = 200; %  Number of iterations, default 100
P2Options.Sigma1 = 1; %  Sigma used to calculate image derivatives, default 10
P2Options.Wline = 0; %  Attraction to lines, if negative to black lines otherwise white lines , default 0.04
P2Options.Wedge = 2; %  Attraction to edges, default 2.0
P2Options.Wterm = 0; % Attraction to terminations of lines (end points) and corners, default 0.01
P2Options.Sigma2 = 1; % Sigma used to calculate the gradient of the edge energy image (which gives the image force), default 20
P2Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors, default 0.2. (Warning setting this to high >0.5 gives an instable Vector Flow)
P2Options.GIterations = 0; % Number of GVF iterations, default 0
P2Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
P2Options.Alpha = .4; % Membrame energy  (first order), default 0.2
P2Options.Beta = 4; %  Thin plate energy (second order), default 0.2
P2Options.Delta = .1; %  Baloon force, default 0.1
P2Options.Kappa = .7; %  Weight of external image force, default 2
P2Options.WidthWeight = .04; %  Weight of distance between snakes, default .04


Out.writeFolder = 'C:/DocsMaracuya/Datos/Snakes/DoubleSnake_SergioAguirre';
Out.x1 = 100;
Out.x2 = 180;
Out.y1 = 120;
Out.y2 = 200;

z = size(V,3);
t = size(V,4);

P1M = cell(z,t);
P2M = cell(z,t);
SpringP1 = cell(z,t);
SpringP2 = cell(z,t);

% SpringP1{9,11} = [133, 185];
SpringP1{9,11} = [147, 179];
% SpringP1{9,11} = [179, 147];

minSize = 50;


%% Function only


[P0, P1, P2, borderClass, Eext] = DoubleSnake(V, segLimits, circle, minSize, P1Options, P2Options, SpringP1, SpringP2);



PlotSnakes( P0, P1, P2, P1M, P2M, SpringP1, SpringP2 ,V, borderClass, Eext, Out, segLimits, true);



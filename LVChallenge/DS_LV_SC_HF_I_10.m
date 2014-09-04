% run DoubleSnake for different datasets
%% Load data
clc;
close all;


% clear all;
% clear V,D;

% Example to read dicom

filename = 'C:/DocsMaracuya/Datos/LVChallenge/ShortAxis/SC-HF-I-10/IM-0024-0001.dcm' ;
[D, info] = ReadData3D(filename);
t = info.CardiacNumberOfImages;
z = size(D,3)/t;
% imlook3d(D);

for i = 1:z
    for j = 1:t
        V(:,:,i,j) = D(:,:,(i-1)*t+(j-1)+1);
    end
end
% for i=1:z
%     imlook3d(V(:,:,i,:));
% end
code = 'SC-HF-I-10';

Out.writeFolder = strcat('C:/DocsMaracuya/Datos/Snakes/DoubleSnake_LV_',code);
fold = 'C:\DocsMaracuya\Datos\LVChallenge\Sunnybrook_Cardiac_MR_Database_ContoursPart1\Sunnybrook Cardiac MR Database ContoursPart1\OnlineDataContours\';
subfolder = strcat(code,'\contours-manual\IRCCI-expert\');
fname = strcat('list-',code,'.txt');
fname2 = strcat(fold,fname);

[P1M, P2M] = LoadManualSeg(fold, subfolder, fname, z, t);

%% Parameters

segLimits.startS=1;
segLimits.endS=10;
segLimits.startF = 7;
segLimits.endF = 20;
segLimits.dir = 1;        % 0 = base to apex

circle.center = [ 133, 124 ] ;
circle.radius = 30;
circle.numPoints = 6*circle.radius;

P1Options.Gamma = 1; % Time step, default 1
P1Options.Iterations = 500; %  Number of iterations, default 100
P1Options.Sigma1 = 1; %  Sigma used to calculate image derivatives, default 10
P1Options.Wline = 0; %  Attraction to lines, if negative to black lines otherwise white lines , default 0.04
P1Options.Wedge = 2; %  Attraction to edges, default 2.0
P1Options.Wterm = 0; % Attraction to terminations of lines (end points) and corners, default 0.01
P1Options.Sigma2 = 2; % Sigma used to calculate the gradient of the edge energy image (which gives the image force), default 20
% options (Gradient Vector Flow)
P1Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors, default 0.2. (Warning setting this to high >0.5 gives an instable Vector Flow)
P1Options.GIterations = 0; % Number of GVF iterations, default 0
P1Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
P1Options.Alpha = .3; % Membrame energy  (first order), default 0.2
P1Options.Beta = 2; %  Thin plate energy (second order), default 0.2
P1Options.Delta = +0.008; %  Baloon force, default 0.1
P1Options.Kappa = .002; %  Weight of external image force, default 2


% Second snake

P2Options.minWidth = 4; %  Distance between two snakes in pixels.
P2Options.maxWidth = 11; %  Distance between two snakes in pixels.
P2Options.medianKernel = 5; % Median filter before clustering kernel
% P2Options.minTopCluster = 3; % First class to be joined, up to numClust
P2Options.minTopCluster = 4; % First class to be joined, up to numClust
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
P2Options.Sigma3 = 2.; %  Sigma used to calculate the laplacian in GVF, default 1.0
P2Options.Alpha = .4; % Membrame energy  (first order), default 0.2
P2Options.Beta = 4; %  Thin plate energy (second order), default 0.2
P2Options.Delta = .1; %  Baloon force, default 0.1
P2Options.Kappa = .7; %  Weight of external image force, default 2
P2Options.WidthWeight = .04; %  Weight of distance between snakes, default .04


Out.x1 = circle.center(1)-50;
Out.x2 = circle.center(1)+50;
Out.y1 = circle.center(2)-50;
Out.y2 = circle.center(2)+50;



minSize = 50;


%% Function only

[P0, P1, P2, borderClass, Eext] = DoubleSnake(V, segLimits, circle, minSize, P1Options, P2Options);

PlotSnakes( P0, P1, P2, P1M, P2M, V, borderClass, Eext, Out, segLimits);
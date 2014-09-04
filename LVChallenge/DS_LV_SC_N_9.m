% run DoubleSnake for different datasets
%% Load data
clc;
close all;


% clear all;

% Example to read dicom

filename = 'C:/DocsMaracuya/Datos/LVChallenge/ShortAxis/SC-N-9/IM-1031-0001.dcm' ;
[D, info] = ReadData3D(filename);
t = info.CardiacNumberOfImages;
z = size(D,3)/t;

for i = 1:z
    for j = 1:t
        V(:,:,i,j) = D(:,:,(i-1)*t+(j-1)+1);
    end
end

fold = 'C:\DocsMaracuya\Datos\LVChallenge\Sunnybrook_Cardiac_MR_Database_ContoursPart1\Sunnybrook Cardiac MR Database ContoursPart1\OnlineDataContours\';
subfolder = 'SC-N-09\contours-manual\IRCCI-expert\';
fname = 'list-SC-N-09.txt';
fname2 = strcat(fold,fname);

[P1M, P2M] = LoadManualSeg(fold, subfolder, fname, z, t);


%% Parameters

segLimits.startS=2;
segLimits.endS=9;
segLimits.startF = 07;
segLimits.endF = 20;
segLimits.dir = 1;        % 0 = base to apex

circle.center = [ 103, 128 ] ;
circle.radius = 19;
circle.numPoints = 6*circle.radius;

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


% Second snake

P2Options.minWidth = 5; %  Distance between two snakes in pixels.
P2Options.maxWidth = 10; %  Distance between two snakes in pixels.
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


Out.writeFolder = 'C:/DocsMaracuya/Datos/Snakes/DoubleSnake_LV_SC_N_9_spring';
Out.x1 = circle.center(1)-50;
Out.x2 = circle.center(1)+50;
Out.y1 = circle.center(2)-50;
Out.y2 = circle.center(2)+50;


minSize = 50;

SpringP1 = cell(z,t);
SpringP1{2,20} = [85, 105];
SpringP1{5,19} = [87, 116];
SpringP1{5,17} = [87, 116];
SpringP1{5,18} = [86, 117];
SpringP1{6,20} = [100, 124; 102, 114];
SpringP1{9,20} = [110, 123; 124, 124];
SpringP1{9,19} = [110, 123; 124, 124];
SpringP1{9,18} = [110, 123; 124, 124];
SpringP1{9,17} = [110, 123; 124, 124];
SpringP1{9,16} = [109, 120; 124, 124];
SpringP1{9,12} = [115, 122; 124, 124];



SpringP2 = cell(z,t);
SpringP2{2,14} = [72, 129];
SpringP2{5,18} = [81, 117];
SpringP2{7,18} = [98, 117];

%% Function only

[P0, P1, P2, borderClass, Eext] = DoubleSnake(V, segLimits, circle, minSize, P1Options, P2Options, SpringP1, SpringP2);



PlotSnakes( P0, P1, P2, P1M, P2M, SpringP1, SpringP2 ,V, borderClass, Eext, Out, segLimits, true);



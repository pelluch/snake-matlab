% Ejemplo de snake 2d leyendo un archivo VTK con funcion ReadData3D y luego
% ajustando parametros de Snake GVF
% Ignacio Prieto 19/6/14

clc;
clear all;
close all;

%Read VTK 3D
for ind = 3:2:11
[V, info] = ReadData3D(strcat('C:/DocsMaracuya/Datos/Snakes/SAguirre/_OriginalReductor_',int2str(ind),'.vtk'));

% for frameIndex = 3:3:18;
for frameIndex = 14;

%Select mid ventricle slice on short axis
Slice = V(:,:,(19-frameIndex));
Slice = fliplr(Slice);
Slice = rot90(Slice);

% Slice = V(:,:,8);

% Generate circle


%center = [ info.Dimensions[1]/2,info.Dimensions[2]/2 ];
center = [ 160, 130 ] ;
radius = 10;
numPoints = 10;

for i=1:numPoints
    % Angles on negative to run circle clockwise
    circ(i,1) = center(1) + radius * cos(- 2 * pi * i / numPoints); 
    circ(i,2) = center(2) + radius * sin(- 2 * pi * i / numPoints);
end

Options.Verbose = true; % : If true show important images, default false
Options.nPoints = 100; %: Number of contour points, default 100
Options.Gamma = 1; % Time step, default 1
Options.Iterations = 200 %  Number of iterations, default 100
%
% options (Image Edge Energy / Image force))
Options.Sigma1 = 3; %  Sigma used to calculate image derivatives, default 10
Options.Wline = .02; %  Attraction to lines, if negative to black lines otherwise white
                     %  lines , default 0.04
Options.Wedge = 2; %  Attraction to edges, default 2.0
Options.Wterm = .01; % Attraction to terminations of lines (end points) and
                     % corners, default 0.01
Options.Sigma2 = 2; % Sigma used to calculate the gradient of the edge energy
                     % image (which gives the image force), default 20
 
% options (Gradient Vector Flow)
Options.Mu = .002; %  Trade off between real edge vectors, and noise vectors,
                 % default 0.2. (Warning setting this to high >0.5 gives
                 % an instable Vector Flow)
Options.GIterations = 0; % Number of GVF iterations, default 0
Options.Sigma3 = 1.; %  Sigma used to calculate the laplacian in GVF, default 1.0
%
% options (Snake)
Options.Alpha = .02; % Membrame energy  (first order), default 0.2
Options.Beta = .2; %  Thin plate energy (second order), default 0.2
Options.Delta = .1 %  Baloon force, default 0.1
Options.Kappa = .001; %  Weight of external image force, default 2    

opengl software;
% folder = 'C:/DocsMaracuya/Dropbox/CIB/Presentaciones/Viernes/V20.6.14/Fotos/Snake2D';
folder = 'C:/DocsMaracuya/Datos/Snakes/Snake2D';
% Snake2DWrite(Slice,circ,Options, frameIndex, ind, folder);
Snake2D(Slice,circ,Options);
end
end
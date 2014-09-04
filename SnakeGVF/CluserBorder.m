% Cluster border
% Classify pixels in the border of first snake using k-means cluster
% jiprieto 27/6/14


clc;
close all;


% clear all;
% 
% %% Read VTK 3D
% for ind = 1:12
% % for ind = 1:24
% [V(:,:,:,ind), info] = ReadData3D(strcat('C:/DocsMaracuya/Datos/Snakes/SAguirre/_OriginalReductor_',int2str(ind-1),'.vtk'));
% % [V(:,:,:,ind), info] = ReadData3D(strcat('C:/epxImages/SAguirre/KarisEjeCorto/_OriginalReductor_',int2str(ind-1),'.vtk'));
% end
% 


frameNumber = 11;

StDevTime = std(V,0,4);

MIPStDev = max(StDevTime,[],3);

% figure;
% imshow(MIPStDev,[]);
% title('MIP Std Dev');


%# Create the gaussian filter with hsize = [5 5] and sigma = 2
kernel = fspecial('gaussian',[10 10],5);
%# Filter it
filteredMip = imfilter(MIPStDev,kernel,'same');

filteredMipGray = mat2gray(filteredMip);
thLevel = multithresh(filteredMipGray, 4)
otsuMask = im2bw (filteredMipGray, thLevel(3));

se5 = strel('disk', 5);
se15 = strel('disk', 15);

closeMask = otsuMask;
closeMask = imerode(closeMask, se5);
closeMask = imdilate(closeMask, se15);


maskedImage = MIPStDev .* closeMask;

% 
% figure;
% imshow(closeMask,[]);
% title('Closed Otsu Mask');
% 
% figure;
% imshow(filteredMip,[]);
% title('Filtered Mip');
% 
% figure;
% imshow(filteredMipGray,[]);
% title('Filtered Mip Gray');
% 
% 
% figure;
% imshow(maskedImage,[]);
% title('Masked Image');
% imshow(filteredOtsu,[]);
%% Generate circle
% 
% forward = true;
forward = false;

if forward
%     center = [ 90, 166 ] ;
    center = [ 160, 130 ] ;
    radius = 6;
else

    center = [ 140, 160 ] ;
    radius = 16;

end

numPoints = 100;

for i=1:numPoints
    % Angles on negative to run circle clockwise
    circ(i,1) = center(1) + radius * cos(- 2 * pi * i / numPoints); 
    circ(i,2) = center(2) + radius * sin(- 2 * pi * i / numPoints);
end


startS=4;
endS=14;

initialP = zeros(numPoints, 2, size(V,3));
if forward
    initialSlice=startS;
    lastSlice=endS;
    initialP(:,:,initialSlice-1) = circ(:,:);
    P1 = initialP(:,:,initialSlice-1);
else
    initialSlice=endS;
    lastSlice=startS;
    initialP(:,:,initialSlice+1) = circ(:,:);
    P1 = initialP(:,:,initialSlice+1);
end
    
%positive if forward and negative if backward
% for sliceIndex = initialSlice:+1:lastSlice
for sliceIndex = initialSlice:-1:lastSlice


%Select mid ventricle slice on short axis
% Slice = V(:,:,(19-frameIndex));
Slice = V(:,:,(sliceIndex),frameNumber);
% Slice = fliplr(Slice);
% Slice = rot90(Slice);
% Slice = Slice .* closeMask;

% Slice = V(:,:,8);


Options.Verbose = true; % : If true show important images, default false
Options.nPoints = 100; %: Number of contour points, default 100
Options.Gamma = .5; % Time step, default 1
Options.Iterations = 500 %  Number of iterations, default 100
%
% options (Image Edge Energy / Image force))
Options.Sigma1 = 3; %  Sigma used to calculate image derivatives, default 10
% Options.Wline = .02; %  Attraction to lines, if negative to black lines otherwise white
Options.Wline = 0; %  Attraction to lines, if negative to black lines otherwise white
                     %  lines , default 0.04
Options.Wedge = 2; %  Attraction to edges, default 2.0
% Options.Wterm = .01; % Attraction to terminations of lines (end points) and
Options.Wterm = 0; % Attraction to terminations of lines (end points) and
                     % corners, default 0.01
Options.Sigma2 = 2; % Sigma used to calculate the gradient of the edge energy
                     % image (which gives the image force), default 20
 
% options (Gradient Vector Flow)
Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors,
                 % default 0.2. (Warning setting this to high >0.5 gives
                 % an instable Vector Flow)
Options.GIterations = 1; % Number of GVF iterations, default 0
Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
%
% options (Snake)
Options.Alpha = .1; % Membrame energy  (first order), default 0.2
Options.Beta = 1; %  Thin plate energy (second order), default 0.2
Options.Delta = +0.005 %  Baloon force, default 0.1
Options.Kappa = .001; %  Weight of external image force, default 2    

opengl software;
% [initialP(:,:,frameIndex),initialJ] = Snake2DMyo(Slice,initialP(:,:,frameIndex-1),Options);
I = Slice;
if forward
    P1 = initialP(:,:,sliceIndex-1);
else
     P1 = initialP(:,:,sliceIndex+1);
end

P10 = P1;

% This function SNAKE implements the basic snake segmentation. A snake is an 
% active (moving) contour, in which the points are attracted by edges and
% other boundaries. To keep the contour smooth, an membrame and thin plate
% energy is used as regularization.
%
% [O,J]=Snake2D(I,P,Options)
%  
% inputs,
%   I : An Image of type double preferable ranged [0..1]
%   P : List with coordinates descriping the rough contour N x 2
%   Options : A struct with all snake options
%   
% outputs,
%   O : List with coordinates of the final contour M x 2
%   J : Binary image with the segmented region
%
% options (general),
%  Option.Verbose : If true show important images, default false
%  Options.nPoints : Number of contour points, default 100
%  Options.Gamma : Time step, default 1
%  Options.Iterations : Number of iterations, default 100
%
% options (Image Edge Energy / Image force))
%  Options.Sigma1 : Sigma used to calculate image derivatives, default 10
%  Options.Wline : Attraction to lines, if negative to black lines otherwise white
%                    lines , default 0.04
%  Options.Wedge : Attraction to edges, default 2.0
%  Options.Wterm : Attraction to terminations of lines (end points) and
%                    corners, default 0.01
%  Options.Sigma2 : Sigma used to calculate the gradient of the edge energy
%                    image (which gives the image force), default 20
%
% options (Gradient Vector Flow)
%  Options.Mu : Trade of between real edge vectors, and noise vectors,
%                default 0.2. (Warning setting this to high >0.5 gives
%                an instable Vector Flow)
%  Options.GIterations : Number of GVF iterations, default 0
%  Options.Sigma3 : Sigma used to calculate the laplacian in GVF, default 1.0
%
% options (Snake)
%  Options.Alpha : Membrame energy  (first order), default 0.2
%  Options.Beta : Thin plate energy (second order), default 0.2
%  Options.Delta : Baloon force, default 0.1
%  Options.Kappa : Weight of external image force, default 2
%
%
% Literature:
%   - Michael Kass, Andrew Witkin and Demetri TerzoPoulos "Snakes : Active
%       Contour Models", 1987
%   - Jim Ivins amd John Porrill, "Everything you always wanted to know
%       about snakes (but wer afraid to ask)
%   - Chenyang Xu and Jerry L. Prince, "Gradient Vector Flow: A New
%       external force for Snakes
%
% Example, Basic:
%
%  % Read an image
%   I = imread('testimage.png');
%  % Convert the image to double data type
%   I = im2double(I); 
%  % Show the image and select some points with the mouse (at least 4)
%  %figure, imshow(I); [y,x] = getpts;
%   y=[182 233 251 205 169];
%   x=[163 166 207 248 210];
%  % Make an array with the clicked coordinates
%   P=[x(:) y(:)];
%  % Start Snake Process
%   Options=struct;
%   Options.Verbose=true;
%   Options.Iterations=300;
%   [O,J]=Snake2D(I,P,Options);
%  % Show the result
%   Irgb(:,:,1)=I;
%   Irgb(:,:,2)=I;
%   Irgb(:,:,3)=J;
%   figure, imshow(Irgb,[]); 
%   hold on; plot([O(:,2);O(1,2)],[O(:,1);O(1,1)]);
%  
% Example, GVF:
%   I=im2double(imread('testimage2.png'));
%   x=[96 51 98 202 272 280 182];
%   y=[63 147 242 262 211 97 59];
%   P=[x(:) y(:)];
%   Options=struct;
%   Options.Verbose=true;
%   Options.Iterations=400;
%   Options.Wedge=2;
%   Options.Wline=0;
%   Options.Wterm=0;
%   Options.Kappa=4;
%   Options.Sigma1=8;
%   Options.Sigma2=8;
%   Options.Alpha=0.1;
%   Options.Beta=0.1;
%   Options.Mu=0.2;
%   Options.Delta=-0.1;
%   Options.GIterations=600;
%   [O,J]=Snake2D(I,P,Options);
%   
% Function is written by D.Kroon University of Twente (July 2010)

% Process inputs
defaultoptions=struct('Verbose',false,'nPoints',100,'Wline',0.04,'Wedge',2,'Wterm',0.01,'Sigma1',10,'Sigma2',20,'Alpha',0.2,'Beta',0.2,'Delta',0.1,'Gamma',1,'Kappa',2,'Iterations',100,'GIterations',0,'Mu',0.2,'Sigma3',1);
if(~exist('Options','var')), 
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))), 
        warning('snake:unknownoption','unknown options found');
    end
end

% Convert input to double
I = double(I);

% If color image convert to grayscale
if(size(I,3)==3), I=rgb2gray(I); end

% The contour must always be clockwise (because of the balloon force)
P1=MakeContourClockwise2D(P1);

% Make an uniform sampled contour description
P1=InterpolateContourPoints2D(P1,Options.nPoints);

% Transform the Image into an External Energy Image
Eext = ExternalForceImage2D(I,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);

% Make the external force (flow) field.
Fx=ImageDerivatives2D(Eext,Options.Sigma2,'x');
Fy=ImageDerivatives2D(Eext,Options.Sigma2,'y');
Fext(:,:,1)=-Fx*2*Options.Sigma2^2;
Fext(:,:,2)=-Fy*2*Options.Sigma2^2;

% Saving original Gradient
Fext0 = Fext;

% Do Gradient vector flow, optimalization
Fext=GVFOptimizeImageForces2D(Fext, Options.Mu, Options.GIterations, Options.Sigma3);

% 
% % Show the image, contour and force field
% if(Options.Verbose)
% 
%     h=figure; set(h,'render','opengl');
%     set(h,'units','normalized','outerposition',[0 0 1 1]);
%      subplot(2,2,1),
%       imshow(I,[]); 
%       hold on; plot(P(:,2),P(:,1),'b.'); 
%       title('The image with initial contour')
%      subplot(2,2,2),
%       imshow(Eext,[]); hold on; plot(P(:,2),P(:,1),'b.'); hold on;
%       title('The external energy');
%      subplot(2,2,3), 
% %      [x,y]=ndgrid(1:10:size(Fext,1),1:10:size(Fext,2)); 
%      [x,y]=ndgrid(1:size(Fext,1),1:size(Fext,2));
% %       imshow(I), hold on; quiver(y,x,Fext0(1:10:end,1:10:end,2),Fext0(1:10:end,1:10:end,1),'b'), hold on;plot(P(:,2),P(:,1),'b.'); 
% %       hold on; quiver(y,x,Fext(1:10:end,1:10:end,2),Fext(1:10:end,1:10:end,1),'r');
%       imshow(I,[]), hold on; 
% %       quiver(y,x,Fext0(1:end,1:end,2),Fext0(1:end,1:end,1),'b'), hold on;
%       plot(P(:,2),P(:,1),'b.'); 
%       hold on; quiver(y,x,Fext(1:end,1:end,2),Fext(1:end,1:end,1),'r');
%       title('The external force field ')
% 
% %      subplot(2,2,3), 
% %       [x,y]=ndgrid(1:3:size(Fext,1),1:3:size(Fext,2));
% %       imshow(I,[]), hold on; quiver(y,x,Fext(1:3:end,1:3:end,2),Fext(1:3:end,1:3:end,1)), hold on;plot(P(:,2),P(:,1),'b.'); 
% %       title('The external force field ')
%      subplot(2,2,4), 
%       imshow(I,[]), hold on; plot(P(:,2),P(:,1),'b.'); 
%       title('Snake movement ')
%     drawnow
% end


% Make the interal force matrix, which constrains the moving points to a
% smooth contour
S=SnakeInternalForceMatrix2D(Options.nPoints,Options.Alpha,Options.Beta,Options.Gamma);
h=[];
for i=1:Options.Iterations
    P1=SnakeMoveIteration2D(S,P1,Fext,Options.Gamma,Options.Kappa,Options.Delta);

%     % Show current contour
%     if(Options.Verbose)
%         if(ishandle(h)), delete(h), end
% %         h=plot(P(:,2),P(:,1),'r.');
%         c=i/Options.Iterations;
% %         plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'-','Color',[c 1-c 0]);  drawnow
%     end
end
%         h=plot(P1(:,2),P1(:,1),'y.');
% if(nargout>1)
%      J=DrawSegmentedArea2D(P1,size(I));
% end


initialP(:,:,sliceIndex) = P1;

end



%% Second snake - myocardium


for sliceIndex = initialSlice:-1:lastSlice
    
% Calculate the distance from inner snake
Pcircle = floor(initialP(:,:,sliceIndex));
curve = zeros(size(V(:,:,sliceIndex,frameNumber)));
circleCell = num2cell(Pcircle,1);

curve(sub2ind(size(curve),circleCell{:})) = 1;
solidCircle = zeros(size(curve));
for i=1:size(solidCircle,1)
    for j=1:size(solidCircle,2)
        solidCircle(i,j) = inpolygon(i,j,Pcircle(:,1),Pcircle(:,2));
    end
end
solidCircle = logical(solidCircle);

dista(:,:,sliceIndex) = bwdist(curve);
dista(:,:,sliceIndex) = dista(:,:,sliceIndex).*~solidCircle;

margin = 25;
% border=I;
% border=I;
% Now trying with smoothed input
kerSigma = 1;
ker = fspecial('gaussian',[3*kerSigma 3*kerSigma] ,kerSigma);
Ismooth=imfilter(V(:,:,sliceIndex,frameNumber), ker, 'replicate');
border(:,:,sliceIndex)=Ismooth;
end


% border=I;
border(dista > margin) = 0;
border(dista == 0) = 0;
% 
% figure;
% imshow(border,[]);
% title(strcat('Border image with margin of   ',int2str(margin)));



%% Clustering

load 'kmeansdata';
nrows = size(border,1);
ncols = size(border,2);
nslices = size(border,3);
borderK = reshape(border,nrows*ncols*nslices,1);
borderK0 = borderK;
borderK0(borderK==0) = [];

% for numClusters = 2:1:4
numClusters = 4;
[clusterIdx, clusterCenter]  = kmeans(borderK0,numClusters,...
    'emptyaction','drop','distance','sqEuclidean', ...
    'Replicates',3);
clusterIdx0 = zeros(nrows*ncols*nslices,1);
j=1;
for i=1:nrows*ncols*nslices
    if(borderK(i)~=0)
        clusterIdx0(i) = clusterIdx(j);
        j=j+1;
    end
end
% Now we have the image classified on clusterIdx0
borderClass = reshape(clusterIdx0, nrows, ncols, nslices);
% figure;
% imshow(borderClass, []);
% title('Classified border');
% 
% imlook3d(border);
% title(strcat('Border for cluster ',int2str(numClusters)) );
% imlook3d(borderClass);

% end

% numClusters = 3;
% [clusterIdx, clusterCenter]  = kmeans(borderK0,numClusters,...
%     'emptyaction','drop','distance','sqEuclidean', ...
%     'Replicates',3);
% clusterIdx0 = zeros(nrows*ncols*nslices,1);
% j=1;
% for i=1:nrows*ncols*nslices
%     if(borderK(i)~=0)
%         clusterIdx0(i) = clusterIdx(j);
%         j=j+1;
%     end
% end
% % Now we have the image classified on clusterIdx0
% borderClass = reshape(clusterIdx0, nrows, ncols, nslices);
% % figure;
% % imshow(borderClass, []);
% % title('Classified border');
% 
% imlook3d(border);
% imlook3d(borderClass);

%% Second Snake



Options.Verbose = true; % : If true show important images, default false
Options.nPoints = 100; %: Number of contour points, default 100
Options.Gamma = .25; % Time step, default 1
Options.Iterations = 200 %  Number of iterations, default 100
%
% options (Image Edge Energy / Image force))
Options.Sigma1 = 1; %  Sigma used to calculate image derivatives, default 10
% Options.Wline = .02; %  Attraction to lines, if negative to black lines otherwise white
Options.Wline = 0; %  Attraction to lines, if negative to black lines otherwise white
                     %  lines , default 0.04
Options.Wedge = 2; %  Attraction to edges, default 2.0
% Options.Wterm = .01; % Attraction to terminations of lines (end points) and
Options.Wterm = 0; % Attraction to terminations of lines (end points) and
                     % corners, default 0.01
Options.Sigma2 = 1; % Sigma used to calculate the gradient of the edge energy
                     % image (which gives the image force), default 20
 
% options (Gradient Vector Flow)
Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors,
                 % default 0.2. (Warning setting this to high >0.5 gives
                 % an instable Vector Flow)
Options.GIterations = 0; % Number of GVF iterations, default 0
Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
%
% options (Snake)
Options.Alpha = .4; % Membrame energy  (first order), default 0.2
Options.Beta = 4; %  Thin plate energy (second order), default 0.2
Options.Delta = 0.1 %  Baloon force, default 0.1
Options.Kappa = .3; %  Weight of external image force, default 2    
Options.Width = 12; %  Distance between two snakes in pixels. Diastole
% Options.Width = 15; %  Distance between two snakes in pixels.
Options.WidthWeight = .04; %  Weight of distance between snakes, default .04

P20=initialP;
P2=initialP;

for sliceIndex = initialSlice:-1:lastSlice



% Transform the Image into an External Energy Image
% Eext = ExternalForceImage2D(I,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);
% Eext = ExternalForceImage2DMyocardium(I,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1,dista,margin);
Eext = ExternalForceImage2D(borderClass(:,:,sliceIndex),Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);


% Make the external force (flow) field.
Fx=ImageDerivatives2D(Eext,Options.Sigma2,'x');
Fy=ImageDerivatives2D(Eext,Options.Sigma2,'y');
Fext(:,:,1)=-Fx*2*Options.Sigma2^2;
Fext(:,:,2)=-Fy*2*Options.Sigma2^2;

% Saving original Gradient
Fext0 = Fext;

% Do Gradient vector flow, optimalization
Fext=GVFOptimizeImageForces2D(Fext, Options.Mu, Options.GIterations, Options.Sigma3);


% Show the image, contour and force field
if(Options.Verbose)

    
%     h2=figure; set(h2,'render','opengl');
%     % Position, [x y width height]
%     set(h2,'units','normalized','outerposition',[1. 0. 1 1]);
%       subplot(2,2,1),
%       imshow(Fx,[]);
%       hold on; plot(P(:,2),P(:,1),'b.'); 
%       title('Fx = d/dx(Ext)')
%      subplot(2,2,2),
%        imshow(Fy,[]); 
%        hold on; plot(P(:,2),P(:,1),'b.');
%       title('Fy = d/dy(Ext)');
%      subplot(2,2,3), 
% %      [x,y]=ndgrid(1:10:size(Fext,1),1:10:size(Fext,2)); 
%      [x,y]=ndgrid(1:size(Fext,1),1:size(Fext,2));
% %       imshow(I), hold on; quiver(y,x,Fext0(1:10:end,1:10:end,2),Fext0(1:10:end,1:10:end,1),'b'), hold on;plot(P(:,2),P(:,1),'b.'); 
% %       hold on; quiver(y,x,Fext(1:10:end,1:10:end,2),Fext(1:10:end,1:10:end,1),'r');
%       imshow(I), hold on; quiver(y,x,Fext0(1:end,1:end,2),Fext0(1:end,1:end,1),'b'), hold on;plot(P(:,2),P(:,1),'b.'); 
%       hold on; quiver(y,x,Fext(1:end,1:end,2),Fext(1:end,1:end,1),'r');
%       title('The external force field ')
% %      subplot(2,2,3), 
% %       [x,y]=ndgrid(1:10:size(Fext,1),1:10:size(Fext,2));
% %       imshow(I), hold on; quiver(y,x,Fext(1:10:end,1:10:end,2),Fext(1:10:end,1:10:end,1)), hold on;plot(P(:,2),P(:,1),'b.'); 
% %       title('The external force field ')
% %      subplot(2,2,4), 
% %       imshow(I), hold on; plot(P(:,2),P(:,1),'b.'); 
% %       title('Snake movement ')
%     drawnow
%  
%     h=figure; set(h,'render','opengl');
%     set(h,'units','normalized','outerposition',[0 0 1 1]);
%      subplot(2,2,1),
%       imshow(borderClass,[]); 
%       hold on; plot(P10(:,2),P10(:,1),'b-'); 
%       title('The image with initial contour')
%      subplot(2,2,2),
%       imshow(Eext,[]); hold on; plot(P20(:,2),P20(:,1),'r-'); hold on;
%       title('The external energy');
%      subplot(2,2,3), 
% %      [x,y]=ndgrid(1:10:size(Fext,1),1:10:size(Fext,2)); 
%      [x,y]=ndgrid(1:size(Fext,1),1:size(Fext,2));
%       imshow(I,[]), hold on; 
% %       quiver(y,x,Fext0(1:end,1:end,2),Fext0(1:end,1:end,1),'b'), hold on;
%       plot(P20(:,2),P20(:,1),'r-'); 
%       hold on; quiver(y,x,Fext(1:end,1:end,2),Fext(1:end,1:end,1),'r');
%       title('The external force field ')

%      subplot(2,2,3), 
%       [x,y]=ndgrid(1:3:size(Fext,1),1:3:size(Fext,2));
%       imshow(I,[]), hold on; quiver(y,x,Fext(1:3:end,1:3:end,2),Fext(1:3:end,1:3:end,1)), hold on;plot(P(:,2),P(:,1),'b.'); 
%       title('The external force field ')
%      subplot(2,2,4), 
%       imshow(I,[]), hold on; 
%       plot(P10(:,2),P10(:,1),'b-'); hold on;
%       plot(P20(:,2),P20(:,1),'r-'); 
%       title('Snake movement ')
%     drawnow
end


% Make the interal force matrix, which constrains the moving points to a
% smooth contour
S=SnakeInternalForceMatrix2D(Options.nPoints,Options.Alpha,Options.Beta,Options.Gamma);
h=[];
for i=1:Options.Iterations
    P2(:,:,sliceIndex) =SnakeMoveIteration2DMyo(S,P2(:,:,sliceIndex)...
        ,initialP(:,:,sliceIndex),Fext,Options.Gamma,Options.Kappa,...
        Options.Delta, Options.Width,Options.WidthWeight);

    % Show current contour
%     if(Options.Verbose)
%         if(ishandle(h)), delete(h), end
%         h=plot(P2(:,2),P2(:,1),'r.');
%         c=i/Options.Iterations;
%         plot([P2(:,2);P2(1,2)],[P2(:,1);P2(1,1)],'-','Color',[c 1-c 0]);  drawnow
%     end
end
%         h=plot(P2(:,2),P2(:,1),'g-');
% if(nargout>1)
%      J=DrawSegmentedArea2D(P2,size(I));
% end
end


    h=figure; set(h,'render','opengl');
    set(h,'units','normalized','outerposition',[0 0 1 1]);
     subplot(2,2,1),
      i=4; imshow(V(:,:,i,frameNumber),[]);hold on;
      plot(P2(:,2,i),P2(:,1,i),'g-');hold on;
      plot(initialP(:,2,i),initialP(:,1,i),'r-');
      title(strcat('Phase: ',int2str(frameNumber),', Slice: ',int2str(i)));
     subplot(2,2,3),
     imshow(borderClass(:,:,i),[]);hold on;
     plot(P2(:,2,i),P2(:,1,i),'g-');
     
     subplot(2,2,2),
      i=7; imshow(V(:,:,i,frameNumber),[]);hold on;
      plot(P2(:,2,i),P2(:,1,i),'g-');hold on;
      plot(initialP(:,2,i),initialP(:,1,i),'r-');
      title(strcat('Phase: ',int2str(frameNumber),', Slice: ',int2str(i)));
     subplot(2,2,4),
     imshow(borderClass(:,:,i),[]);hold on;
     plot(P2(:,2,i),P2(:,1,i),'g-');
     
    drawnow
    
    
    h=figure; set(h,'render','opengl');
    set(h,'units','normalized','outerposition',[1 0 1 1]);
     subplot(2,2,1),
      i=10; imshow(V(:,:,i,frameNumber),[]);hold on;
      plot(P2(:,2,i),P2(:,1,i),'g-');hold on;
      plot(initialP(:,2,i),initialP(:,1,i),'r-');
      title(strcat('Phase: ',int2str(frameNumber),', Slice: ',int2str(i)));
     subplot(2,2,3),
     imshow(borderClass(:,:,i),[]);hold on;
     plot(P2(:,2,i),P2(:,1,i),'g-');
     
     subplot(2,2,2),
      i=13; imshow(V(:,:,i,frameNumber),[]);hold on;
      plot(P2(:,2,i),P2(:,1,i),'g-');hold on;
      plot(initialP(:,2,i),initialP(:,1,i),'r-');
      title(strcat('Phase: ',int2str(frameNumber),', Slice: ',int2str(i)));
     subplot(2,2,4),
     imshow(borderClass(:,:,i),[]);hold on;
     plot(P2(:,2,i),P2(:,1,i),'g-');
     
    drawnow
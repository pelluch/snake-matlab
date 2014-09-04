% Now trying to segment the middle slice of a short axis image, and use the
% final contour as the initial seg for the next slice
% Ignacio Prieto 20/6/14

clc;
close all;

%Read VTK 3D
clear all;
for ind = 1:12
[V(:,:,:,ind), info] = ReadData3D(strcat('C:/DocsMaracuya/Datos/Snakes/SAguirre/_OriginalReductor_',int2str(ind-1),'.vtk'));
end


StDevTime = std(V,0,4);

MIPStDev = max(StDevTime,[],3);

figure;
imshow(MIPStDev,[]);
title('MIP Std Dev');


%# Create the gaussian filter with hsize = [5 5] and sigma = 2
kernel = fspecial('gaussian',[10 10],5);
%# Filter it
filteredMip = imfilter(MIPStDev,kernel,'same');

filteredMipGray = mat2gray(filteredMip);

thLevel = multithresh(filteredMipGray, 4)



% maxi = max(max(MIPStDev));
% otsu4Mask = zeros(size(MIPStDev));
% otsu4Mask(:,:) = true;

% for i = 1:(size(filteredMip,1))
%     for j = 1:(size(filteredMip,1))
%         if filteredMip(i,j) > thLevel(3)
%             otsuMask(i,j) = true;
%         else
%             otsuMask(i,j) = false;
%         end
%     end
% end
otsuMask = im2bw (filteredMipGray, thLevel(3));

se5 = strel('disk', 5);
se15 = strel('disk', 15);

closeMask = otsuMask;
closeMask = imerode(closeMask, se5);
closeMask = imdilate(closeMask, se15);





% otsu4Mask = im2bw(MIPStDev,  thLevel(3));
% otsu4Mask = im2bw(MIPStDev, (thLevel(4)/ maxi));


maskedImage = MIPStDev .* closeMask;
% maskedImage = MIPStDev .* otsuMask;
% maskedImage = MIPStDev .* filteredMip;

figure;
imshow(otsuMask,[]);
title('Otsu Mask');

figure;
imshow(closeMask,[]);
title('Closed Otsu Mask');

figure;
imshow(filteredMip,[]);
title('Filtered Mip');

figure;
imshow(filteredMipGray,[]);
title('Filtered Mip Gray');


figure;
imshow(maskedImage,[]);
title('Masked Image');
% imshow(filteredOtsu,[]);
%%
% Generate circle
center = [ 160, 130 ] ;
radius = 10;
numPoints = 100;

for i=1:numPoints
    % Angles on negative to run circle clockwise
    circ(i,1) = center(1) + radius * cos(- 2 * pi * i / numPoints); 
    circ(i,2) = center(2) + radius * sin(- 2 * pi * i / numPoints);
end

initialFrame = 8;
initialP = zeros(numPoints, 2, size(V,3));
initialP(:,:,initialFrame-1) = circ(:,:);

for frameIndex = initialFrame:17;
% for frameIndex = 14;

%Select mid ventricle slice on short axis
Slice = V(:,:,(19-frameIndex));
Slice = fliplr(Slice);
Slice = rot90(Slice);
Slice = Slice .* closeMask;

% Slice = V(:,:,8);


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
Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors,
                 % default 0.2. (Warning setting this to high >0.5 gives
                 % an instable Vector Flow)
Options.GIterations = 1; % Number of GVF iterations, default 0
Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
%
% options (Snake)
Options.Alpha = .02; % Membrame energy  (first order), default 0.2
Options.Beta = 2; %  Thin plate energy (second order), default 0.2
Options.Delta = .1 %  Baloon force, default 0.1
Options.Kappa = .001; %  Weight of external image force, default 2    

opengl software;
% folder = 'C:/DocsMaracuya/Dropbox/CIB/Presentaciones/Viernes/V20.6.14/Fotos/Snake2D';
folder = 'C:/DocsMaracuya/Datos/Snakes/Snake2D';
% Snake2DWrite(Slice,circ,Options, frameIndex, ind, folder);
% [previous, initialJ] = Snake2D(Slice,initialP(:,:,frameIndex-1),Options);
% initialP(:,:,frameIndex) = previous;
[initialP(:,:,frameIndex),initialJ] = Snake2D(Slice,initialP(:,:,frameIndex-1),Options);
end
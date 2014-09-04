function K=InterpolateContourPoints2DConvHull(P,xsize,ysize, minSize)
% This function resamples a few points describing a countour , to a smooth
% contour of uniform sampled points.
%
% K=InterpolateContourPoints(P,nPoints)
%
% input,
%  P : Inpute Contour, size N x 2  (with N>=4)
%  nPoints : Number of Contour points as output
%  xsize: n rows
%  ysize: n cols
%  minSize: minimum area where the snake is re-initialized for the next
%    slice/frame
%
% output,
%  K : Uniform sampled Contour points, size nPoints x 2
%
% Function is written by D.Kroon University of Twente (July 2010)
%
% example,
%  % Show an image
%   figure, imshow(imread('moon.tif'));
%  % Select some points with the mouse
%   [y,x] = getpts;
%  % Make an array with the clicked coordinates
%   P=[x(:) y(:)];
%  % Interpolate inbetween the points
%   Pnew=InterpolateContourPoints2D(P,100)
%  % Show the result
%   hold on; plot(P(:,2),P(:,1),'b*');
%   plot(Pnew(:,2),Pnew(:,1),'r.');
%
% Function is written by D.Kroon University of Twente (July 2010)
%
% % Interpolate points inbetween
% O(:,1)=interp([P(end-3:end,1);P(:,1);P(:,1);P(1:4,1)],10);
% O(:,2)=interp([P(end-3:end,2);P(:,2);P(:,2);P(1:4,2)],10);
% O=O(41:end-39,:);
%
% % Calculate distance between points
% dis=[0;cumsum(sqrt(sum((O(2:end,:)-O(1:end-1,:)).^2,2)))];
%
% % Resample to make uniform points
% K(:,1) = interp1(dis,O(:,1),linspace(0,dis(end),nPoints*2));
% K(:,2) = interp1(dis,O(:,2),linspace(0,dis(end),nPoints*2));
% K=K(round(end/4):round(end/4)+nPoints-1,:);
%
%
solidPolygon = poly2mask(P(:,2),P(:,1), xsize, ysize);

% if the circle collapsed, create a new circle

if(sum(sum(solidPolygon)) < minSize)
    
    %     [xpoints, ypoints] = find( solidPolygon == 1);
    xpoints = P(:,2);
    ypoints = P(:,1);
    center = [mean( xpoints) , mean( ypoints)];
    numPoints = minSize;

    radius = sqrt(minSize/pi);
    
    circ =zeros(numPoints,2);
    for i =1 : numPoints
        
        % Angles on negative to run circle clockwise
        circ(i,1) = center(1) + radius * cos(- 2 * pi * i / numPoints);
        circ(i,2) = center(2) + radius * sin(- 2 * pi * i / numPoints);
        
    end
    
    K = circ;
else
    
    bord = bwboundaries(solidPolygon);
    try
        K = cell2mat(bord(1));
    catch
        K = P;
    end
end
% if (max(max(solidPolygon)) ==0)
%     pause on;
%     pause;
% end

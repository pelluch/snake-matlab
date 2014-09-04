function [area2, x, diffRadius] = CalculateInterSnakeDistance( P2, area1 )
%CalculateInterSnakeDistance Gets the area of the inner and outer snake and
% assumes them as circles, returning the difference of the radius.
%   x = radius difference
%   BorderImage = m,n,s segmented image
%   P2 = OuterSnake, NPoints x 2


area2 = zeros(size(area1));

for i = 1:size(P2,1)
    
    try
        Pcircle = round(P2{i,1});
        solidPoly = poly2mask(Pcircle(:,2),Pcircle(:,1),max(max(Pcircle)),max(max(Pcircle)));
        area2(i) = sum(sum(solidPoly == 1));
    catch
        area2(i) = 0;
    end

end

radius1 = (area1./pi).^(.5);
radius2 = (area2./pi).^(.5);


diffRadius = radius2-radius1;

diffRadiusNon0 = diffRadius;
diffRadiusNon0(diffRadius==0) = [];

x = mean(diffRadiusNon0);


        % Area calculation for inter snake distance
    

end


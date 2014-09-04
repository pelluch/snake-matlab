function [area2, diffRadius] = CalculateInterSnakeDistancePerSlice( P2, area1 )
%CalculateInterSnakeDistance Gets the area of the inner and outer snake and
% assumes them as circles, returning the difference of the radius.
%   x = radius difference
%   BorderImage = m,n,s segmented image
%   P2 = OuterSnake, NPoints x 2


try
    Pcircle = round(P2);
    solidPoly = poly2mask(Pcircle(:,2),Pcircle(:,1),max(max(Pcircle)),max(max(Pcircle)));
    area2 = sum(sum(solidPoly == 1));
catch
    area2 = 0;
end

radius1 = (area1/pi)^(.5);
radius2 = (area2/pi)^(.5);


diffRadius = radius2-radius1;

end


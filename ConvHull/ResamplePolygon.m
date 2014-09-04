function y = ResamplePolygon(x, minDistance)
% x = nx2
% y = mx2
k_y = 1;
y = zeros(0,2);
for i=1:size(x,1)-1

    dist = sqrt( (x(i+1,1)-x(i,1))^2 + (x(i+1,2)-x(i,2))^2); 
    
    if(dist < minDistance/2 )
    elseif( dist <= minDistance )
        
        y(k_y, :) = x(i,:);
        k_y = k_y + 1;
        
    else
        % set a point in the middle and call recursively
        z   = [x(i,:); (x(i,:)+x(i+1,:))/2 ;x(i+1,:)];
        z2  = ResamplePolygon(z, minDistance);
        y   = [y; z2];
        k_y = k_y + size(z2,1);
    end
end
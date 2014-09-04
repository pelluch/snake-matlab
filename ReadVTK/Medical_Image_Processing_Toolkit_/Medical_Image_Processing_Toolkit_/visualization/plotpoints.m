function h = plotpoints(p,varargin)
% plots points using plot3
% p is a N x 3 array

    if (nargin-1)
        h=plot3(p(:,1),p(:,2),p(:,3),varargin{1:end});
    else
        h=plot3(p(:,1),p(:,2),p(:,3));
    end
end
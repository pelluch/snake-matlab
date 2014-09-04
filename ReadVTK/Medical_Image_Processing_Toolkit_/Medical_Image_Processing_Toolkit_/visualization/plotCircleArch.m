function h = plotCircleArch( circle,varargin)
%PLOTCIRCLE Plot a circle or circle arch as lines
%   circle must be a structure or an array of structures with the
%   foillowing fields:
%   circle.centre
%   circle.radius
%   circle.direction
%   circle.theta

resolution = 35;
color=[0 0 0];
linewidth=1;
linestyle = '-';
i = 1;
while (i <= size(varargin,2))
    if (strcmp( varargin{i} , 'LineWidth'))
            linewidth = varargin{i+1};
            i = i+1;
    elseif(strcmp( varargin{i} , 'Color'))
            color = varargin{i+1};
            i = i+1;
    elseif(strcmp( varargin{i} , 'LineStyle'))
            linestyle= varargin{i+1};
            i = i+1;
    end
    i = i+1;
end


ncircs = numel(circle);
for i=1:ncircs
m  = circle3DMesh(circle(i).centre,circle(i).radius,circle(i).direction,...
                            'resolution',resolution,'theta',circle(i).theta);
hold on;
h(i)=plotpoints(m.points',linestyle,'LineWidth',linewidth,'Color',color); 
hold off

end

end


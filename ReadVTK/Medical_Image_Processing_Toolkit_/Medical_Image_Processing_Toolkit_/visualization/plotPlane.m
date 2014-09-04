function handle = plotPlane(point, normalVector,varargin)
% handle = plotPlane(point, normalVector,size)
% handle = plotPlane(point, normalVector,size,varargin)
% options are
%
%   'color'     =   [r g b] 
%   'opacity'   =   int
%   'noplot'
%   'pointIsAtBorder'


color = [ 1 1 1];
opacity =1;
pointIsAtBorder = false;
s=1;
if (size(varargin,2)>0)
    i = 1;
    while (i <= size(varargin,2))
        if (strcmp( varargin{i} , 'color'))
                color=  varargin{i+1};
                i = i+1;
        elseif (strcmp( varargin{i} , 'scale'))
                s=  varargin{i+1};
                i = i+1;
        elseif(strcmp( varargin{i} , 'opacity'))
                 opacity= varargin{i+1};
                 i = i+1;
        elseif(strcmp( varargin{i} , 'noplot'))
                noplot = true;
        elseif(strcmp( varargin{i} , 'pointIsAtBorder'))
                pointIsAtBorder = true;
        end
        i = i+1;
    end
end

m = planeMesh(point, normalVector,'scale',s);

viewMesh(m,'color',color,'opacity',opacity);

 end
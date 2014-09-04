classdef AttributeType
%ATTRIBUTETYPE  Base class to handle attributed associated to meshes
%
%
%   See also IMAGETYPE, MESHTYPE, VIEWMESH.

%   Written by Alberto Gomez 2011
%   King's College London
%   OpenSource code under the BSD-2 license
%   This software is distributed with no warranties.
%   This software was designed for research purposes and not for clinical use.
    
     properties(GetAccess = 'public', SetAccess = 'public')
         attribute;
         name;
         numComp=1; % between 1 and 4
         type;
         lookup_table;
         nelements; %number of tuples
         
         attribute_array=[];
    end
    
    methods(Access = public)
        %constructor
        function obj = AttributeType(npoints)
            obj.lookup_table='default';
            if (nargin==1 && strcmp(class(npoints),'AttributeType'))
                % copy constructor
                obj.attribute_array = npoints.attribute_array;
                obj.nelements = npoints.nelements;
                obj.name = npoints.name;
                obj.attribute=npoints.attribute;
                obj.type=npoints.attribute;                
            elseif (nargin>0)
                obj.attribute_array = zeros(npoints,1);
                obj.nelements = npoints;
                obj.name = 'default';
                obj.attribute='field';
                obj.type='float';                
            end
        end
        
  
 
       
    end % methods
    
    
        
    
end



function write_vtkMesh(filename, mesh)
% Function for writing a mesh in a Visualization Toolkit (VTK) format
% 
% mesh  = write_vtkMesh(filename);
%
% examples:
%   mesh=write_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType


fid=fopen(filename,'wb');
if(fid<0)
    fprintf('could not open file %s for writing\n',filename);
    return
end



% write header
dlmwrite(filename, [ '# vtk DataFile Version 3.0' ],'delimiter','');
dlmwrite(filename, [ 'DataEntityType: Surface mesh' ],'delimiter','','-append');
dlmwrite(filename, [ 'ASCII' ],'delimiter','','-append');
dlmwrite(filename, [ 'DATASET POLYDATA' ],'delimiter','','-append');

% write points
dlmwrite(filename, [ 'POINTS ' num2str(mesh.npoints) ' float' ],'delimiter','','-append');
dlmwrite(filename, mesh.points,'delimiter',' ','-append');

%write triangles
dlmwrite(filename, [ 'POLYGONS ' num2str(mesh.ntriangles) ' ' num2str(4*mesh.ntriangles)  ],'delimiter','','-append');
dlmwrite(filename, [ ones(mesh.ntriangles,1)*3 mesh.triangles-1] ,'delimiter',' ','-append');

% write attributes

nattributes = numel(mesh.attributes);

nfielddatas = 0;
if (nattributes>0)
%if (0)
    dlmwrite(filename, [ ' ' ],'delimiter','','-append');
    dlmwrite(filename, [ 'POINT_DATA ' num2str(mesh.attributes(1).nelements) ],'delimiter','','-append');
    for j=1:nattributes
        if ~numel(lower(mesh.attributes(j).attribute))
            continue;
        end
        switch (lower(mesh.attributes(j).attribute))
            case 'scalars'
                dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' ' mesh.attributes(j).name ' ' mesh.attributes(j).type ],'delimiter','','-append');
                dlmwrite(filename, [ 'LOOKUP_TABLE ' mesh.attributes(j).lookup_table ],'delimiter','','-append');
                dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter','','-append');
            case 'field'
                nfielddatas = nfielddatas+1;
                dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' FieldData ' num2str(nfielddatas)  ],'delimiter',' ','-append');
                dlmwrite(filename, [ mesh.attributes(j).name ' ' num2str(mesh.attributes(j).numComp) ' ' num2str(mesh.attributes(j).nelements) '' mesh.attributes(j).type ],'delimiter','','-append');
                dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter',' ','-append');
             case 'vectors'
                 nfielddatas = nfielddatas+1;
                 %dlmwrite(filename, [ upper(mesh.attributes(j).attribute) ' FieldData ' num2str(nfielddatas) ],'delimiter','','-append');
                 dlmwrite(filename, [ upper(mesh.attributes(j).attribute) '  ' mesh.attributes(j).name   ' ' mesh.attributes(j).type ],'delimiter','','-append');
                 dlmwrite(filename, mesh.attributes(j).attribute_array ,'delimiter',' ','-append');
            otherwise
                % do nothing
        end


    end

end
  

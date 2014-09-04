function m =read_CHeartMesh(filenameT,filenameX)
% Function for reading a mesh in a Visualization Toolkit (VTK) format
%
% mesh  = read_vtkMesh(filename);
%
% examples:
%   mesh=read_vtkMesh('volume.vtk');
%   mesh would be of the class MeshType


fid = fopen(filenameX,'r');
M = fscanf(fid,'%i',[2 1]);
Nodes = fscanf(fid, '%f',[3 M(1)]);
fclose(fid);
Nodes = Nodes';

  fid = fopen(filenameT,'r');
MT = fscanf(fid,'%i',[2 1]);
%LETTER = fscanf(fid,'%s',[3 1]);
T = fscanf(fid, '%i',[4 MT(1)]);
fclose(fid);
 T = T';

m.Nodes = Nodes;
m.T= T;

end


function [v_out] = transformMatrix( v, mode,  m_f)
%[vout] = transformMatrix(v, mode,  m_f ) 
%
% computes the transfomation from a matrix or a dof file, taking into
% account the geometry (spacing, origin, etc.)
%   v is the volume and must be of the class ImageType
%   mode is the transformation source: 'mat' is a matrix, 'mfl' is a matrix
%   from a file and 'dof' is DOF from a file. In all cases, m_f is the
%   matrix or the filename.
%   Linear interpolation is used for the data
  method = 'nearest'; % changing this causes undesired changes in the histogram!
 % method = 'linear';
  
  v.spacing = v.spacing(:);
     
if (mode == 'mfl')
    % Read from a matrix
    [M(:,1) M(:,2) M(:,3) M(:,4)] = textread(m_f,'%f %f %f %f',4);
elseif (mode == 'dof')
    % Read from dof file
    M = dof2matrix(m_f);
elseif (mode == 'mat')
    % Use the input matrix
    M=m_f;
else
    disp('mode must be "mat", "mfl" or "dof"');
    return;
end
%show the matrix
M;

% for the interpolation
    firstVoxel_mm = v2mm([1 1 1 1]', v);
    lastVoxel_mm = v2mm([size(v.data) 1]', v);
    % center of each voxel of the input volume
    X=(firstVoxel_mm(1):v.spacing(2):lastVoxel_mm(1))';
    Y=(firstVoxel_mm(2):v.spacing(1):lastVoxel_mm(2))';
    Z=(firstVoxel_mm(3):v.spacing(3):lastVoxel_mm(3))';
   
%initialize volume   
    v_out.data = zeros(size(v.data));
    
    v_out.spacing = v.spacing; % Spacing remains the same
 % we go through the *result* image, whose size is identical to the
 % original
 newIndexX = zeros(size(v.data,1),size(v.data,2),size(v.data,3));
 newIndexY = zeros(size(v.data,1),size(v.data,2),size(v.data,3));
 newIndexZ = zeros(size(v.data,1),size(v.data,2),size(v.data,3));
 
  for i=1:size(v.data,1)
    for j=1:size(v.data,2)
            for k=1:size(v.data,3)
     
            coords = v2mm([i j k 1]',v); % voxel space to world coordinates
            txCoords=M*coords;
          
            newIndexX(i,j,k)=txCoords(1);
            newIndexY(i,j,k)=txCoords(2);
            newIndexZ(i,j,k)=txCoords(3);
                  
            end
     end
  end
 % disp('Interpolate -')

  v_out.data=interp3(X,Y,Z,v.data, newIndexX,newIndexY,newIndexZ,method);
  aux = [v.origin ; 1];
  tx_aux = M*aux;
  v_out.origin = [tx_aux(1) tx_aux(2) tx_aux(3)]';%TODO to be changed!  
  v_out.data(v_out.data~=v_out.data)=0;

    
end

function M = dof2matrix(dof_file)
%Build a matrix from the dof file
 [dummy1 dummy2 dof] = textread(dof_file,'%s %f %f',7);

T = eye(4);
T(1,4) = -dof(2);
T(2,4) = -dof(3);
T(3,4) = -dof(4);
R1 = eye(4);
R2 = R1;
R3 = R1;
dof(5:7) = dof(5:7) *pi/180;
R1(2,2) = cos(dof(5));
R1(2,3) = -sin(dof(5));
R1(3,2) = sin(dof(5));
R1(3,3) = cos(dof(5));

R2(1,1) = cos(dof(6));
R2(3,1) = -sin(dof(6));
R2(1,3) = sin(dof(6));
R2(3,3) = cos(dof(6));

R3(1,1) = cos(dof(7));
R3(1,2) = -sin(dof(7));
R3(2,1) = sin(dof(7));
R3(2,2) = cos(dof(7));

%M = R3*R2*R1*T; TODO verify which of those is the good one
M= R1*R2*R3*T
end

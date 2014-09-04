

%% Test script for Medical Image Processing Toolkit

% This is not a comprehensive test of all functionalities. I will upload
% more tests little by little.
% Feel free to make suggestions of examples.
% Queries: alberto.gomez@kcl.ac.uk

%% Basic operations
% 1. Image creation

imsiz = [20,20,20];
imspcg = [0.7,0.9,0.86];
imorig = -(imsiz-1)/2.*imspcg;
imgorient = eye(3);

img = ImageType(imsiz,imorig,imspcg,imgorient);
img.data(3:10,5:8,2:17)=rand(8,4,16)*256;

% 2. Image copy

img2 = ImageType(img);

% 3. Image IO

write_mhd('img.mhd',img); %write image
img3 = read_mhd('img.mhd'); % read image

sum(abs(img3.data(:)-img2.data(:)))

if sum(abs(img3.data(:)-img2.data(:)))==0
    disp('Both images are equal, test passed');
end

% 4. Retrieve the value of a couple of points point, using different interpolation schemes

point = [-2.3 -3 -2; -2.8 -3.2 -2]';

vNN1 = img.GetValue(point);
vNN2 = img.GetValue(point,'NN');
vlinear = img.GetValue(point,'linear');
vcubic = img.GetValue(point,'spline');


disp(sprintf('px\tpy\tpz\tdef\tNN\tlinear\tcubic'));
disp( num2str([point' vNN1' vNN2' vlinear' vcubic'])  );

% 5. Get the continuous index of a point

disp('COntinuous indices')
disp(img3.GetContinuousIndex(point))


%% Create a vector image which goes as flow

v=flow();
imv = VectorImageType(size(v),[0 0 0]',[1 1 1]',eye(3));
[gx,gy,gz]=gradient(v,imv.spacing(1),imv.spacing(2),imv.spacing(3));
imv.datax=gx.*v;
imv.datay=gy.*v;
imv.dataz=gz.*v;
imv_scalar = ImageType(imv);
imv_scalar.data=v;


write_mhd('flow_scalar.mhd',imv_scalar); %write image
write_mhd('flow_vector.mhd',imv); %write image

% resample an image

imv_scalar_downsampled= resampleImage(imv_scalar,[],'spacing',[2 2.1 2.2]','interpolation','linear');
imv_scalar_upsampled = resampleImage(imv_scalar_downsampled,imv_scalar,'interpolation','linear');
write_mhd('flow_scalar_downsampled.mhd',imv_scalar_downsampled); %write image
write_mhd('flow_scalar_upsampled.mhd',imv_scalar_upsampled); %write image

% reslice an image with an arbitrary plane
 slicing_plane_normal=[0.490443 -0.500059 -0.713727]';
 slicing_plane_point=[12.3662 24.1266 11.4671]';
imv_scalar_slice = resliceImage(imv_scalar,'plane',slicing_plane_normal, slicing_plane_point,'interpolation','linear');
    % imv_scalar_slice is a 2D image oriented accorfinly to the plane 
write_mhd('flow_scalar_slice.mhd',imv_scalar_slice); %write image









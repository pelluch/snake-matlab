function Eextern = ExternalForceImage2DMyocardium(I,Wline, Wedge, Wterm,Sigma,Distance,margin)
% Eextern = ExternalForceImage2D(I,Wline, Wedge, Wterm,Sigma)
% 
% inputs, 
%  I : The image
%  Sigma : Sigma used to calculated image derivatives 
%  Wline : Attraction to lines, if negative to black lines otherwise white
%          lines
%  Wedge : Attraction to edges
%  Wterm : Attraction to terminations of lines (end points) and corners
%
% outputs,
%  Eextern : The energy function described by the image
%
% Function is written by D.Kroon University of Twente (July 2010)

Ix=ImageDerivatives2D(I,Sigma,'x');
Iy=ImageDerivatives2D(I,Sigma,'y');
Ixx=ImageDerivatives2D(I,Sigma,'xx');
Ixy=ImageDerivatives2D(I,Sigma,'xy');
Iyy=ImageDerivatives2D(I,Sigma,'yy');

Eline = imgaussian(I,Sigma);
Eterm = (Iyy.*Ix.^2 -2*Ixy.*Ix.*Iy + Ixx.*Iy.^2)./((1+Ix.^2 + Iy.^2).^(3/2));
Eedge = sqrt(Ix.^2 + Iy.^2);

Distance = Distance - margin;
Distance = heaviside(Distance).*Distance;

figure;
imshow(Eline,[]);
title('Eline');


figure;
imshow(Eterm,[]);
title('Eterm');

figure;
imshow(Eedge,[]);
title('Eedge');

figure;
imshow(Distance,[]);
title('Distance from first snake');

Eextern= (Wline*Eline - Wedge*Eedge -Wterm * Eterm); 
figure;
imshow(Eextern,[]);
title('Eextern');

Eextern2= (Wline*Eline - Wedge*Eedge -Wterm * Eterm)./((Distance+1)); 
figure;
imshow(Eextern2,[]);
title('Eextern2');

Eextern3= (Wline*Eline - Wedge*Eedge -Wterm * Eterm).^(1./(Distance+1)); 
figure;
imshow(Eextern3,[]);
title('Eextern3');

Eextern4= (Wline*Eline - Wedge*Eedge -Wterm * Eterm)./((Distance+1)^2); 
figure;
imshow(Eextern4,[]);
title('Eextern4');

Eextern5= (Wline*Eline - Wedge*Eedge -Wterm * Eterm).^(1./(Distance+1)); 
figure;
imshow(Eextern5,[]);
title('Eextern5');

Eextern6= (Wline*Eline - Wedge*Eedge -Wterm * Eterm).*(Distance+1); 
figure;
imshow(Eextern6,[]);
title('Eextern6');

Eextern7= (Wline*Eline - Wedge*Eedge -Wterm * Eterm)+(Distance+1).^1.5; 
figure;
imshow(Eextern7,[]);
title('Eextern7');
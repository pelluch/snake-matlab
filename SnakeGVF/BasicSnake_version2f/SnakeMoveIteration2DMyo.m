function [ExtSnake,  P] = SnakeMoveIteration2DMyo(B,P,P0,Fext,gamma,kappa,delta, width,ww, r, springs)
% This function will calculate one iteration of contour Snake movement
%
% P=SnakeMoveIteration2D(B,P,P0,Fext,gamma,kappa,delta,width)
%
% inputs,
%   B : Internal force (smoothness) matrix
%   P : The contour points N x 2;
%   P0: The original contour points N x 2;
%   Fext : External vector field (from image)
%   gamma : Time step
%   kappa : External (image) field weight
%   delta : Balloon Force weight
%   width : Distance between two snakes.
%
% outputs,
%   P : The (moved) contour points N x 2;
%
% Function is written by D.Kroon University of Twente (July 2010)

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext,1));
P(:,2)=min(max(P(:,2),1),size(Fext,2));

% Get image force on the contour points
Fext1(:,1)=kappa*interp2(Fext(:,:,1),P(:,2),P(:,1));
Fext1(:,2)=kappa*interp2(Fext(:,:,2),P(:,2),P(:,1));
% Interp2, can give nan's if contour close to border
Fext1(isnan(Fext1))=0;

try
    N0=GetContourNormals2D(P0);
catch
    pause;
end


Fext0=ww*(P0+width*N0-P);

ExtSnake = P0 + width*N0;


% Calculate the baloonforce on the contour points
N=GetContourNormals2D(P);
Fext2=(delta)*N;

% Get Spring force to points defined by user
FSpring = SpringForce(P, springs, r);

% Update contour positions
ssx = gamma*P(:,1) + Fext1(:,1) + Fext2(:,1) + Fext0(:,1) + FSpring(:,1);
ssy = gamma*P(:,2) + Fext1(:,2) + Fext2(:,2) + Fext0(:,2) + FSpring(:,2);
P(:,1) = B * ssx;
P(:,2) = B * ssy;

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext,1));
P(:,2)=min(max(P(:,2),1),size(Fext,2));

    

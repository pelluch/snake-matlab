function P1Corrected = UpdateNumberOfPoints(P2, P1)
% This function will adjust the number of points of P1 to match P2
%
% P1Corrected = UpdateNumberOfPoints(P2, P1)
% 
% inputs,
%   P2 : External snake
%   P1 : Internal snake
%
% outputs,
%   P1Corrected : P1 with P2 size

P1Corrected = P2;

for i=1:size(P2,1)
    distances = ( P2(i,1)-P1(:,1) ).^2 +  ( P2(i,2)-P1(:,2) ).^2 ;
    [minDist, minIdx] = min(distances);
    
    P1Corrected(i,:) = [ P1(minIdx,1), P1(minIdx,2)];
end

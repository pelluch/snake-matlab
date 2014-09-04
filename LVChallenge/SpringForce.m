function FSpring = SpringForce(P, Springs, wS)
% Calculate spring force according to Kass 1988
%
% FSpring = SpringForce(P, Springs, wS)
% 
%           wS * (P - Springs[i])^2
% input,
%     P: m,2 Snake points
%     Springs: n,2 N Springs 
%     wS: spring force weight
% output,
%     F: m,2 Output force


if(size(Springs,1)>0)

PSpringsDist = zeros(size(P));
% Calculate spring force
for i=1:size(Springs,1)
PSpringsDist(:,1,i) =  Springs(i,1) - P(:,1);
PSpringsDist(:,2,i) =  Springs(i,2) - P(:,2);
end

PSpringsDist2 = PSpringsDist(:,1,:).^2 + PSpringsDist(:,2,:).^2;
[minDist, minIdx] = min(PSpringsDist2);
minDist = squeeze(minDist);
minIdx = squeeze(minIdx);

MinSpringDist = PSpringsDist(minIdx,:,:);
MinSpringDist2 = MinSpringDist.^2;
MinSpringDist3 = MinSpringDist2 .* sign(MinSpringDist);



% FSpring = P(minIdx,:);
FSpring = zeros(size(P));

for j = 1:size(minIdx);
    FSpring(minIdx(j),:) = wS * MinSpringDist3(j,:,j);
end
springLimit = 10;
FSpring(:,1) = max(min(FSpring(:,1),springLimit),-springLimit);
FSpring(:,2) = max(min(FSpring(:,2),springLimit),-springLimit);

else
    FSpring = zeros(size(P));
end
    


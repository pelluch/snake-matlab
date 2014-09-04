function [P0, P1, P2, borderClass, Eext2] = DoubleSnake(V, segLimits, circle, minSize, P1Options, P2Options, SpringsP1, SpringsP2)
% Generate both snakes and save the result on png files.
% Classify pixels in the border of first snake using k-means cluster, segment from diastole to systole
% input,
%    V = Data, 4D x,y,z,t
%    segLimits.startS = first Z axis slice
%    segLimits.endS = last Z axis slice
%    segLimits.startF = first time frame
%    segLimits.endF = last time frame
%    segLimits.dir = segLimits.dir: apex->base;~segLimits.dir:base->apex
%    circle = struct with initial circle parameters: x,y,r
%    minSize = minimum area of polygon where snake is reinitialized
%    P1Options = struct with snake 1 params:        
%         P1Options.Gamma = 1; % Time step, default 1
%         P1Options.Iterations = 500 %  Number of iterations, default 100
%         P1Options.Sigma1 = 3; %  Sigma used to calculate image derivatives, default 10
%         P1Options.Wline = 0; %  Attraction to lines, if negative to black lines otherwise white lines , default 0.04
%         P1Options.Wedge = 2; %  Attraction to edges, default 2.0
%         P1Options.Wterm = 0; % Attraction to terminations of lines (end points) and corners, default 0.01
%         P1Options.Sigma2 = 2; % Sigma used to calculate the gradient of the edge energy image (which gives the image force), default 20
%         % options (Gradient Vector Flow)
%         P1Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors, default 0.2. (Warning setting this to high >0.5 gives an instable Vector Flow)
%         P1Options.GIterations = 1; % Number of GVF iterations, default 0
%         P1Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
%         P1Options.Alpha = .1; % Membrame energy  (first order), default 0.2
%         P1Options.Beta = 1; %  Thin plate energy (second order), default 0.2
%         P1Options.Delta = +0.009 %  Baloon force, default 0.1
%         P1Options.Kappa = .001; %  Weight of external image force, default 2
%    P2Options = struct with snake 2 params:        
%         P2Options.nPoints = 100; %: Number of contour points, default 100
%         P2Options.Gamma = .25; % Time step, default 1
%         P2Options.Iterations = 200 %  Number of iterations, default 100
%         P2Options.Sigma1 = 1; %  Sigma used to calculate image derivatives, default 10
%         P2Options.Wline = 0; %  Attraction to lines, if negative to black lines otherwise white lines , default 0.04
%         P2Options.Wedge = 2; %  Attraction to edges, default 2.0
%         P2Options.Wterm = 0; % Attraction to terminations of lines (end points) and corners, default 0.01
%         P2Options.Sigma2 = 1; % Sigma used to calculate the gradient of the edge energy image (which gives the image force), default 20
%         P2Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors, default 0.2. (Warning setting this to high >0.5 gives an instable Vector Flow)
%         P2Options.GIterations = 0; % Number of GVF iterations, default 0
%         P2Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
%         P2Options.Alpha = .4; % Membrame energy  (first order), default 0.2
%         P2Options.Beta = 4; %  Thin plate energy (second order), default 0.2
%         P2Options.Delta = .1 %  Baloon force, default 0.1
%         P2Options.Kappa = .7; %  Weight of external image force, default 2
%         P2Options.minWidth = 10; %  Distance between two snakes in pixels.
%         P2Options.maxWidth = 16; %  Distance between two snakes in pixels.
%         P2Options.medianKernel = 4; % Median filter before clustering kernel
%         P2Options.minTopCluster = 3 % First class to be joined, up to numClust
%         P2Options.numClust = 4 % Number of kmeans clusters
%         P2Options.WidthWeight = .04; %  Weight of distance between snakes, default .04
%         



% jiprieto 1/7/14







%% Generate circle


circ =zeros(circle.numPoints,2);

for i=1:circle.numPoints
    % Angles on negative to run circle clockwise
    circ(i,1) = circle.center(1) + circle.radius * cos(- 2 * pi * i / circle.numPoints);
    circ(i,2) = circle.center(2) + circle.radius * sin(- 2 * pi * i / circle.numPoints);
end



P0 = cell(size(V,3),size(V,4));
P1 = cell(size(V,3),size(V,4));
P2 = cell(size(V,3),size(V,4));

if segLimits.dir
    initialSlice=segLimits.startS;
    lastSlice=segLimits.endS;
    if(initialSlice ==1)
        P1{initialSlice,segLimits.endF} = circ(:,:);
    else
        P1{initialSlice-1,segLimits.endF} = circ(:,:);
    end
    
    index=1; % loop orientation
else
    initialSlice=segLimits.endS;
    lastSlice=segLimits.startS;
    P1{initialSlice+1,segLimits.endF} = circ(:,:);
    index=-1; % loop orientation
end

for fNum = segLimits.endF:-1:segLimits.startF
    for sIdx = initialSlice:index:lastSlice
        
        %Select mid ventricle slice on short axis
        I = V(:,:,sIdx,fNum);
        
        if (fNum ~= segLimits.endF)
            P1{sIdx,fNum} = P1{sIdx,fNum+1};
        else
            if (segLimits.dir && sIdx ==1)
                P1{sIdx,fNum+1} = P1{sIdx,fNum};
            elseif (segLimits.dir && sIdx ~=1)
                P1{sIdx,fNum} = P1{sIdx-1,fNum};
                P1{sIdx,fNum+1} = P1{sIdx-1,fNum};
            else 
                P1{sIdx,fNum} = P1{sIdx+1,fNum};
                P1{sIdx,fNum+1} = P1{sIdx+1,fNum};
            end
        end
        
%         if(sIdx == 8 && fNum == 17 )
%         pause on;
%         pause;
%         end

        if ( max(max( isnan (P1{sIdx,fNum}))))
            pause on;
            pause;
        end
        % Make an uniform sampled contour description
        P1{sIdx,fNum} = InterpolateContourPoints2DConvHull(P1{sIdx,fNum}...
            , size(V,1),size(V,2),minSize);
        
        if ( max(max( isnan (P1{sIdx,fNum}))))
            pause on;
            pause;
        end
        
        
        if (min(min(P1{sIdx,fNum})) == 1)
            pause on;
            pause;
        end
        % Transform the Image into an External Energy Image
        Eext = ExternalForceImage2D(I,P1Options.Wline, P1Options.Wedge, P1Options.Wterm,P1Options.Sigma1);
        
        % Make the external force (flow) field.
        Fx=ImageDerivatives2D(Eext,P1Options.Sigma2,'x');
        Fy=ImageDerivatives2D(Eext,P1Options.Sigma2,'y');
        Fext(:,:,1)=-Fx*2*P1Options.Sigma2^2;
        Fext(:,:,2)=-Fy*2*P1Options.Sigma2^2;
        
        % Saving original Gradient
        % Fext0 = Fext;
        
        % Do Gradient vector flow, optimalization
        % Fext=GVFOptimizeImageForces2D(Fext, P1Options.Mu, P1Options.GIterations, P1Options.Sigma3);
        
%         spring = [85, 105; 200, 200; 100, 100];
        
        
        % Make the interal force matrix, which constrains the moving points to a
        % smooth contour
        % S=SnakeInternalForceMatrix2D(P1Options.nPoints,P1Options.Alpha,P1Options.Beta,P1Options.Gamma);
        S=SnakeInternalForceMatrix2D(size(P1{sIdx,fNum},1),P1Options.Alpha,P1Options.Beta,P1Options.Gamma);
        h=[];
%         figure;
        for i=1:P1Options.Iterations
%             P1{sIdx,fNum} = SnakeMoveIteration2D(S,P1{sIdx,fNum},Fext,P1Options.Gamma,P1Options.Kappa,P1Options.Delta);
            P1{sIdx,fNum} = SnakeMoveIteration2DInterac(S,P1{sIdx,fNum},Fext,P1Options.Gamma,P1Options.Kappa,P1Options.Delta, P1Options.r, SpringsP1{sIdx,fNum});
            %     imshow(I,[]); hold on; plot(P1{sIdx,fNum}(:,2),P1{sIdx,fNum}(:,1),'r.');
%             if(sIdx == 8 && fNum == 17 )
%                 imshow(I,[]); hold on; plot(P1{sIdx,fNum}(:,2),P1{sIdx,fNum}(:,1),'r-');
%                 pause on;
%                 pause;
%             end
%         if (min(min(P1{sIdx,fNum})) == 1)
%             pause on;
%             pause;
%         end
            
        end
        close all;
        
    end
end



%% Myocardium Clustering

area1 = zeros(size(V,3),size(V,4));
distances = zeros(size(V,3),size(V,4));

for fNum = segLimits.endF:-1:segLimits.startF
    for sIdx = initialSlice:index:lastSlice
        
        % Calculate the distance from inner snake
        Pcircle = round(P1{sIdx,fNum});
        solidPoly = poly2mask(Pcircle(:,2),Pcircle(:,1),size(V,2),size(V,1));
        dista(:,:,sIdx,fNum) = bwdist(solidPoly);
        
        % Area calculation for inter snake distance
        area1(sIdx,fNum) = sum(sum(solidPoly == 1));
        
        margin = 1.8*P2Options.maxWidth;
%         
%         % Now trying with smoothed input
%         kerSigma = 1;
%         ker = fspecial('gaussian',[3*kerSigma 3*kerSigma] ,kerSigma);
%         border(:,:,sIdx,fNum)=imfilter(V(:,:,sIdx,fNum), ker, 'replicate');

        % Now trying with median filter
        kerSigma = P2Options.medianKernel;
        border(:,:,sIdx,fNum)= medfilt2(V(:,:,sIdx,fNum),[kerSigma kerSigma]);
        

    end
end

border(dista > margin) = 0;
border(dista == 0) = 0;



load 'kmeansdata';
nrows = size(border,1);
ncols = size(border,2);
nslices = size(border,3);
nframes = size(border,4);
borderK = reshape(border,nrows*ncols*nslices*nframes,1);
borderK0 = borderK;
borderK0(borderK==0) = [];


numClusters = P2Options.numClust;
[clusterIdx, clusterCenter]  = kmeans(borderK0,numClusters,...
    'emptyaction','drop','distance','sqEuclidean');

% Changing indexes so cluster output is intensity sorted.
[clusterCenterSorted , clusterRearrangement] = sort(clusterCenter);
clusterIdxSort = clusterIdx;
for i=1:numClusters
    clusterIdxSort(clusterIdx == clusterRearrangement(i)) = i;
end
%

for i=P2Options.minTopCluster:P2Options.numClust;
    clusterIdxSort(clusterIdxSort==i) = P2Options.minTopCluster;
end
% clusterIdxSort(clusterIdxSort==3) = 2;

clusterIdx0 = zeros(nrows*ncols*nslices*nframes,1);
j=1;
for i=1:nrows*ncols*nslices*nframes
    if(borderK(i)~=0)
        %         clusterIdx0(i) = clusterIdx(j);
        clusterIdx0(i) = clusterIdxSort(j);
        j=j+1;
    end
end
% Now we have the image classified on clusterIdx0
borderClass = reshape(clusterIdx0, nrows, ncols, nslices, nframes);

% 
% borderClass2 = V;
% mask = false(size(V));%,1),size(V,2),size(V,3),size(V,4));
% 
% for i=1:numClusters-1
%     minCluster(i) = min(borderK(clusterIdx0==i));
%     maxCluster(i) = max(borderK(clusterIdx0==i));
%    
%     mask(V <= maxCluster(i)) = true;
%     mask2=mask;
%     
%     mask2(V > minCluster(i)) = true;
%     mask = bitand(mask, mask2);
%     
% %     mask = 
%     borderClass2(mask) = i;
% end
% borderClass3(V == 0) = 0;

% figure;
% imshow(borderClass, []);
% title('Classified border');

% imlook3d(mask(:,:,:,11));
% title(strcat('MaskCluster3'));%,int2str(numClusters)) );
% % imlook3d(borderClass);
% 
% % end
% 



%% Second Snake
minDistance = P2Options.minWidth;
maxDistance = P2Options.maxWidth;
distances(1:size(V,3)+1,1:size(V,4)+1) = minDistance;


for fNum = segLimits.endF:-1:segLimits.startF

    if(fNum == segLimits.endF)
        distances(:,fNum:fNum+1) = (minDistance+maxDistance)/2;
    end
    
    if(fNum ~= segLimits.endF)
        P2(:,fNum) = P2(:,fNum+1);
    else
        P2(:,fNum) = P1(:,fNum);
    end
    
    
    
    
    
    for sIdx = initialSlice:index:lastSlice
        
        % Make an uniform sampled contour description
        P2{sIdx,fNum} = InterpolateContourPoints2DConvHull(P2{sIdx,fNum}...
            , size(V,1),size(V,2), minSize);
        
        P1Corrected = UpdateNumberOfPoints(P2{sIdx,fNum},P1{sIdx,fNum});
        

        distances2=distances;
        distances2(distances2 < minDistance) = minDistance;
        distances2(distances2 > maxDistance) = maxDistance;
        
        distanceF(sIdx,fNum) = distances2(sIdx,fNum+1);
  
        if( sIdx == 1 )
            myoWidth = distances(sIdx,fNum);
        else
            distanceS(sIdx,fNum) = distances2(sIdx+((-1)*index),fNum);
            myoWidth = distanceS(sIdx, fNum); %  Same distance only per frame
        end
%         distanceFS(sIdx,fNum) = (distanceF(sIdx,fNum)+distanceS(sIdx,fNum))/2;
        
        
        
%         myoWidth = distanceFS(sIdx, fNum); %  Distance between two snakes in pixels. Diastole
        
        
        
        % Transform the Image into an External Energy Image
        Eext2(:,:,sIdx,fNum) = ExternalForceImage2D(borderClass(:,:,sIdx,fNum),P2Options.Wline, P2Options.Wedge, P2Options.Wterm,P2Options.Sigma1);
        E2 = Eext2(:,:,sIdx,fNum);
        E2(E2 < -.8) = -.8; % threshold for avoiding excesively black borders
%         Eext2(:,:,sIdx,fNum) = double(edge(borderClass(:,:,sIdx,fNum)));
        Eext2(:,:,sIdx,fNum) = E2;
        
        Fx=ImageDerivatives2D(Eext2(:,:,sIdx,fNum),P2Options.Sigma2,'x');
        Fy=ImageDerivatives2D(Eext2(:,:,sIdx,fNum),P2Options.Sigma2,'y');
        Fext(:,:,1)=-Fx*2*P2Options.Sigma2^2;
        Fext(:,:,2)=-Fy*2*P2Options.Sigma2^2;

        
        %         E2 = Eext2(:,:,sIdx,fNum);
%         E2(E2>0) = 1;
%         Eext2(:,:,sIdx,fNum) = double(edge(borderClass(:,:,sIdx,fNum)));
%         Eext2(:,:,sIdx,fNum) = double(Eext2(:,:,sIdx,fNum));
        
        % Make the external force (flow) field.
        
%         kerSigma = P2Options.Sigma2;
%         kerX = fspecial('gaussian',[3*kerSigma 1] ,kerSigma);
%         kerY = fspecial('gaussian',[1 3*kerSigma] ,kerSigma);
%         Fext(:,:,1)=-2*P2Options.Sigma2^2*imfilter(Eext2(:,:,sIdx,fNum), kerX, 'replicate');
%         Fext(:,:,2)=-2*P2Options.Sigma2^2*imfilter(Eext2(:,:,sIdx,fNum), kerY, 'replicate');
%         [Fext(:,:,1), Fext(:,:,2)] = imgradientxy(double(edge(borderClass(:,:,sIdx,fNum))));
%         =-2*P2Options.Sigma2^2*imfilter(Eext2(:,:,sIdx,fNum), kerY, 'replicate');

%         Eext2(:,:,sIdx,fNum) = Fext(:,:,1).^2 + Fext(:,:,2).^2;

%        
% 
        
        % Saving original Gradient
        % Fext0 = Fext;
        
        % Do Gradient vector flow, optimalization
        % Fext=GVFOptimizeImageForces2D(Fext, P2Options.Mu, P2Options.GIterations, P2Options.Sigma3);
        
        
        
        % Make the interal force matrix, which constrains the moving points to a
        % smooth contour
        S=SnakeInternalForceMatrix2D(size(P2{sIdx,fNum},1),P2Options.Alpha,P2Options.Beta,P2Options.Gamma);
        h=[];
        %         P20 = P1{sIdx,fNum};
        P20 = P1Corrected;
        
        for i=1:P2Options.Iterations
            [P0{sIdx,fNum}, P2{sIdx,fNum}] =SnakeMoveIteration2DMyo(S,P2{sIdx,fNum}...
                ,P20,Fext,P2Options.Gamma,P2Options.Kappa,...
                P2Options.Delta, myoWidth,P2Options.WidthWeight,...
                P1Options.r, SpringsP2{sIdx,fNum});
            
            
        end
        
        [~ , distances(sIdx,fNum)] = CalculateInterSnakeDistancePerSlice(...
            P2{sIdx,fNum},area1(sIdx,fNum));
        
        
    end
    P2(:,segLimits.endF+1)=P2(:,segLimits.endF);
end

% %% Plot Results
% close all;
% 
% 
% 
% 
% %Zoom
% x1 = OutputParam.x1;
% x2 = OutputParam.x2;
% y1 = OutputParam.y1;
% y2 = OutputParam.y2;
% P2(:,fNum+1) = P2(:,fNum);
% try
%     rmdir(OutputParam.writeFolder,'s');
% catch
% end
% pause on;
% pause(.001);
% 
% mkdir(OutputParam.writeFolder);
% %
% for fNum = segLimits.endF:-1:segLimits.startF
%     for nFig =segLimits.startS:1:segLimits.endS
%         
%         framenum = sprintf('%02d',fNum);
%         slicenum = sprintf('%02d',nFig);
%         
%         
%         h=figure('name',strcat('Frame ',int2str(fNum),', Fig',int2str(nFig)));
%         %         set(h,'render','opengl');
%         
%         %         set(h,'units','normalized','outerposition',[(floor(nFig/10)) 0 1 1]);
%         
%         %         offset = 3*mod(nFig,2);
%         offset = 0;
%         subplot(1,3,1+offset),
%         i=nFig; imshow(V(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1  x2],'Xdata',[y1  y2]);hold on;
%         P21 = plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');hold on;
%         P22 = plot(P2{i+(-1)*index,fNum}(:,2),P2{i+(-1)*index,fNum}(:,1),'y--');hold on;
%         P23 = plot(P2{i,fNum+1}(:,2),P2{i,fNum+1}(:,1),'c--');hold on;
%         P11 = plot(P1{i,fNum+1}(:,2),P1{i,fNum+1}(:,1),'b--');hold on;
%         P12 = plot(P1{i,fNum}(:,2),P1{i,fNum}(:,1),'r-');hold on;
%         P10 = plot(P0{i,fNum}(:,2),P0{i,fNum}(:,1),'w--');hold on;
%         title(['Frame: ',int2str(fNum),', Slice: ',int2str(i)]);
% 
%         
%         subplot(1,3,2+offset),
%         imshow(borderClass(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1 x2],'Xdata',[y1 y2]);hold on;
% %         imshow(borderClass2(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1 x2],'Xdata',[y1 y2]);hold on;
%         plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');
%         subplot(1,3,3+offset),
%         imshow(Eext2(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1 x2],'Xdata',[y1 y2]);hold on;
%         plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');
% %             ,['P2 S',int2str(i),', F',int2str(fNum+1)]...
% % hL = legend([P21, P22, P11, P12]...        
%         hL = legend([P21, P22, P23, P11, P12, P10]...
%             ,['P2 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
%             ,['P2 S',sprintf('%02d',i+(-1)*index),', F',sprintf('%02d',fNum)]...
%             ,['P2 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]...
%             ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]...
%             ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
%             ,['P0 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]);
%             
%         % Construct a Legend with the data from the sub-plots
%         newPosition = [0.5 0.2 0.1 0.1];
%         newUnits = 'normalized';
%         set(hL,'Position', newPosition,'Units', newUnits,'Orientation','horizontal', 'fontsize', 6);
%         
%         
%         
%         
%         export_fig(h, strcat(OutputParam.writeFolder,'/Plot_Frame',framenum,'_Slice',slicenum,'.bmp'));
%         export_fig(h, strcat(OutputParam.writeFolder,'/Plot_Slice',slicenum,'_Frame',framenum,'.bmp'));
% %         
% %         h2=figure('name',strcat('Frame ',int2str(fNum),', Fig',int2str(nFig)));
% %         imshow(V(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1  x2],'Xdata',[y1  y2]);hold on;
% %         P21 = plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');hold on;
% %         P11 = plot(P1{i,fNum+1}(:,2),P1{i,fNum+1}(:,1),'r-');hold on;
% %         for j=15:15
% %             P21M{i} = plot(P2Manual{j}(:,1),P2Manual{j}(:,2),'g--');hold on;
% %             P11M{i} = plot(P1Manual{j}(:,1),P1Manual{j}(:,2),'r--');hold on;
% %         end
% %         title(['Frame: ',int2str(fNum),', Slice: ',int2str(i)]);
% %         export_fig(h2, strcat(OutputParam.writeFolder,'/PlotManual_Slice',slicenum,'_Frame',framenum,'.bmp'));
% 
%         
% %         close all;
%         
%     end
%     %     drawnow;
% end
% distances
% distanceS
% distanceF
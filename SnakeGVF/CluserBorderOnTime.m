% Cluster border on time
% Classify pixels in the border of first snake using k-means cluster
% Then segment from diastole to systole
% jiprieto 1/7/14


clc;
close all;


clear all;
% Read VTK 3D
for ind = 1:12
% for ind = 1:24
[V(:,:,:,ind), info] = ReadData3D(strcat('C:/DocsMaracuya/Datos/Snakes/SAguirre/_OriginalReductor_',int2str(ind-1),'.vtk'));
% [V(:,:,:,ind), info] = ReadData3D(strcat('C:/epxImages/SAguirre/KarisEjeCorto/_OriginalReductor_',int2str(ind-1),'.vtk'));
end




%% Generate circle

%
% forward = true;
forward = false;

if forward
    %     center = [ 90, 166 ] ;
    center = [ 160, 130 ] ;
    radius = 6;
else
    
    center = [ 140, 160 ] ;
    radius = 16;
    
end

numPoints = 100;

for i=1:numPoints
    % Angles on negative to run circle clockwise
    circ(i,1) = center(1) + radius * cos(- 2 * pi * i / numPoints);
    circ(i,2) = center(2) + radius * sin(- 2 * pi * i / numPoints);
end


startS=4;
endS=14;
startF = 3;
endF = 11;


P0 = cell(size(V,3),size(V,4));
P1 = cell(size(V,3),size(V,4));
P2 = cell(size(V,3),size(V,4));
% P1 = zeros(numPoints, 2, size(V,3),size(V,4));
% P2 = zeros(numPoints, 2, size(V,3),size(V,4));

if forward
    initialSlice=startS;
    lastSlice=endS;
    %     P1(:,:,initialSlice-1,endF) = circ(:,:);
    P1{initialSlice-1,endF} = circ(:,:);
else
    initialSlice=endS;
    lastSlice=startS;
    %     P1(:,:,initialSlice+1,endF) = circ(:,:);
    P1{initialSlice+1,endF} = circ(:,:);
end

for fNum = endF:-1:startF
    
    
    %positive if forward and negative if backward
    % for sliceIndex = initialSlice:+1:lastSlice
    for sIdx = initialSlice:-1:lastSlice
        
        
        
        Options.Verbose = true; % : If true show important images, default false
        Options.nPoints = 500; %: Number of contour points, default 100
        Options.Gamma = 1; % Time step, default 1
        Options.Iterations = 500 %  Number of iterations, default 100
        %
        % options (Image Edge Energy / Image force))
        Options.Sigma1 = 3; %  Sigma used to calculate image derivatives, default 10
        % Options.Wline = .02; %  Attraction to lines, if negative to black lines otherwise white
        Options.Wline = 0; %  Attraction to lines, if negative to black lines otherwise white
        %  lines , default 0.04
        Options.Wedge = 2; %  Attraction to edges, default 2.0
        % Options.Wterm = .01; % Attraction to terminations of lines (end points) and
        Options.Wterm = 0; % Attraction to terminations of lines (end points) and
        % corners, default 0.01
        Options.Sigma2 = 2; % Sigma used to calculate the gradient of the edge energy
        % image (which gives the image force), default 20
        
        % options (Gradient Vector Flow)
        Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors,
        % default 0.2. (Warning setting this to high >0.5 gives
        % an instable Vector Flow)
        Options.GIterations = 1; % Number of GVF iterations, default 0
        Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
        %
        % options (Snake)
        Options.Alpha = .1; % Membrame energy  (first order), default 0.2
        Options.Beta = 1; %  Thin plate energy (second order), default 0.2
        Options.Delta = +0.009 %  Baloon force, default 0.1
        Options.Kappa = .001; %  Weight of external image force, default 2
        
        opengl software;
        
        
        %Select mid ventricle slice on short axis
        I = V(:,:,sIdx,fNum);
        
        if (fNum ~= endF)
            P1{sIdx,fNum} = P1{sIdx,fNum+1};
        else
            if forward
                P1{sIdx,fNum} = P1{sIdx-1,fNum};
                P1{sIdx,fNum+1} = P1{sIdx-1,fNum};
            else
                P1{sIdx,fNum} = P1{sIdx+1,fNum};
                P1{sIdx,fNum+1} = P1{sIdx+1,fNum};
            end
        end
        
        
        % Make an uniform sampled contour description
        P1{sIdx,fNum} = InterpolateContourPoints2DConvHull(P1{sIdx,fNum}...
            , size(V,1),size(V,2));
        
        
        % Transform the Image into an External Energy Image
        Eext = ExternalForceImage2D(I,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);
        
        % Make the external force (flow) field.
        Fx=ImageDerivatives2D(Eext,Options.Sigma2,'x');
        Fy=ImageDerivatives2D(Eext,Options.Sigma2,'y');
        Fext(:,:,1)=-Fx*2*Options.Sigma2^2;
        Fext(:,:,2)=-Fy*2*Options.Sigma2^2;
        
        % Saving original Gradient
        % Fext0 = Fext;
        
        % Do Gradient vector flow, optimalization
        % Fext=GVFOptimizeImageForces2D(Fext, Options.Mu, Options.GIterations, Options.Sigma3);
        
        
        % Make the interal force matrix, which constrains the moving points to a
        % smooth contour
        % S=SnakeInternalForceMatrix2D(Options.nPoints,Options.Alpha,Options.Beta,Options.Gamma);
        S=SnakeInternalForceMatrix2D(size(P1{sIdx,fNum},1),Options.Alpha,Options.Beta,Options.Gamma);
        h=[];
        % figure;
        for i=1:Options.Iterations
            %     P1(:,:,sIdx,fNum) = SnakeMoveIteration2D(S,P1(:,:,sIdx,fNum),Fext,Options.Gamma,Options.Kappa,Options.Delta);
            P1{sIdx,fNum} = SnakeMoveIteration2D(S,P1{sIdx,fNum},Fext,Options.Gamma,Options.Kappa,Options.Delta);
            %     imshow(I,[]); hold on; plot(P1{sIdx,fNum}(:,2),P1{sIdx,fNum}(:,1),'r.');
        end
        
    end
end



%% Myocardium Clustering

area1 = zeros(size(V,3),size(V,4));
distances = zeros(size(V,3),size(V,4));

for fNum = endF:-1:startF
    
    for sIdx = initialSlice:-1:lastSlice
        
        % Calculate the distance from inner snake
        Pcircle = round(P1{sIdx,fNum});
        solidPoly = poly2mask(Pcircle(:,2),Pcircle(:,1),size(V,2),size(V,1));
        dista(:,:,sIdx,fNum) = bwdist(solidPoly);
        
        % Area calculation for inter snake distance
        area1(sIdx,fNum) = sum(sum(solidPoly == 1));
        
        margin = 25;
%         
%         % Now trying with smoothed input
%         kerSigma = 1;
%         ker = fspecial('gaussian',[3*kerSigma 3*kerSigma] ,kerSigma);
%         border(:,:,sIdx,fNum)=imfilter(V(:,:,sIdx,fNum), ker, 'replicate');

        % Now trying with median filter
        kerSigma = 4;
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


numClusters = 4;
[clusterIdx, clusterCenter]  = kmeans(borderK0,numClusters,...
    'emptyaction','drop','distance','sqEuclidean', ...
    'Replicates',3);

% Changing indexes so cluster output is intensity sorted.
[clusterCenterSorted , clusterRearrangement] = sort(clusterCenter);
clusterIdxSort = clusterIdx;
for i=1:numClusters
    clusterIdxSort(clusterIdx == clusterRearrangement(i)) = i;
end
%
clusterIdxSort(clusterIdxSort==4) = 3;
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


borderClass2 = V;
mask = false(size(V));%,1),size(V,2),size(V,3),size(V,4));

for i=1:numClusters-1
    minCluster(i) = min(borderK(clusterIdx0==i));
    maxCluster(i) = max(borderK(clusterIdx0==i));
   
    mask(V <= maxCluster(i)) = true;
    mask2=mask;
    
    mask2(V > minCluster(i)) = true;
    mask = bitand(mask, mask2);
    
%     mask = 
    borderClass2(mask) = i;
end
% borderClass3(V == 0) = 0;

% figure;
% imshow(borderClass, []);
% title('Classified border');

imlook3d(mask(:,:,:,11));
title(strcat('MaskCluster3'));%,int2str(numClusters)) );
% imlook3d(borderClass);

% end




%% Second Snake
minDistance = 10;
maxDistance = 18;
distances(1:size(V,3)+1,1:size(V,4)+1) = minDistance;


for fNum = endF:-1:startF
    
    %     P2(:,fNum) = P1(:,fNum);
    %
    if(fNum == endF)
        distances(:,fNum:fNum+1) = 12;
        
    end
    
    if(fNum ~= endF)
        P2(:,fNum) = P2(:,fNum+1);
    else
        P2(:,fNum) = P1(:,fNum);
    end
    
    
    
    
    
    for sIdx = initialSlice:-1:lastSlice
        
        % Make an uniform sampled contour description
        P2{sIdx,fNum} = InterpolateContourPoints2DConvHull(P2{sIdx,fNum}...
            , size(V,1),size(V,2));
        
        P1Corrected = UpdateNumberOfPoints(P2{sIdx,fNum},P1{sIdx,fNum});
        
        
        
        Options.Verbose = true; % : If true show important images, default false
        Options.nPoints = 100; %: Number of contour points, default 100
        Options.Gamma = .25; % Time step, default 1
        Options.Iterations = 200 %  Number of iterations, default 100
        %
        % options (Image Edge Energy / Image force))
        Options.Sigma1 = 1; %  Sigma used to calculate image derivatives, default 10
        % Options.Wline = .02; %  Attraction to lines, if negative to black lines otherwise white
        Options.Wline = 0; %  Attraction to lines, if negative to black lines otherwise white
        %  lines , default 0.04
        Options.Wedge = 2; %  Attraction to edges, default 2.0
        % Options.Wterm = .01; % Attraction to terminations of lines (end points) and
        Options.Wterm = 0; % Attraction to terminations of lines (end points) and
        % corners, default 0.01
        Options.Sigma2 = 1; % Sigma used to calculate the gradient of the edge energy
        % image (which gives the image force), default 20
        
        % options (Gradient Vector Flow)
        Options.Mu = .00005; %  Trade off between real edge vectors, and noise vectors,
        % default 0.2. (Warning setting this to high >0.5 gives
        % an instable Vector Flow)
        Options.GIterations = 0; % Number of GVF iterations, default 0
        Options.Sigma3 = 5.; %  Sigma used to calculate the laplacian in GVF, default 1.0
        %
        % options (Snake)
        Options.Alpha = .4; % Membrame energy  (first order), default 0.2
        Options.Beta = 4; %  Thin plate energy (second order), default 0.2
        Options.Delta = .1 %  Baloon force, default 0.1
        Options.Kappa = .7; %  Weight of external image force, default 2
        %     Options.Width = distances(sIdx,fNum+1); %  Distance between two snakes in pixels. Diastole
        % Options.Width = 12; %  Distance between two snakes in pixels.
        Options.WidthWeight = .04; %  Weight of distance between snakes, default .04
        
        df = min([max([12, (distances(sIdx,fNum+1))/2]), 22]);
        ds = min([max([12, (distances(sIdx+1,fNum))/2]), 22]);
        
        
        distances2=distances;
        distances2(distances2 < minDistance) = minDistance;
        distances2(distances2 > maxDistance) = maxDistance;
        
        distanceF(sIdx,fNum) = distances2(sIdx,fNum+1);
        distanceS(sIdx,fNum) = distances2(sIdx+1,fNum);
        
        distanceFS(sIdx,fNum) = (distanceF(sIdx,fNum)+distanceS(sIdx,fNum))/2;
        
        
        %     Options.Width = min([max([8, (df+ds)/2]), 18]); %  Distance between two snakes in pixels. Diastole
        Options.Width = distanceFS(sIdx, fNum); %  Distance between two snakes in pixels. Diastole
        
        
        
        
        % Transform the Image into an External Energy Image
        % Eext = ExternalForceImage2D(I,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);
        % Eext = ExternalForceImage2DMyocardium(I,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1,dista,margin);
%         Eext2(:,:,sIdx,fNum) = ExternalForceImage2D(borderClass(:,:,sIdx,fNum),Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);
        Eext2(:,:,sIdx,fNum) = ExternalForceImage2D(borderClass2(:,:,sIdx,fNum),Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);
        
        % Make the external force (flow) field.
        Fx=ImageDerivatives2D(Eext2(:,:,sIdx,fNum),Options.Sigma2,'x');
        Fy=ImageDerivatives2D(Eext2(:,:,sIdx,fNum),Options.Sigma2,'y');
        Fext(:,:,1)=-Fx*2*Options.Sigma2^2;
        Fext(:,:,2)=-Fy*2*Options.Sigma2^2;
        
        % Saving original Gradient
        % Fext0 = Fext;
        
        % Do Gradient vector flow, optimalization
        % Fext=GVFOptimizeImageForces2D(Fext, Options.Mu, Options.GIterations, Options.Sigma3);
        
        
        
        % Make the interal force matrix, which constrains the moving points to a
        % smooth contour
        S=SnakeInternalForceMatrix2D(size(P2{sIdx,fNum},1),Options.Alpha,Options.Beta,Options.Gamma);
        h=[];
        %         P20 = P1{sIdx,fNum};
        P20 = P1Corrected;
        
        for i=1:Options.Iterations
            %     P2(:,:,sIdx,fNum) =SnakeMoveIteration2DMyo(S,P2(:,:,sIdx,fNum)...
            [P0{sIdx,fNum}, P2{sIdx,fNum}] =SnakeMoveIteration2DMyo(S,P2{sIdx,fNum}...
                ,P20,Fext,Options.Gamma,Options.Kappa,...
                Options.Delta, Options.Width,Options.WidthWeight);
            
        end
        
        [~ , distances(sIdx,fNum)] = CalculateInterSnakeDistancePerSlice(...
            P2{sIdx,fNum},area1(sIdx,fNum));
        
        
    end
    P2(:,endF+1)=P2(:,endF);
end

%% Plot Results
close all;

% Zoom

%     x1 = 1;
%     x2 = 288;
%     y1 = 1;
%     y2 = 288;

x1 = 100;
x2 = 180;
y1 = 120;
y2 = 200;
P2(:,fNum+1) = P2(:,fNum);
%
for fNum = endF:-1:startF
    for nFig =4:1:13
        
        framenum = sprintf('%02d',fNum);
        slicenum = sprintf('%02d',nFig);
        
        
        h=figure('name',strcat('Frame ',int2str(fNum),', Fig',int2str(nFig)));
        %         set(h,'render','opengl');
        
        %         set(h,'units','normalized','outerposition',[(floor(nFig/10)) 0 1 1]);
        
        %         offset = 3*mod(nFig,2);
        offset = 0;
        subplot(1,3,1+offset),
        i=nFig; imshow(V(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1  x2],'Xdata',[y1  y2]);hold on;
        P21 = plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');hold on;
        P22 = plot(P2{i+1,fNum}(:,2),P2{i+1,fNum}(:,1),'y--');hold on;
        P23 = plot(P2{i,fNum+1}(:,2),P2{i,fNum+1}(:,1),'c--');hold on;
        P11 = plot(P1{i,fNum+1}(:,2),P1{i,fNum+1}(:,1),'b--');hold on;
        P12 = plot(P1{i,fNum}(:,2),P1{i,fNum}(:,1),'r-');hold on;
        P10 = plot(P0{i,fNum}(:,2),P0{i,fNum}(:,1),'w--');hold on;
        title(['Frame: ',int2str(fNum),', Slice: ',int2str(i)]);

        
        subplot(1,3,2+offset),
%         imshow(borderClass(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1 x2],'Xdata',[y1 y2]);hold on;
        imshow(borderClass2(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1 x2],'Xdata',[y1 y2]);hold on;
        plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');
        subplot(1,3,3+offset),
        imshow(Eext2(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1 x2],'Xdata',[y1 y2]);hold on;
        plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');
%             ,['P2 S',int2str(i),', F',int2str(fNum+1)]...
% hL = legend([P21, P22, P11, P12]...        
        hL = legend([P21, P22, P23, P11, P12, P10]...
            ,['P2 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P2 S',sprintf('%02d',i+1),', F',sprintf('%02d',fNum)]...
            ,['P2 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]...
            ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]...
            ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P0 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]);
            
        % Construct a Legend with the data from the sub-plots
        newPosition = [0.5 0.2 0.1 0.1];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits,'Orientation','horizontal', 'fontsize', 6);
        
        folder = 'C:/DocsMaracuya/Datos/Snakes/Snake2DPropagateFNUM';

        
        export_fig(h, strcat(folder,'/Plot_Frame',framenum,'_Slice',slicenum,'.png'));
        export_fig(h, strcat(folder,'/Plot_Slice',slicenum,'_Frame',framenum,'.png'));
        close all;
        
    end
    %     drawnow;
end

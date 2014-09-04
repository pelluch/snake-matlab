% Cluster border on time
% Classify pixels in the border of first snake using k-means cluster
% Then segment from diastole to systole
% jiprieto 1/7/14


clc;
close all;


% clear all;
% % Read VTK 3D
% for ind = 1:12
% % for ind = 1:24
% [V(:,:,:,ind), info] = ReadData3D(strcat('C:/DocsMaracuya/Datos/Snakes/SAguirre/_OriginalReductor_',int2str(ind-1),'.vtk'));
% % [V(:,:,:,ind), info] = ReadData3D(strcat('C:/epxImages/SAguirre/KarisEjeCorto/_OriginalReductor_',int2str(ind-1),'.vtk'));
% end




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
        Options.Alpha = .05; % Membrame energy  (first order), default 0.2
        Options.Beta = .5; %  Thin plate energy (second order), default 0.2
        Options.Delta = .35 %  Baloon force, default 0.1
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
% P2(:,fNum+1) = P2(:,fNum);
%
for fNum = endF:-1:startF
    for nFig =4:1:13
        
        framenum = sprintf('%02d',fNum);
        slicenum = sprintf('%02d',nFig);
        
        
        h=figure('name',strcat('Frame ',int2str(fNum),', Fig',int2str(nFig)),'units','normalized','position',[.1 .1 .5 .5]);
%         set(h,'render','opengl');
%         set(h,'units','normalized','outerposition',[(floor(nFig/10)) 0 1 1]);
        offset = 0;
%         subplot(1,3,1+offset),
        i=nFig; imshow(V(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1  x2],'Xdata',[y1  y2],'InitialMagnification',600);hold on;
        P11 = plot(P1{i,fNum+1}(:,2),P1{i,fNum+1}(:,1),'b--');hold on;
        P12 = plot(P1{i,fNum}(:,2),P1{i,fNum}(:,1),'r-');hold on;
        title(['Frame: ',int2str(fNum),', Slice: ',int2str(i)]);
        hL = legend([ P11, P12]...
            ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]...
            ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]);
            
        % Construct a Legend with the data from the sub-plots
%         newPosition = [0.5 0.2 0.1 0.1];
%         newUnits = 'normalized';
        set(hL,'Location', 'SouthOutside','Orientation','horizontal');%,'Units');%, newUnits,'Orientation','horizontal', 'fontsize', 6);
        
        
        folder = 'C:/DocsMaracuya/Datos/Snakes/Snake2DP1';
        mkdir(folder);
        
        export_fig(h, strcat(folder,'/Plot_Frame',framenum,'_Slice',slicenum,'.bmp'));
%         export_fig(h, strcat(folder,'/Plot_Slice',slicenum,'_Frame',framenum,'.bmp'));
        
%         export_fig(h, strcat(folder,'/Plot_Frame',framenum,'_Slice',slicenum,'.png'));
%         export_fig(h, strcat(folder,'/Plot_Slice',slicenum,'_Frame',framenum,'.png'));
        close all;
        
    end
    %     drawnow;
end

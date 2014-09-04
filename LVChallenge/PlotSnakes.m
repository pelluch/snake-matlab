function PlotSnakes(P0, P1, P2, P1M, P2M, SpringP1, SpringP2, V, border, Eext, OutParam, segLimits, verbose)

% Plot Results
close all;


%Zoom
x1 = OutParam.x1;
x2 = OutParam.x2;
y1 = OutParam.y1;
y2 = OutParam.y2;

if (segLimits.dir == 1)
    index = 1;
else
    index = -1;
end

P2(:,segLimits.endF+1) = P2(:,segLimits.endF);
try
    rmdir(OutParam.writeFolder,'s');
catch
end

mkdir(OutParam.writeFolder);
%
for fNum = segLimits.endF:-1:segLimits.startF
    for nFig =segLimits.startS:1:segLimits.endS
        
        framenum = sprintf('%02d',fNum);
        slicenum = sprintf('%02d',nFig);
        
        
        h=figure('name',strcat('Frame ',int2str(fNum),', Fig',int2str(nFig)));
        %         set(h,'render','opengl');
        
        %         set(h,'units','normalized','outerposition',[(floor(nFig/10)) 0 1 1]);

        subplot(1,3,1),
        
        i=nFig; imshow(V(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1  x2],'Xdata',[y1  y2]);hold on;
        Plot21 = plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');hold on;
%         Plot22 = plot(P2{i+(-1)*index,fNum}(:,2),P2{i+(-1)*index,fNum}(:,1),'y--');hold on;
%         Plot23 = plot(P2{i,fNum+1}(:,2),P2{i,fNum+1}(:,1),'c--');hold on;
        Plot12 = plot(P1{i,fNum}(:,2),P1{i,fNum}(:,1),'r-');hold on;
        Plot10 = plot(P0{i,fNum}(:,2),P0{i,fNum}(:,1),'w--');hold on;
        Plot11 = plot(P1{i,fNum+1}(:,2),P1{i,fNum+1}(:,1),'b--');hold on;
        if(size(P1M{i,fNum},1)>0)
            Plot1M = plot(P1M{i,fNum}(:,1),P1M{i,fNum}(:,2),'c-');hold on;
        end
        if(size(P2M{i,fNum},1)>0)
            Plot2M = plot(P2M{i,fNum}(:,1),P2M{i,fNum}(:,2),'m-');hold on;
        end
        if(size(SpringP1{i,fNum},1)>0)
            PlotS1 = plot(SpringP1{i,fNum}(:,2),SpringP1{i,fNum}(:,1),'yx');hold on;
        end
        if(size(SpringP2{i,fNum},1)>0)
            PlotS2 = plot(SpringP2{i,fNum}(:,2),SpringP2{i,fNum}(:,1),'yx');hold on;
        end
        
        
        title(['Frame: ',int2str(fNum),', Slice: ',int2str(i)]);

        
        subplot(1,3,2),
        imshow(border(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1 x2],'Xdata',[y1 y2]);hold on;
        plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');
        subplot(1,3,3),
        imshow(Eext(x1:x2,y1:y2,i,fNum),[],'Ydata',[x1 x2],'Xdata',[y1 y2]);hold on;
        plot(P2{i,fNum}(:,2),P2{i,fNum}(:,1),'g-');
        
        if(size(P2M{i,fNum},1)>0 && size(P2M{i,fNum},1)>0) 
        hL = legend([Plot21, Plot12, Plot10, Plot1M, Plot2M, Plot11]...
            ,['P2 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P0 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P1M S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P2M S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]);
%             ,['P1 S',sprintf('%02d',i+(-1)*index),', F',sprintf('%02d',fNum)]);
%             ,['P2 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]...
            
        else
            
        hL = legend([Plot21, Plot12, Plot10, Plot11]...
            ,['P2 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P0 S',sprintf('%02d',i),', F',sprintf('%02d',fNum)]...
            ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]);
%             ,['P2 S',sprintf('%02d',i+(-1)*index),', F',sprintf('%02d',fNum)]...
%             ,['P2 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]...
%             ,['P1 S',sprintf('%02d',i),', F',sprintf('%02d',fNum+1)]...

        end
            
        % Construct a Legend with the data from the sub-plots
        newPosition = [0.5 0.2 0.1 0.1];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits,'Orientation','horizontal', 'fontsize', 6);
       
        export_fig(h, strcat(OutParam.writeFolder,'/Plot_Frame',framenum,'_Slice',slicenum,'.bmp'));
        export_fig(h, strcat(OutParam.writeFolder,'/Plot_Slice',slicenum,'_Frame',framenum,'.bmp'));
        if ~verbose
            close all;
        end
        
    end
end
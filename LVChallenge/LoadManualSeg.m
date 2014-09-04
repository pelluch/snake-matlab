function [P1Manual, P2Manual] = LoadManualSeg(folder, subfolder, fileName, nZ, nT)
% Load manual segmentation for LV challenge. The code represents the dicom
% slice where the polygon is situated

fileName = strcat(folder, fileName);
fileList = importdata(fileName);
numContours = size(fileList,1);
arrayPosition = cell(numContours);


for i=1:size(fileList,1)
    index = strfind(fileList{i,1}, '/');
    fname = fileList{i,1}(index(4):end);
    index2 = strfind(fname,'-');
    sliceNum(i) = str2num(fname(index2(2)+1:index2(3)-1));
    arrayPosition{numContours} = [floor(mod(sliceNum(i),nT)), floor(sliceNum(i)/nT)];
    snakeStr = fname(index2(3)+1:index2(4)-1);
    if (strcmp(snakeStr, 'icontour') || strcmp(snakeStr, 'p1contour'))
        snakeType(i) = 1;
    else
        snakeType(i) = 2;
    end
    zindex(i) = floor( (sliceNum(i)-1)/nT )+1;
%     tindex(i) = sliceNum(i) - (zindex(i)-1)*nT;
%     mod( sliceNum(i), nT ) ;
    tindex(i) = mod( sliceNum(i), nT );
    if(tindex (i) ==0)
        tindex (i) = nT;
    end
    
end

P1Manual = cell(nZ, nT);
P2Manual = cell(nZ, nT);

for i=1:numContours
    pfilename = strcat(folder, fileList{i,1}(3:end));
    pfilename(pfilename == '/') = '\';
    if( snakeType(i) == 1)
        P1Manual{zindex(i), tindex(i)} = importdata(pfilename);
    else
        P2Manual{zindex(i), tindex(i)} = importdata(pfilename);
    end
end




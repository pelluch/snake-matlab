function out = resampleImage(im,ref_,varargin)
% out = resampleImage(im,ref)
% out = resampleImage(im,ref,'interpolation',interp)
% out = resampleImage(im,0,'spacing',spacing)
%
% Fast resample an image using a reference image
% Works for ImageType and VectorImageType
%
%   interp is 'NN' (default), 'linear'

interp='NN';
userSpecifiedSpacing=false;
spacing = [-1 -1 -1];
MAX_CHUNK_SIZE = 50;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'interpolation'))
        interp=varargin{i+1};
    elseif (strcmp(varargin{i},'spacing'))
        spacing=varargin{i+1};
        userSpecifiedSpacing=true;
    elseif (strcmp(varargin{i},'maxChunkSize'))
        MAX_CHUNK_SIZE = varargin{i+1};
        i=i+1;
    end
    
end
%----------------------------

if (userSpecifiedSpacing)
    newSize = floor(im.size(:).*im.spacing(:)./spacing(:));
    % positions matrix
    ref = ImageType(newSize,im.origin,spacing,im.orientation);
else
    ref = ref_;
end


% -----------------------------Interpolate ----------------
% create the grid of the reference image
ndims = numel(ref.size);

NCHUNKS = ceil(ref.size/MAX_CHUNK_SIZE);

if (ndims==4)
    
    if (isa(im,'VectorImageType') || isfield(im,'datax'))
        out = VectorImageType(ref);
    else
        out = ImageType(ref);
    end
    
    chunked_size = ceil(out.size./NCHUNKS)';
    
    [ix, iy, iz, it]= ndgrid(0:NCHUNKS(1)-1,0:NCHUNKS(2)-1,0:NCHUNKS(3)-1,0:NCHUNKS(4)-1);
    intervals = [ix(:) iy(:) iz(:) it(:)];
    clear ix iy iz it;
    for i=1:size(intervals,1)
        ranges([1 3 5 7]) = intervals(i,:).*chunked_size+1;
        ranges([2 4 6 8]) = min([(intervals(i,:)+[1 1 1 1]).*chunked_size ; out.size']);
        
        ranges_size = ranges([2 4 6 8])-ranges([1 3 5 7])+[1 1 1 1];
        % generate all the indexes of the target image
        [x, y, z, t] = ndgrid( ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8));
        
        positions = ref.GetPosition([x(:) y(:) z(:) t(:)]');
        clear x y z t;
        datas = im.GetValue(positions,interp);
        
        if isa(im,'VectorImageType') || any(strcmp(properties(im), 'datax'))
            out.datax(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = reshape(datas(1,:),ranges_size);
            out.datay(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = reshape(datas(2,:),ranges_size);
            out.dataz(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = reshape(datas(3,:),ranges_size);
            
            data_ = im.GetValue(positions,interp,'data');
            out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = ...
                reshape(data_,ranges_size);
            clear data positions;
        else
            out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6),ranges(7):ranges(8)) = ...
                reshape(datas,ranges_size);
        end
        clear datas;
        
        
    end
    
elseif (ndims==3)
    %[X Y Z]=ndgrid(1:ref.size(1),1:ref.size(2),1:ref.size(3));
    % Retrieve the world coordinates on the reference image
    
    if (isa(im,'VectorImageType') || isfield(im,'datax'))
        out = VectorImageType(ref);
    else
        out = ImageType(ref);
    end
    
    chunked_size = ceil(out.size./NCHUNKS)';
    
    [ix iy iz]= ndgrid(0:NCHUNKS(1)-1,0:NCHUNKS(2)-1,0:NCHUNKS(3)-1);
    intervals = [ix(:) iy(:) iz(:)];
    clear ix iy iz;
    for i=1:size(intervals,1)
        ranges([1 3 5]) = intervals(i,:).*chunked_size+1;
        ranges([2 4 6]) = min([(intervals(i,:)+[1 1 1]).*chunked_size ; out.size']);
        
        ranges_size = ranges([2 4 6])-ranges([1 3 5])+[1 1 1];
        % generate all the indexes of the target image
        [x y z] = ndgrid( ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6));
        
        positions = ref.GetPosition([x(:) y(:) z(:)]');
        clear x y z;
        datas = im.GetValue(positions,interp);
        
        if isa(im,'VectorImageType') || any(strcmp(properties(im), 'datax'))
            out.datax(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = reshape(datas(1,:),ranges_size);
            out.datay(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = reshape(datas(2,:),ranges_size);
            out.dataz(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = reshape(datas(3,:),ranges_size);
            
            data_ = im.GetValue(positions,interp,'data');
            out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = ...
                reshape(data_,ranges_size);
            clear data positions;
        else
            out.data(ranges(1):ranges(2),ranges(3):ranges(4),ranges(5):ranges(6)) = ...
                reshape(datas,ranges_size);
        end
        clear datas;
        
        
    end
elseif(ndims==2)
    [X, Y]=ndgrid(1:ref.size(1),1:ref.size(2));
    % Retrieve the world coordinates on the reference image
    positions = ref.GetPosition([X(:) Y(:)]');
    
    if isa(im,'VectorImageType') ||  any(strcmp(properties(obj), 'datax'))
        out = VectorImageType(ref);
        datas = im.GetValue(positions,interp);
        out.datax = reshape(datas(1,:),ref.size(1),ref.size(2));
        out.datay = reshape(datas(2,:),ref.size(1),ref.size(2));
        % out.dataz = reshape(datas(3,:),ref.size(1),ref.size(2));
    else
        out = ImageType(ref);
        out.data = reshape(im.GetValue(positions,interp),ref.size(1),ref.size(2));
    end
    
end
out.data(out.data~=out.data)=0;

end
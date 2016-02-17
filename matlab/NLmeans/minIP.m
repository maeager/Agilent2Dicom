function [mIP] = minIP(mag, slab, dim)
%%minimum intensity projection
% minimum projection of 3D sections slab  
% 
% -(C) 2015  Michael Eager (michael.eager@monash.edu)
mIP=zeros(size(mag));
if nargin < 3 || dim == 3
    for i =1:size(mag,3)
        start = i;endp = min(size(mag,3), i+slab);
        mIP(:,:,i) = min(mag(:,:,start:endp),[],3);
    end
else
    if dim == 1
        
        for i =1:size(mag,1)
            start = i;endp = min(size(mag,1), i+slab);
            mIP(:,:,i) = squeeze(min(mag(start:endp,:,:),[],1));
        end
    elseif dim == 2
        for i =1:size(mag,2)
            start = i;endp = min(size(mag,2), i+slab);
            mIP(:,:,i) = squeeze(min(mag(:,start:endp,:),[],2));
        end
    end
end

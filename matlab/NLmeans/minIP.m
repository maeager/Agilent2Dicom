function [mIP]=minIP(mag,slab)
%%minimum intensity projection
% 2D minimum of 3D slab  
% 
% -(C) 2015  Michael Eager (michael.eager@monash.edu)


for i =1:size(mag,3)
    start = i;endp = min(size(mag,3), i+slab);
    mIP(:,:,i) = min(mag(:,:,start:endp),[],3);
end

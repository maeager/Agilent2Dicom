function [img, Range] = NormaliseImage2(img,Range)
%% Simple normalisation of N-dimensional image
%  
%
%   - Michael Eager,   (michael.eager@monash.edu) 
%   - (c) 2012, Monash Biomedical Imaging, Monash University, Australia

if nargin == 1
  Range(1) = min(img(:));
  Range(2) = max(img(:));
else
if length(Range) ~= 2 || Range(2) <= Range(1)
    disp 'NormaliseImage2 error - range incorrect'
    return
end
  img(img<Range(1))=Range(1);
  img(img>Range(2))=Range(2);
end
% TODO :  make the process more effieicent - do we need to create another img

%if min_image ~= 0 || max_image ~= 1
  img = (img - Range(1))./(Range(2) - Range(1));
%end

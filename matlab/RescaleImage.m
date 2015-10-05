function img = RescaleImage(img,Range)
%% Simple reverse normalisation of N-dimensional image
%  
%
%   - Michael Eager,   (michael.eager@monash.edu) 
%   - (c) 2012, Monash Biomedical Imaging, Monash University, Australia


  img = img .*(Range(2) - Range(1)) + Range(1);

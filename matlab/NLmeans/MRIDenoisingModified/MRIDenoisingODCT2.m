function [ imaODCT] = MRIDenoisingODCT2(ima, sigma, rician, verbose)
%
%   Description: Denoising of a 3D MRI image using the ODCT filter
%               
%   Usage:      MRIDenoisingODCT(ima, sigma, rician, verbose)
%
%   ima:        3D MR image
%   sigma:          std of noise
%   rician:     0: Gaussian noise 1: Rician noise 
%   verbose:    0 no display of graph
%               1 display of graph
%
%                          Details on ODCT filter
%**************************************************************************
%* The ODCT filter is described in:                                       *
%*                                                                        *
%* JV Manjon, P. Coupe, A Buades, DL Collins, M. Robles                   *
%* New methods for MRI denoising based on sparseness and self-similarity. * 
%* Medical Image Analysis, 16(1):18-27                                    *
%* 2012                                                                   *
%**************************************************************************
%
% Pierrick Coupe - pierrick.coupe@labri.fr
% Jose V. Manjon - jmanjon@fis.upv.es
%
% Universite de Bordeaux
% LaBRI, UMR 5800, 33400 Talence
%
% Copyright (C) 2011 Pierrick Coupe and Jose V. Manjon

s=size(ima);

if (size(s,2)~=3)
    disp('your image is not a 3D image')
    imaODCT = 0;
    return
end
    if nargin < 2 || isempty(sigma) || (sigma==0)
        sigma=1;
    end
    
    if  nargin < 3 ||  isempty(rician) || (rician ~= 0)
        rician=1;
    end
    if  nargin < 4
        verbose=0;
    end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Denoising using ODCT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
disp('.')
disp('Denoising using ODCT')
disp('.')


disp('           Details on ODCT filter')
disp('**************************************************************************')
disp('* The ODCT filter is described in:                                       *')
disp('*                                                                        *')
disp('* JV Manjon, P. Coupe, A Buades, DL Collins, M. Robles                   *')
disp('* New methods for MRI denoising based on sparseness and self-similarity. *')
disp('* Medical Image Analysis, 16(1):18-27                                    *')
disp('* Avril 2012                                                             *')
disp('**************************************************************************')
 
  %}  
    
    imaODCT=myODCT3d(single(ima),single(sigma), rician);
    map = find(imaODCT<0);
    imaODCT(map) =0;

   
    
    if(verbose==1)
       figure;
        mini = min(ima(:));
        maxi = max(imaODCT(:));
        
       
            subplot(1,3,1)
            imagesc(ima(:,:,floor(s(3)/2)))%,[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Noisy Image');
            title(tit)
            subplot(1,3,2)
            imagesc(imaODCT(:,:,floor(s(3)/2)))%,[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Denoised Image');
            title(tit)
            subplot(1,3,3)
            imagesc(abs(imaODCT(:,:,floor(s(3)/2)) - ima(:,:,floor(s(3)/2))))
            axis image;
            axis off;
            tit = sprintf('Removed Noise');
            title(tit)
            colormap(gray)
            drawnow
        
    end

end


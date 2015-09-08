function [ MRORNLM] = MRIDenoisingAONLM2(ima, patchsize, ...
                                         searcharea, rician, coil, verbose)
%
%   Description: Denoising of a 3D MRI image using the multiresolution AONLM
%   filter.. Using HSM wavelet mixing
%         
%   Usage:      MRIDenoisingAONLM(ima, sigma, beta, patchsize, searcharea, rician, verbose)
%
%   ima:        3D MR image
%   sigma:          std of noise
%   beta:       smooting parameter (default 1)
%   patchsize:  radius of patch in voxel
%   searcharea: radius of search area in voxel
%   rician:     0: Gaussian noise 1: Rician noise 
%   verbose:    0 no display of graph
%               1 display of graph
%
%                          Details on MRONLM filter
%**************************************************************************
%* The AONLM filter is described in:                                     *
%*                                                                        *
%* J. V. Manjon, P. Coupe, L. Marti-Bonmati, D. L. Collins, M. Robles.    *
%* Adaptive Non-Local Means Denoising of MR Images with Spatially         *
%* varying Noise Levels.                                                  *
%* Journal of Magnetic Resonance Imaging, 31(1):192-203, 2010             *
% 
%* P. Coupe, J.V. Manjon, R. Robles, D.L. Collins                         *
%* Adaptive Multiresolution Non-Local Means Filter for 3D MR Image        *
%* Denoising. IET Image Processing.                                       *
%* 2012                                                                   *
%**************************************************************************

% Pierrick Coupe - pierrick.coupe@labri.fr
% Jose V. Manjon - jmanjon@fis.upv.es
%
% Universit� de Bordeaux
% LaBRI, UMR 5800, 33400 Talence
%
% Copyright (C) 2011 Pierrick Coupe and Jose V. Manjon
%
% Modifications: Michael Eager, Monash Biomedical Imaging
% 

s=size(ima);

if (size(s,2)~=3)
    disp('your image is not a 3D image')
    AORNLM = 0;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Denoising each image using MRORNLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
disp('.')
disp('Denoising using Multi-resolution AONLM')
disp('.')

disp('                          Please cite')
disp('**************************************************************************')
disp('* The AONLM filter is described in:                                      *')
disp('*                                                                        *')
disp('* J. V. Manjon, P. Coupe, L. Marti-Bonmati, D. L. Collins, M. Robles.    *')
disp('* Adaptive Non-Local Means Denoising of MR Images with Spatially         *')
disp('* varying Noise Levels.                                                  *')
disp('* Journal of Magnetic Resonance Imaging, 31(1):192-203, 2010             *')
disp('**************************************************************************')
disp('.')
  %}  
    if nargin < 2 || isempty(patchsize) || (patchsize==0)
        patchsize=1;
    end
    
    if  nargin < 3 ||  isempty(searchsize) || (searcharea==0)
        searcharea=3;
    end
 
    if  nargin < 4 ||  isempty(rician) || (rician ~= 0)
        rician=1;
    end
    if  nargin < 5 || isempty(coil)
        %% create coil-sensitivity matrix here
        coil = ones(size(ima));
    end
    if  nargin < 6 
        verbose=0;
    end
    
    
    ORNLMu=myMBONLM3d(single(ima),searcharea,patchsize,sigma, ...
                    rician,coil);
    
    ORNLMo=myMBONLM3d(single(ima),searcharea,patchsize+1,sigma, ...
                    rician,coil);
    
    AORNLM = hsm(ORNLMu, ORNLMo); 
    map = find(AORNLM<0);
    AORNLM(map) =0;
    
    
    if(verbose==1)
       figure;
        mini = min(ima(:));
        maxi = max(AORNLM(:));
        
       
            subplot(1,3,1)
            imagesc(ima(:,:,floor(s(3)/2)),[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Noisy Image');
            title(tit)
            subplot(1,3,2)
            imagesc(AORNLM(:,:,floor(s(3)/2)),[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Denoised Image');
            title(tit)
            subplot(1,3,3)
            imagesc(abs(AORNLM(:,:,floor(s(3)/2))-ima(:,:,floor(s(3)/2))))
            axis image;
            axis off;
            tit = sprintf('Denoised Image');
            title(tit)
            colormap(gray)
            drawnow
        
    end

end


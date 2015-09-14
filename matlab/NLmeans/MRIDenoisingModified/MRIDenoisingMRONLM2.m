function [ MRORNLM] = MRIDenoisingMRONLM2(ima, sigma, patchsize, ...
                                         searcharea, rician, coil, verbose)
%
%   Description: Denoising of a 3D MRI image using the multiresolution ONLM
%   filter
%         
%   Usage:      MRIDenoisingMRONLM(ima, sigma, beta, patchsize, searcharea, rician, verbose)
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
%* The MRONLM filter is described in:                                     *
%*                                                                        *
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

s=size(ima);

if (size(s,2)~=3)
    disp('your image is not a 3D image')
    MRORNLM = 0;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Denoising each image using MRORNLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
disp('.')
disp('Denoising using Multi-resolution ONLM')
disp('.')

disp('                          Please cite')
disp('**************************************************************************')
disp('* The multi-resolution ONLM filter is described in:                      *')
disp('*                                                                        *')
disp('* P. Coupe, J. V Manjon, M. Robles, D. L. Collins                        *')
disp('* Adaptive Multiresolution Non-Local Means Filter for three-dimensional  *')
disp('* magnetic resonance image denoising.                                    *')
disp('* IET Image Processing. 6(5): 558-568, 2012                              *')
disp('**************************************************************************')
disp('.')
  %}  
    if (patchsize==0)
        patchsize=1;
    end
    
    if (searcharea==0)
        searcharea=3;
    end
    
    if isempty(coil)
        %% create coil-sensitivity matrix here
        coil = ones(size(ima));
    end
    
    ORNLMu=myMBONLM3D(single(ima),searcharea,patchsize,sigma, ...
                    rician,coil);
    
    ORNLMo=myMBONLM3D(single(ima),searcharea,patchsize+1,sigma, ...
                    rician,coil);
    
    MRORNLM = ascm(ima, ORNLMu, ORNLMo, sigma); 
    map = find(MRORNLM<0);
    
    MRORNLM(map) =0;
    
    
    if(verbose==1)
       figure;
        mini = min(ima(:));
        maxi = max(MRORNLM(:));
        
       
            subplot(1,3,1)
            imagesc(ima(:,:,floor(s(3)/2)))%,[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Noisy Image');
            title(tit)
            subplot(1,3,2)
            imagesc(MRORNLM(:,:,floor(s(3)/2)))%,[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Denoised Image');
            title(tit)
            subplot(1,3,3)
            imagesc(abs(MRORNLM(:,:,floor(s(3)/2))-ima(:,:,floor(s(3)/2))))
            axis image;
            axis off;
            tit = sprintf('Denoised Image');
            title(tit)
            colormap(gray)
            drawnow
        
    end

    




end


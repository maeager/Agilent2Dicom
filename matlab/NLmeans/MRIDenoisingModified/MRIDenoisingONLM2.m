function [ ORNLM] = MRIDenoisingONLM2(ima, sigma, patchsize, searcharea, rician, verbose)
%
%   Description:  Denoising of a 3D MRI image using the ONLM filter
%               
%
%   Usage:      MRIDenoisingONLM(ima, h, patchsize, searcharea, rician, verbose)
%
%   ima:        3D MR image
%   sigma:          std of rician noise
%   patchsize:  radius of patch in voxel
%   searcharea: radius of search area in voxel
%   rician:     0: Gaussian noise 1: Rician noise 
%   verbose:    0 no display of graph
%               1 display of graph
%
%                          Details on ONLM filter
%**************************************************************************
%* The ONLM filter is described in:                                       *
%*                                                                        *
%* P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
%* An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
%* Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, *
%* Avril 2008                                                             *
%**************************************************************************
%                      Details on Rician adaptation
%**************************************************************************
%* The adaptation to Rician noise is described in:                        *
%*                                                                        *
%* N. Wiest-Daessle, S. Prima, P. Coupe, S.P. Morrissey, C. Barillot.     *
%* Rician noise removal by non-local means filtering for low              *
%* signal-to-noise ratio MRI: Applications to DT-MRI. In 11th             *
%* International Conference on Medical Image Computing and                *
%* Computer-Assisted Intervention, MICCAI'2008,                           *
%* Pages 171-179, New York, Etats-Unis, Septembre 2008                    *
%**************************************************************************
%
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
    ORNLM = 0;
    return
end

    if nargin < 2 || isempty(sigma) || (sigma==0)
        sigma=1;
    end

    if nargin < 3 || isempty(patchsize) || (patchsize==0)
        patchsize=1;
    end
    
    if  nargin < 4 ||  isempty(searcharea) || (searcharea==0)
        searcharea=3;
    end
    

    if  nargin < 5 ||  isempty(rician) || (rician ~= 0)
        rician=1;
    end
    if  nargin < 6
        verbose=0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Denoising each image using ORNLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
disp('.')
disp('Denoising using ORNLM')
disp('.')

disp('                          Please cite')
disp('**************************************************************************')
disp('* The ONLM filter is described in:                                       *')
disp('*                                                                        *')
disp('* P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *')
disp('* An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*')
disp('* Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, *')
disp('* Avril 2008                                                             *')
disp('**************************************************************************')
disp('                              and')
disp('**************************************************************************')
disp('* The adaptation to Rician noise is described in:                        *')
disp('*                                                                        *')
disp('* N. Wiest-Daessle, S. Prima, P. Coupe, S.P. Morrissey, C. Barillot.     *')
disp('* Rician noise removal by non-local means filtering for low              *')
disp('* signal-to-noise ratio MRI: Applications to DT-MRI. In 11th             *')
disp('* International Conference on Medical Image Computing and                *')
disp('* Computer-Assisted Intervention, MICCAI''2008,                           *')
disp('* Pages 171-179, New York, Etats-Unis, Septembre 2008                    *')
disp('**************************************************************************')
disp('.')
%}  

%        if isempty(coil)
        %% create coil-sensitivity matrix here
        coil = single(ones(size(ima)));
%    end
    
    ORNLM=myMBONLM3D(single(ima),searcharea,patchsize,single(sigma),rician,coil);
    map = find(ORNLM<0);
    ORNLM(map) =0;

   
    
    if(verbose==1)
       figure;
        mini = min(ima(:));
        maxi = max(ORNLM(:));
        
       
            subplot(1,3,1)
            imagesc(ima(:,:,floor(s(3)/2)))%,[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Noisy Image');
            title(tit)
            subplot(1,3,2)
            imagesc(ORNLM(:,:,floor(s(3)/2)))%,[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Denoised Image');
            title(tit)
             subplot(1,3,3)
            imagesc(abs(ORNLM(:,:,floor(s(3)/2))-ima(:,:,floor(s(3)/2))))
            axis image;
            axis off;
            tit = sprintf('Removed noise');
            title(tit)
            colormap(gray)
            drawnow
        
    end

    




end


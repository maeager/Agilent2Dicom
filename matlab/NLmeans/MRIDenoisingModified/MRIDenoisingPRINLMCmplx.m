function [ imaPRINLM] = MRIDenoisingPRINLMCmplx(cmplxima, sigma, patchsize, ...
                                            searcharea, beta, rician, verbose)
%
%   Description:  Denoising of a 3D MRI image using the PRINLM filter
%               
%
%   Usage:      MRIDenoisingPRINLM(ima, h, rician, verbose)
%
%   ima:        3D MR image
%   h:          std of noise
%   rician:     0: Gaussian noise 1: Rician noise
%   verbose:    0 no display of graph
%               1 display of graph
%
%                          Details on PRINLM filter
%**************************************************************************
%* The PRINLM filter is described in:                                     *
%*                                                                        *
%*JV Manjon, P. Coupe, A Buades, DL Collins, M. Robles                    *
% New methods for MRI denoising based on sparseness and self-similarity.  * 
% Medical Image Analysis, 16(1):18-27                                     *
%* Avril 2012                                                             *
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
    imaPRINLM = 0;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Denoising using ODCT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (patchsize==0)
        patchsize=1;
    end
    
    if (searcharea==0)
        searcharea=3;
    end
    
    if (beta==0)
        beta=1;
    end
    

 
%{    
    disp('.')
disp('Denoising using ODCT')
disp('.')

%}
    imaODCTr=myODCT3d(single(real(cmplxima)),sigma,rician);
    map = find(imaODCTr<0);
    imaODCTr(map) =0;  
    if ~isreal(cmplxima))
    imaODCTi=myODCT3d(single(imag(cmplxima)),sigma,rician);
    map = find(imaODCTi<0);
    imaODCTi(map) =0;
    else
    imaODCTi=imaODCTr*0;
    end
    
%{
    disp('.')
disp('Denoising using PRINLM')
disp('.')

%}
    imaPRINLM=myRINLM3d(single(abs(cmplxima)),searcharea,patchsize,sigma,abs(complex(imaODCTr,imaODCTi)),rician); 
    map = find(imaPRINLM<0);
    imaPRINLM(map) = 0;

%{
disp('                     Please cite')
disp('**************************************************************************')
disp('* The PRINLM filter is described in:                                     *')
disp('*                                                                        *')
disp('* JV Manjon, P. Coupe, A Buades, DL Collins, M. Robles                   *')
disp('* New methods for MRI denoising based on sparseness and self-similarity. *')
disp('* Medical Image Analysis, 16(1):18-27                                    *')
disp('* Avril 2012                                                             *')
disp('**************************************************************************')
 
  %}  
   
    
    if(verbose==1)
       figure;
        mini = min(ima(:));
        maxi = max(imaPRINLM(:));
        
       
            subplot(1,3,1)
            imagesc(ima(:,:,floor(s(3)/2))) %,[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Noisy Image');
            title(tit)
            subplot(1,3,2)
            imagesc(imaPRINLM(:,:,floor(s(3)/2))) %,[mini maxi-0.25*maxi])
            axis image;
            axis off;
            tit = sprintf('Denoised Image');
            title(tit)
             subplot(1,3,3)
            imagesc(abs(imaPRINLM(:,:,floor(s(3)/2))-ima(:,:,floor(s(3)/2))))
            axis image;
            axis off;
            tit = sprintf('Removed noise');
            title(tit)
            colormap(gray)
            drawnow
        
    end

    




end


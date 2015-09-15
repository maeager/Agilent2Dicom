function [ imaPRINLM, imaODCT, coil] = MRIDenoisingPRINLMCmplx(cmplxima, sigma, patchsize, searcharea, rician,coil, verbose)
%
%   Description:  Denoising of a 3D MRI image using the PRINLM filter
%               
%
%   Usage:      MRIDenoisingPRINLM(ima, sigma, rician, verbose)
%
%   ima:        3D MR image
%   sigma:          std of noise
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

if nargin < 2 || isempty(sigma) || (sigma==0)
    sigma=1;
end

if nargin < 3 || isempty(patchsize) || (patchsize==0)
    patchsize=1;
end

if  nargin < 4 ||  isempty(searchsize) || (searcharea==0)
    searcharea=3;
end

if  nargin < 6 ||  isempty(rician) || (rician ~= 0)
    rician=1;
end
if  nargin < 7 || isempty(coil)
coil = single(ones(size(ima)));
end
if  nargin < 8 
    verbose=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Denoising using ODCT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{    
    disp('.')
    disp('Denoising using ODCT')
    disp('.')

%}

    imaODCTr=myODCT3d(single(real(cmplxima)),sigma,rician); 
    if ~isreal(cmplxima))
    imaODCTi=myODCT3d(single(imag(cmplxima)),single(sigma),rician);

    else
        imaODCTi=imaODCTr*0;
    end
    imaODCT = abs(complex(imaODCTr,imaODCTi))
    map = find(imaODCT<0);
    imaODCT(map) =0; 
%{
    disp('.')
    disp('Denoising using PRINLM')
    disp('.')

%}

imaPRINLM=myRINLM3d(single(abs(cmplxima)),searcharea,patchsize,single(sigma),single(imaODCT),rician, single(coil)); 
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



function [ imaPRINLM, imaODCT, B1coil,imaPRINLMc] = MRIDenoisingB1PRINLM2(real_ima,imag_ima, real_sigma,imag_sigma, patchsize, searcharea, rician,coil, verbose)
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

s=size(real_ima);

if (size(s,2)~=3)
    disp('Real image is not a 3D image')
    imaPRINLM = 0;
    return
end

s=size(imag_ima);

if (size(s,2)~=3)
    disp('Imag image is not a 3D image')
    imaPRINLM = 0;
    return
end


if nargin < 3 || isempty(real_sigma) || (real_sigma==0)
    real_sigma=1;
end

if nargin < 4 || isempty(imag_sigma) || (imag_sigma==0)
    imag_sigma=1;
end

if nargin < 5 || isempty(patchsize) || (patchsize==0)
    patchsize=1;
end

if  nargin < 6 ||  isempty(searchsize) || (searcharea==0)
    searcharea=3;
end

if  nargin < 7 ||  isempty(rician) || (rician ~= 0)
    rician = 0;
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
 
 
 % Step (1): Scale real and imag
  [norm_real,realRange] = NormaliseImage2(real(cmplxima));
  [norm_imag,imagRange] = NormaliseImage2(imag(cmplxima));


 
imaODCTr=myODCT3d(single(real_ima),single(real_sigma),rician); 
if ~isempty(imag_ima)
    min_imag = min(imag_ima(:));
    imag_ima = imag_ima - min_imag;
    imaODCTi=myODCT3d(single(imag_ima),single(imag_sigma),rician);
else
    imaODCTi=imaODCTr*0;
end
imaODCT = abs(complex(imaODCTr,imaODCTi+min_imag));
map = find(imaODCT<0);
imaODCT(map) =0; 



%{
disp('.')
disp('Denoising using PRINLM')
disp('.')
%}

% imaPRINLM=myRINLM3d(single(abs(complex(real_ima,imag_ima))),searcharea,patchsize,single(sqrt(real_sigma^2+imag_sigma^2)),single(imaODCT),rician, single(B1)); 
% map = find(imaPRINLM<0);
% imaPRINLM(map) = 0;

imaPRINLMr=myRINLM3d(single(real_ima),searcharea,patchsize,single(real_sigma),single(imaODCT), rician, single(coil)); 
display 'Starting imaginary RINLM3D'
imaPRINLMi=myRINLM3d(single(imag_ima),searcharea,patchsize,single(imag_sigma),single(imaODCT), rician, single(coil)); 

imaPRINLMc=complex(imaPRINLMr,imaPRINLMi+min_imag);


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



function [denoised_img,swi_n2,swi_n1,swi_nOrig,swi_n0,sigma_real,sigma_imag] = cmplx_pipeline1(mag,phase,NLfilter,searcharea,patchsize,rician,hfactor,sigma,noisepc,B1bias,fgamma)
%% cmplx_pipeline1: Complex non-local means denoising Option 1
%  Option 1 calls the automatic noise estimate before running the
%  NLmeans filter.  This function normalises real and imag components of
%  complex images before processing
%
%% Input Args:
%      mag           Magnitude of input image
%      phase         Phase of input image
%      NLfilter      Set NL denoiseing method
%      searcharea    NL search area
%      patchsize     Patch size
%      rician        Bool flag to enable Rician distribution in estimates
%      hfactor       Factor to scale sigma
%      noisepc       Add noise to input image as percentage (0.05 to 0.5
%                    appropriate)
%      B1bias        Use B1 bias correction (3D image same size as mag)
%      fgamma        kernel gamma function (0=MINMAX,1=MULT,2=EXPDR,3=DR,4=DIFF)
%
% Output Args:
% denoised_img  Complex denoised image
% swi_n2        SWI of the denoised mag and denoised phase
% swi_n1        SWI of the denoised mag and the Orig phase (if noise added, Orig phase will have changed)
% swi_n0        SWI of the input mag and phase (if noise added, mag and phase will have changed)
% swi_nOrig     SWI of the input mag and phase (no effect of added noise)
% sigma_real    Estimated noise in the real component of the input image
% sigma_imag    Estimated noise in the imag component of the input image
%
% Example usage:
%  mag = single(t2swiMag);
%  pha_max = nextpow2(max(abs(t2swiPha)));
%  phase = single(t2swiPha*pi/pha_max);
%  [denoised_img,swi_denoised,swi_withOrigphase] = cmplx_pipeline1(mag,phase,-1,5,1,1,75,5)
%
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
% -     Monash Biomedical Imaging

narginchk(2,11);

if max(abs(phase(:))) > pi
    display('Correcting phase to nearest 2^power');
    pha_max = nextpow2(max(abs(phase(:))));
    phase = single(phase*pi/(2^pha_max));
end

img = mag.*exp(1i.*phase);
[swi_nOrig, swi_pOrig,pmaskOrig] = swi_nohdyne(mag,-phase,4);
real_img=real(img);imag_img=imag(img);

sigma_real =[];
sigma_imag=[];
if nargin >= 8 && ~isempty(sigma) && numel(sigma) == 2
    sigma_real = sigma(1);
    sigma_imag = sigma(2);
end

if nargin < 7
    hfactor=[];
end

if nargin < 3 || isempty(NLfilter)
    NLfilter =-1;
end

if nargin < 4 || isempty(searcharea)
    if NLfilter == -5 || NLfilter == 0 ||NLfilter == -2 || NLfilter == 2
        searcharea = 3;
    else
        searcharea = 5;
    end
end
if nargin < 5 || isempty(patchsize)
    patchsize = 1;
end
if nargin < 6 || isempty(rician)
    rician = 0;
end



% add some noise
if nargin >= 9 && ~isempty(noisepc) && noisepc > 0
    real_img = real_img + (randn(size(real_img))*noisepc*mean(mag(:)));
    imag_img = imag_img + (randn(size(imag_img))*noisepc*mean(mag(:)));
    %Adjust centre of histogrma in real and imag images
    %real_img=real_img-mode(real_img(:));
    %imag_img=imag_img-mode(imag_img(:));
    %reset img, mag and phase with noise
    
    img = complex(real_img,imag_img);
    mag = abs(img);
    phase = angle(img);
    
end

if nargin < 10
    B1bias=single(ones(size(img))); 
end
if nargin < 11 || isempty(fgamma)
   fgamma = 0; 
end
%ima=t2swi;
% Step (2): Scale real and imag

[norm_real_img,realRange] = NormaliseImage2(real_img);
[norm_imag_img,imagRange] = NormaliseImage2(imag_img);

norm_real_img=norm_real_img*256.0;
norm_imag_img=norm_imag_img*256.0;  

% Denoise using pipeline 
tic()
[MRIdenoisedReal,sigma_real,filtername] = pipeline1(norm_real_img, NLfilter,...
    sigma_real, hfactor, searcharea, patchsize, rician,B1bias,fgamma);toc
[MRIdenoisedImag,sigma_imag,filtername] = pipeline1(norm_imag_img, NLfilter,...
    sigma_imag, hfactor, searcharea, patchsize, rician,B1bias,fgamma);
toc()

% Step (4) rescale
MRIdenoisedRealRescaled = MRIdenoisedReal*(realRange(2)-realRange(1))/256.0  + realRange(1);
display (['Real Rescaled mode ' num2str(mode(MRIdenoisedRealRescaled(:)))])
MRIdenoisedRealRescaled = MRIdenoisedRealRescaled - mode(MRIdenoisedRealRescaled(:));
MRIdenoisedImagRescaled = MRIdenoisedImag*(imagRange(2) - imagRange(1))/256.0 + imagRange(1);
display (['Imag Rescaled mode ' num2str(mode(MRIdenoisedImagRescaled(:)))])
MRIdenoisedImagRescaled = MRIdenoisedImagRescaled - mode(MRIdenoisedImagRescaled(:));

denoised_img = complex(MRIdenoisedRealRescaled, MRIdenoisedImagRescaled);
% Step (5) Save denoised images 

% Save denoised SWI
swi_n0= swi_nohdyne(mag, -phase, 4);
swi_n1 = swi_nohdyne(abs(denoised_img), -phase, 4);
swi_n2 = swi_nohdyne(abs(denoised_img), -angle(denoised_img), 4);

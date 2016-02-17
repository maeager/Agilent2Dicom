function [MRIdenoised,sigma, filtername] = pipeline1(img,NLfilter,hfinal,hfactor,searcharea,patcharea,rician,B1bias,fgamma)
%% Non-local means denoising Option 1
%  Option 1 calls the automatic noise estimate before running the
%  NLmeans filter
%
% Input Args:
%      img           Input image (magn or complex, 3-5 dimensions)
%      NLfilter      Set NL denoiseing method
%      hfinal        Sigma, estimate of noise
%      hfactor       Factor to scale sigma
%      searcharea    NL search area
%      patchsize     Patch size
%      rician        Bool flag to enable Rician distribution in estimates
%      B1bias        Use B1 bias correction (3D image same size as mag)
%      fgamma        kernel gamma function (0=MINMAX,1=MULT,2=EXPDR,3=DR,4=DIFF)
%
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
% -     Monash Biomedical Imaging

%% parse inputs
narginchk(2, 9)
beta_=1;


if nargin < 9 || isempty(fgamma)
   fgamma = 0; 
end
if nargin < 8 || isempty(B1bias)
   B1bias = 0; 
end
if nargin < 7 || isempty(rician)
    rician = single(1);
end
if nargin < 6 || isempty(patcharea)
    patcharea = 1;
end

if nargin < 5 || isempty(searcharea)
  if NLfilter == -5 || NLfilter == 0 
    searcharea = 3;
  else
    searcharea = 5;
  end
end
if nargin < 4 || isempty(hfactor)
    hfactor=100;
end
if nargin < 3 || isempty(hfinal)
    display 'Calculating noise estimate.'
    [hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(img,rician,0);
end

hfinal = hfinal * (hfactor/100);
display(['HFINAL: ' num2str(hfinal)]);

%% run filter
display(['Noise estimate: ' num2str(hfinal)])

sigma=hfinal;
filtername='';

switch NLfilter
 case 0 
 display('Processing denoised image - MRONLM')
 tic(),MRIdenoised = MRIDenoisingMRONLM(img, sigma, beta_,patcharea,...
        searcharea, rician, B1bias); toc()
        				filtername='MRONLM';  
 case 1
 display('Processing denoised image - PRINLM')
 tic(),MRIdenoised = MRIDenoisingPRINLM(img, sigma,patcharea,...
				       searcharea, beta_,rician, B1bias);toc()
					filtername='PRINLM';
 case 2
 display('Processing denoised image - AONLM')
 tic(),MRIdenoised = MRIDenoisingAONLM(img, beta_, patcharea, ...
				       searcharea, rician, B1bias);toc()
				       filtername='AONLM';
 case 3
 display('Processing denoised image - ONLM ')
 tic(),MRIdenoised = MRIDenoisingONLM(img, sigma, ...
				      beta_, patcharea, searcharea, ...
				      rician , B1bias);toc()  
				      filtername='ONLM';
 case 4
 display('Processing denoised image - ODCT ')
 tic(),MRIdenoised = MRIDenoisingODCT(img, ...
				      sigma, ...
				      beta_,rician,B1bias);toc()   
				      filtername='ODCT';  
 case -5 
 display('Processing denoised image - MRONLM2')
 tic(),MRIdenoised = MRIDenoisingMRONLM2(img, sigma, ...
					 patcharea, searcharea, ...
					 rician, B1bias,fgamma);toc() 
					 filtername='MRONLM2';
 case -1
 display('Processing denoised image - PRINLM2')
 tic(),MRIdenoised = MRIDenoisingPRINLM2(img, sigma, ...
                                         patcharea, searcharea, ...
					 rician, B1bias,fgamma);toc()
					 filtername='PRINLM2';
 case -2
 display('Processing denoised image - AONLM2')
 tic(),MRIdenoised = MRIDenoisingAONLM2(img, patcharea, ...
					searcharea, rician, B1bias,fgamma);toc()
					filtername='AONLM2';
 case -3
 display('Processing denoised image - ONLM2 ')
 tic(),MRIdenoised = MRIDenoisingONLM2(img, sigma, ...
				       patcharea, searcharea, ...
				       rician , B1bias,fgamma);toc()  
				       filtername='ONLM2';
 case -4
 display('Processing denoised image - ODCT2 ')
 tic(),MRIdenoised = MRIDenoisingODCT2(img, ...
				       sigma, ...
				       rician,B1bias,fgamma);toc()   
				       filtername='ODCT2';

 otherwise
display('Processing Real denoised image - MRONLM')
tic(),MRIdenoised = MRIDenoisingMRONLM(img,sigma,...
				       beta_, patcharea, searcharea, ...
				       rician,B1bias);toc()
				       filtername='MRONLM';
end

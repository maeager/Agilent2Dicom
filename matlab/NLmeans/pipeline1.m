function [MRIdenoised,sigma, filtername] = pipeline1(img,NLfilter,hfinal,hfactor,searcharea,patcharea,rician,B1bias)
%% Non-local means denoising Option 1
%  Option 1 calls the automatic noise estimate before running the
%  NLmeans filter
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
% -     Monash Biomedical Imaging

%% parse inputs

if nargin < 4 || isempty(hfactor)
    hfactor=100;
end
hfinal = hfinal * (hfactor/100)

if nargin < 5 || isempty(searcharea)
  if NLfilter == -5 || NLfilter == 0 
    searcharea = 3;
  else
    searcharea = 5;
  end
end
if nargin < 6 || isempty(patcharea)
    patcharea = 1;
end
if nargin < 7 || isempty(rician)
    rician = 1;
end
beta_=1;
if nargin < 8 || isempty(B1bias)
   B1bias = 0; 
end
if nargin < 3 || isempty(hfinal)
    display 'Calculating noise estimate.'
    [hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(img,rician,0)
end


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
					 rician, B1bias);toc() 
					 filtername='MRONLM2';
 case -1
 display('Processing denoised image - PRINLM2')
 tic(),MRIdenoised = MRIDenoisingPRINLM2(img, sigma,  ...
					 rician, B1bias);toc()
					 filtername='PRINLM2';
 case -2
 display('Processing denoised image - AONLM2')
 tic(),MRIdenoised = MRIDenoisingAONLM2(img, patcharea, ...
					searcharea, rician, B1bias);toc()
					filtername='AONLM2';
 case -3
 display('Processing denoised image - ONLM2 ')
 tic(),MRIdenoised = MRIDenoisingONLM2(img, sigma, ...
				       patcharea, searcharea, ...
				       rician , B1bias);toc()  
				       filtername='ONLM2';
 case -4
 display('Processing denoised image - ODCT2 ')
 tic(),MRIdenoised = MRIDenoisingODCT2(img, ...
				       sigma, ...
				       rician,B1bias);toc()   
				       filtername='ODCT2';

 otherwise
display('Processing Real denoised image - MRONLM')
tic(),MRIdenoised = MRIDenoisingMRONLM(img,sigma,...
				       beta_, patcharea, searcharea, ...
				       rician,B1bias);toc()
				       filtername='MRONLM';
end

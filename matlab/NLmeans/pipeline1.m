function [MRIdenoised,sigma, filtername] = pipeline1(img,NLfilter,hfinal,hfactor,searcharea,patcharea,rician)
%% Non-local means denoising Option 1
%  Option 1 calls the automatic noise estimate before running the
%  NLmeans filter
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
% -     Monash Biomedical Imaging

%% parse inputs

ima1 = NormaliseImage2(abs(img))*256.0;

if nargin < 3 || isempty(hfinal)
   [hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(ima1,1,0)
end
if nargin < 4 || isempty(hfactor)
    hfactor=100;
end
hfinal = hfinal * (hfactor/100)

if nargin < 5 || isempty(searcharea)
    searcharea=3;
end
if nargin < 6 || isempty(patcharea)
    patcharea=1;
end
if nargin < 7 || isempty(rician)
    rician=1;
end
beta_=1;

%% run filter
display(['Noise estimate: ' num2str(hfinal)])

sigma=hfinal;
filtername='';

switch NLfilter
 case 0 
 display('Processing denoised image - MRONLM')
 tic(),MRIdenoised = MRIDenoisingMRONLM(ima1, sigma, beta_,patcharea, searcharea, rician, 0); toc()
					filtername='MRONLM';  
 case 1
 display('Processing denoised image - PRINLM')
 tic(),MRIdenoised = MRIDenoisingPRINLM(ima1, sigma, beta_,...
					rician, 0);toc()
					filtername='PRINLM';
 case 2
 display('Processing denoised image - AONLM')
 tic(),MRIdenoised = MRIDenoisingAONLM(ima1, beta_, patcharea, ...
				       searcharea, rician, 0);toc()
				       filtername='AONLM';
 case 3
 display('Processing denoised image - ONLM ')
 tic(),MRIdenoised = MRIDenoisingONLM(ima1, sigma, ...
				      beta_, patcharea, searcharea, ...
				      rician , 0);toc()  
				      filtername='ONLM';
 case 4
 display('Processing denoised image - ODCT ')
 tic(),MRIdenoised = MRIDenoisingODCT(ima1, ...
				      sigma, ...
				      beta_,rician,0);toc()   
				      filtername='ODCT';  
 case -5 
 display('Processing denoised image - MRONLM2')
 tic(),MRIdenoised = MRIDenoisingMRONLM2(ima1, sigma, ...
					 patcharea, searcharea, ...
					 rician, 0);toc() 
					 filtername='MRONLM2';
 case -1
 display('Processing denoised image - PRINLM2')
 tic(),MRIdenoised = MRIDenoisingPRINLM2(ima1, sigma,  ...
					 rician, 0);toc()
					 filtername='PRINLM2';
 case -2
 display('Processing denoised image - AONLM2')
 tic(),MRIdenoised = MRIDenoisingAONLM2(ima1, patcharea, ...
					searcharea, rician, 0);toc()
					filtername='AONLM2';
 case -3
 display('Processing denoised image - ONLM2 ')
 tic(),MRIdenoised = MRIDenoisingONLM2(ima1, sigma, ...
				       patcharea, searcharea, ...
				       rician , 0);toc()  
				       filtername='ONLM2';
 case -4
 display('Processing denoised image - ODCT2 ')
 tic(),MRIdenoised = MRIDenoisingODCT2(ima1, ...
				       sigma, ...
				       rician,0);toc()   
				       filtername='ODCT2';

 otherwise
display('Processing Real denoised image - MRONLM')
tic(),MRIdenoised = MRIDenoisingMRONLM(ima1,sigma,...
				       beta_, patcharea, searcharea, ...
				       rician,0);toc()
				       filtername='MRONLM';
end

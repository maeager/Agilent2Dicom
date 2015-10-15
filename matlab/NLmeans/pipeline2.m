function [MRIdenoised,sigma,filtername] = pipeline2(img1,img2,NLfilter,hfinal,hfactor,searcharea,patcharea,rician)
%% Non-local means denoising Option 2
% this method calculates the noise estimate from two images 
% and applies NL means to average image
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
% -     Monash Biomedical Imaging


if nargin<3
    NLfilter = 0;
end

%ima1=NormaliseImage2(abs(img1))*256.0;
%ima2=NormaliseImage2(abs(img2))*256.0;

avg=(img1+img2)/2;
if nargin < 4 || isempty(hfinal)
  diffimg=(img1-img2);
				%est_mean = mean(diffimg(:));
  est_std = std(diffimg(:))/sqrt(2);
  hfinal = est_std/sqrt(2);
else
  est_std = hfinal*sqrt(2);
end
if nargin < 5 || isempty(hfactor)
    hfactor=100;
end
hfinal = hfinal * (hfactor/100)

if nargin < 6 || isempty(searcharea)
  if NLfilter == -5 || NLfilter == 0 
    searcharea=3;
  else
    searcharea=5;
  end
end
if nargin < 7 || isempty(patcharea)
    patcharea=1;
end
if nargin < 8 || isempty(rician)
    rician=1;
end
beta_=1;





display(['Noise estimate: ' num2str(est_std)])
display(['Noise estimate of average image: ' num2str(hfinal)])

sigma=hfinal;
filtername='';

switch NLfilter
  case 0
    display('Processing denoised image - MRONLM')
    tic(),MRIdenoised = MRIDenoisingMRONLM(avg, hfinal, beta_, ...
                                           patcharea, searcharea, ...
                                           rician, 0); toc()
					   filtername='MRONLM';  
  case 1
    display('Processing denoised image - PRINLM')
    tic(),MRIdenoised = MRIDenoisingPRINLM(avg, hfinal, beta_, ...
                                           rician, 0);toc()
					   filtername='PRINLM';
  case 2
    display('Processing denoised image - AONLM')
    tic(),MRIdenoised = MRIDenoisingAONLM(avg, beta_, patcharea, ...
                                          searcharea, rician, 0);toc()
					  filtername='AONLM';
  case 3
    display('Processing denoised image - ONLM ')
    tic(),MRIdenoised = MRIDenoisingONLM(avg, hfinal, beta_, ...
                                         patcharea, searcharea, ...
                                         rician, 0);toc()
					 filtername='ONLM';
  case 4
    display('Processing denoised image - ODCT ')
    tic(),MRIdenoised = MRIDenoisingODCT(avg, hfinal, ...
                                         beta_, rician, 0);toc()
					 filtername='ODCT';  
  case -5
    display('Processing denoised image - MRONLM2')
    tic(),MRIdenoised = MRIDenoisingMRONLM2(avg, hfinal, ...
                                           patcharea, searcharea, ...
                                           rician);toc()
					   filtername='MRONLM2';
      
  case -1
    display('Processing denoised image - PRINLM2')
    tic(),MRIdenoised = MRIDenoisingPRINLM2(avg, hfinal,rician);toc()
      filtername='PRINLM2';
  case -2
    display('Processing denoised image - AONLM2')
    tic(),MRIdenoised = MRIDenoisingAONLM2(avg, hfinal, patcharea,searcharea, rician);toc()
					  filtername='AONLM2';
  case -3
    display('Processing denoised image - ONLM2 ')
    tic(),MRIdenoised = MRIDenoisingONLM2(avg, hfinal, ...
                                         patcharea, searcharea, ...
                                         rician, 0);toc()
					 filtername='ONLM2'
  case -4
    display('Processing denoised image - ODCT2 ')
    tic(),MRIdenoised = MRIDenoisingODCT2(avg, hfinal, ...
                                          rician, 0);toc()  
					 filtername='ODCT2';
  otherwise
    display('Processing denoised image - MRONLM')
    tic(),MRIdenoised = MRIDenoisingMRONLM(avg, hfinal, beta_,...
                                           patacharea, searcharea, rician, 0);toc()
					   filtername='MRONLM';
end



function MRIdenoised = pipeline2(img1,img2,NLfilter,hfinal,hfactor,searcharea,patacharea,rician)
%% Non-local means denoising Option 2
% this method calculates the noise estimate from two images 
% and applies NL means to average image
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
% -     Monash Biomedical Imaging


if nargin<3
    NLfilter = 0;
end

ima1=NormaliseImage2(abs(img1))*256;
ima2=NormaliseImage2(abs(img2))*256;

avg=(ima1+ima2)/2;
if nargin < 4 || isempty(hfinal)
diffima=(ima1-ima2);
%est_mean = mean(diffima(:));
est_std = std(diffima(:))/sqrt(2);
hfinal = est_std/sqrt(2);
else
est_std = hfinal*sqrt(2);
end
if nargin < 5 || isempty(hfactor)
    hfactor=100;
end
hfinal = hfinal * (hfactor/100)

if nargin < 6 || isempty(searcharea)
    searcharea=3;
end
if nargin < 7 || isempty(patcharea)
    patcharea=1;
end
if nargin < 8 || isempty(rician)
    rician=1;
end
beta=1;





display(['Noise estimate: ' num2str(est_std)])
display(['Noise estimate of average image: ' num2str(hfinal)])



switch NLfilter
  case 0
    display('Processing denoised image - MRONLM')
    tic(),MRIdenoised = MRIDenoisingMRONLM(avg, hfinal, beta, ...
                                           patcharea, searcharea, ...
                                           rician, 0);toc()
      
  case 1
    display('Processing denoised image - PRINLM')
    tic(),MRIdenoised = MRIDenoisingPRINLM(avg, hfinal, beta, ...
                                           rician, 0);toc()
      
  case 2
    display('Processing denoised image - AONLM')
    tic(),MRIdenoised = MRIDenoisingAONLM(avg, beta, patcharea, ...
                                          searcharea, rician, 0);toc()
    
  case 3
    display('Processing denoised image - ONLM ')
    tic(),MRIdenoised = MRIDenoisingONLM(avg, hfinal, beta, ...
                                         patcharea, searcharea, ...
                                         rician, 0);toc()
  case 4
    display('Processing denoised image - ODCT ')
    tic(),MRIdenoised = MRIDenoisingODCT(avg, hfinal, ...
                                         beta, rician, 0);toc()  
    
  otherwise
    display('Processing denoised image - MRONLM')
    tic(),MRIdenoised = MRIDenoisingMRONLM(avg, hfinal, beta,...
                                           patacharea, searcharea, rician, 0);toc()
      
end



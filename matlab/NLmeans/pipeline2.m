function MRIdenoised = pipeline2(img1,img2,NLfilter)
%% Non-local means denoising Option 2
% this method calculates the noise estimate from two images 
% and applies NL means to average image
%
% - (C) Michael Eager 2015


ima1=NormaliseImage2(abs(img1))*256;
ima2=NormaliseImage2(abs(img2))*256;

avg=(ima1+ima2)/2;
diffima=(ima1-ima2);
%est_mean = mean(diffima(:));
est_std = std(diffima(:))/sqrt(2);

display(['Noise estimate: ' num2str(est_std)])
display(['Noise estimate of average image: ' num2str(est_std/sqrt(2))])



switch NLfilter
  case 0
    display('Processing denoised image - MRONLM')
    tic(),MRIdenoised = MRIDenoisingMRONLM(avg,est_std/sqrt(2),1,1,3,1,0);toc()
      
  case 1
    display('Processing  denoised image - PRINLM')
    tic(),MRIdenoised = MRIDenoisingPRINLM(avg,est_std/sqrt(2),1,1,0);toc()
      
  case 2
    display('Processing Real denoised image - AONLM')
    tic(),MRIdenoised = MRIDenoisingAONLM(avg,1,1,3,1,0);toc()
    
  case 3
    display('Processing  denoised image -ONLM ')
    tic(),MRIdenoised = MRIDenoisingONLM(avg,est_std/sqrt(2),1,1,3,1,0);toc()
      
  otherwise
    display('Processing denoised image - MRONLM')
    tic(),MRIdenoised = MRIDenoisingMRONLM(avg,est_std/sqrt(2),1,1,3,1,0);toc()
      
end



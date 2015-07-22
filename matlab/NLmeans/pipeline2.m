function MRIdenoised = pipeline2(img1,img2)
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

MRIdenoised = MRIDenoisingMRONLM(avg,est_std/sqrt(2),1,1,3,1,0);





function MRIdenoised = pipeline1(img)
%% Non-local means denoising Option 1
%
%
% - (C) Michael Eager 2015


%%  Option 1
% this method calls the automatic noise estimate before running the
% NLmeans filter

ima1=NormaliseImage2(abs(img))*256;
[hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(ima1,1,1)
display(['Noise estimate: ' num2str(hfinal)])
MRIdenoised = MRIDenoisingMRONLM(ima1,hfinal,1,1,3,1,0);





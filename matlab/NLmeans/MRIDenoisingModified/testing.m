%% Testing MRI denoising

root_path='/home/vnmr1/src/Agilent2Dicom.git/matlab'
addpath(fullfile(root_path,'../matlab'))
addpath(fullfile(root_path,'../matlab/NIFTI'))
addpath(fullfile(root_path, '../matlab/Agilent/'))
addpath(fullfile(root_path, '../matlab/NLmeans'))
% add recursive directories in MRI denoise package
addpath(genpath(fullfile(root_path, '../matlab/NLmeans/MRIDenoisingPackage')))
addpath(fullfile(root_path, '../matlab/NLmeans/MRIDenoisingModified'))
run (fullfile(root_path,'../matlab/NLmeans/vlfeat/toolbox/vl_setup.m'))
img1 = readfid('~/NLmeans/New-joint-data/Knee1/ge3d-rep1_01.fid/');
img2 = readfid('~/NLmeans/New-joint-data/Knee1/ge3d-rep2_01.fid/');
cd matlab/NLmeans/MRIDenoisingModified/
mex myRINLM3d.cpp
mex -g myRINLM3d.cpp
coil=single(ones(size(img1)));
ima = single(NormaliseImage2(abs(img1)))*256.0;
[hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(ima,1,1)

tic;PRINL=myRINLM3d(ima,3,1,single(23.3298),single(ODCT),1,coil);toc



tic;ODCT=MRIDenoisingODCT2(ima,single(23.3298),1,1);toc
tic;ONLM=MRIDenoisingONLM2(ima,single(23.3298),1,3,1,1);toc
tic;AONLM=MRIDenoisingAONLM2(ima,single(23.3298),1,3,1,coil,1);toc
tic;MRNLM=MRIDenoisingMRONLM2(ima,single(23.3298),1,3,1,coil,1);toc
tic;imaPRINLM = MRIDenoisingPRINLM2(ima,single(23.3298),1,3,1,coil,1);toc



%% Pipeline 2
ima1=NormaliseImage2(abs(img1))*256.0;
ima2=NormaliseImage2(abs(img2))*256.0;
avg=(ima1+ima2)/2;
diffima=(ima1-ima2);
%est_mean = mean(diffima(:));
est_std = std(diffima(:))/sqrt(2);
hfinal = est_std/sqrt(2)
% 12.216
[hfinalOLD, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(avg,1,1)
% 12.5423
tic;ODCT=MRIDenoisingODCT2(avg,single(12.216),1,1);toc
tic;ONLM=MRIDenoisingONLM2(avg,single(12.216),1,3,1,1);toc
tic;AONLM=MRIDenoisingAONLM2(avg,single(12.216),1,3,1,coil,1);toc
tic;MRNLM=MRIDenoisingMRONLM2(avg,single(12.216),1,3,1,coil,1);toc
tic;imaPRINLM = MRIDenoisingPRINLM2(avg,single(12.216),1,3,1,coil,1);toc


%% Pipeline 3

[img1,hdr,ksp1] = readfid('~/NLmeans/New-joint-data/Knee1/ge3d-rep1_01.fid/');
[img2,hdr,ksp2] = readfid('~/NLmeans/New-joint-data/Knee1/ge3d-rep2_01.fid/');




[pha1, swi_n1, swi_p1, mag1] = phaserecon_v1(ksp1,ksp1,0.4,1,0.05);
hdyne1 = mag1.*exp(-1i*pha1);
[pha2, swi_n2, swi_p2, mag2] = phaserecon_v1(ksp2,ksp2,0.4,1,0.05);
hdyne2 = mag2.*exp(-1i*pha2);

hdyne1=flipdim(flipdim(flipdim(hdyne1,1),2),3);
hdyne2=flipdim(flipdim(flipdim(hdyne2,1),2),3);
hdyne1=circshift(hdyne1,[1,1,1]);
hdyne2=circshift(hdyne2,[1,1,1]);

%ima1=NormaliseImage2(abs(hdyne1))*256;
%ima2=NormaliseImage2(abs(hdyne2))*256;

rima1=real(hdyne1);
rima2=real(hdyne2);
iima1=imag(hdyne1);
iima2=imag(hdyne2);

avg_real=(rima1+rima2)/2;
avg_imag=(iima1+iima2)/2;
    % Calculate the noise estimate in the real image
    diffima_real=(rima1-rima2);
    est_std_real = std(diffima_real(:))/sqrt(2);
    % Calculate the noise estimate in the imaginary image
    diffima_imag=(iima1-iima2);
    est_std_imag = std(diffima_imag(:))/sqrt(2);
    hfinal_real = est_std_real/sqrt(2);
    hfinal_imag = est_std_imag/sqrt(2);
    hfactor=100;
hfinal_real = hfinal_real * (hfactor/100)
hfinal_imag = hfinal_imag * (hfactor/100)
    
rician=1;beta_=1;searcharea=3;patcharea=1

    display('Processing Real denoised image - MRONLM2')
    tic(),MRIdenoised_real = MRIDenoisingMRONLM2(avg_real, hfinal_real,beta_, ...
                                           patcharea, searcharea, ...
                                           rician);toc()
    display('Processing Imag denoised image - MRONLM2')
    MRIdenoised_imag = MRIDenoisingMRONLM2(avg_imag, hfinal_imag, beta_, ...
                                           patcharea, searcharea, ...
                                           rician);toc()
  
    display('Processing Real denoised image - PRINLM2')
    tic(),MRIdenoised_real = MRIDenoisingPRINLM2(avg_real,hfinal_real, ...
                                           rician);toc()
    display('Processing Imag denoised image - PRINLM2')
    MRIdenoised_imag = MRIDenoisingPRINLM2(avg_imag,hfinal_imag, ...
                                           rician);toc()
   
    display('Processing Real denoised image - AONLM2')
    tic(),MRIdenoised_real = MRIDenoisingAONLM2(avg_real, hfinal, patcharea,    searcharea, rician);toc()
    display('Processing Imag denoised image - AONLM2 ')
    MRIdenoised_imag = MRIDenoisingAONLM2(avg_imag, hfinal, patcharea,  searcharea, rician);toc()
   
    display('Processing Real denoised image -ONLM2 ')
    tic(),MRIdenoised_real = MRIDenoisingONLM2(avg_real,hfinal_real, patcharea, searcharea,rician);toc()
    display('Processing Imag denoised image -ONLM2 ')
    MRIdenoised_imag = MRIDenoisinONLM2(avg_imag, hfinal_imag,  patcharea, searcharea, rician);toc()
  
    display('Processing Real denoised image - ODCT2')
    tic(),MRIdenoised_real = MRIDenoisingODCT2(avg_real,hfinal_real,  rician);toc()
    display('Processing Imag denoised image - ODCT2')
    MRIdenoised_imag = MRIDenoisingODCT2(avg_imag,hfinal_imag,  rician); toc()
                         
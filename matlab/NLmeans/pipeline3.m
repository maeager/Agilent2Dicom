function MRIdenoised = pipeline3(ksp1,ksp2)
%%  Non-local means denoising Option 3
% this method calculates the noise estimate from two images 
% and applies NL means to the complex average image
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
% Monash Biomedical Imaging


[pha1, swi_n1, swi_p1, mag1] = phaserecon_v1(ksp1,ksp1,0.4,1,0.05);
hdyne1 = mag1.*exp(-1i*pha1);
[pha2, swi_n2, swi_p2, mag2] = phaserecon_v1(ksp2,ksp2,0.4,1,0.05);
hdyne2 = mag2.*exp(-1i*pha2);

% Necessary translations to match FDF images
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

display(['Noise estimate (real and imag): ' num2str(est_std_real) ' ' num2str(est_std_imag)])
display(['Noise estimate of 2-average image: ' num2str(est_std_real/sqrt(2)) ' ' num2str(est_std_imag/sqrt(2))])

if est_std_real < eps(single(1))*1000 || est_std_imag < ...
        eps(single(1))*1000
    display(['Noise estimate too small.  Make sure the two images ' ...
             'are not the same.'])
    return
end

display('Processing Real denoised image')
tic(),MRIdenoised_real = MRIDenoisingMRONLM(avg_real,est_std_real/sqrt(2),1,1,3,1,0);toc()
display('Processing Imag denoised image')
MRIdenoised_imag = MRIDenoisingMRONLM(avg_imag,est_std_imag/sqrt(2),1,1,3,1,0);toc()
MRIdenoised = complex(MRIdenoised_real,MRIdenoised_imag);





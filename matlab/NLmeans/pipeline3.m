function MRIdenoised = pipeline3(ksp1,ksp2,NLfilter,hfinal,hfactor,searcharea,patacharea,rician)
%%  Non-local means denoising Option 3
% this method calculates the noise estimate from two images 
% and applies NL means to the complex average image
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
% Monash Biomedical Imaging

if nargin<3
    NLfilter = 0;
end


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

if nargin < 4 || isempty(hfinal)
    % Calculate the noise estimate in the real image
    diffima_real=(rima1-rima2);
    est_std_real = std(diffima_real(:))/sqrt(2);
    % Calculate the noise estimate in the imaginary image
    diffima_imag=(iima1-iima2);
    est_std_imag = std(diffima_imag(:))/sqrt(2);
    
    
    if est_std_real < eps(single(1))*1000 || est_std_imag < ...
            eps(single(1))*1000
        display(['Noise estimate too small.  Make sure the two images ' ...
                 'are not the same.'])
        return
    end
    
else
    display(['pipeline 3 cannot be used with user derived hfinal. ' ...
             'Terminating process.'])
    return
end
if nargin < 5 || isempty(hfactor)
    hfactor=100;
end
hfinal_real = hfinal_real * (hfactor/100)
hfinal_imag = hfinal_imag * (hfactor/100)

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




display(['Noise estimate (real and imag): ' num2str(est_std_real) ' ' num2str(est_std_imag)])
display(['Noise estimate of 2-average image: ' num2str(hfinal_real) ' ' num2str(hfinal_imag)])




switch NLfilter
  case 0
    display('Processing Real denoised image - MRONLM')
    tic(),MRIdenoised_real = MRIDenoisingMRONLM(avg_real, hfinal_real,beta, ...
                                           patcharea, searcharea, ...
                                           rician, 0);toc()
    display('Processing Imag denoised image - MRONLM')
    MRIdenoised_imag = MRIDenoisingMRONLM(avg_imag, hfinal_imag, beta, ...
                                           patcharea, searcharea, ...
                                           rician, 0);toc()
  case 1
    display('Processing Real denoised image - PRINLM')
    tic(),MRIdenoised_real = MRIDenoisingPRINLM(avg_real,hfinal_real, beta, ...
                                           rician, 0);toc()
    display('Processing Imag denoised image - PRINLM')
    MRIdenoised_imag = MRIDenoisingPRINLM(avg_imag,hfinal_imag, beta, ...
                                           rician, 0);toc()
  case 2
    display('Processing Real denoised image - AONLM')
    tic(),MRIdenoised_real = MRIDenoisingAONLM(avg_real,beta, patcharea, ...
                                          searcharea, rician, 0);toc()
    display('Processing Imag denoised image - AONLM ')
    MRIdenoised_imag = MRIDenoisingAONLM(avg_imag,beta, patcharea, ...
                                          searcharea, rician, 0);toc()
  case 3
    display('Processing Real denoised image -ONLM ')
    tic(),MRIdenoised_real = MRIDenoisingONLM(avg_real,hfinal_real, beta, ...
                                         patcharea, searcharea, ...
                                         rician, 0);toc()
    display('Processing Imag denoised image -ONLM ')
    MRIdenoised_imag = MRIDenoisinONLM(avg_imag, hfinal_imag, beta, ...
                                         patcharea, searcharea, ...
                                         rician, 0);toc()
  case 4
    display('Processing Real denoised image - ODCT')
    tic(),MRIdenoised_real = MRIDenoisingODCT(avg_real,hfinal_real,...
                                         beta, rician, 0);toc()
    display('Processing Imag denoised image - ODCT')
    MRIdenoised_imag = MRIDenoisingODCT(avg_imag,hfinal_imag,...
                                         beta, rician, 0);toc()
  otherwise
    display('Processing Real denoised image - MRONLM')
    tic(),MRIdenoised_real = MRIDenoisingMRONLM(avg_real,hfinal_real,beta,...
                                           patacharea, searcharea, rician, 0);toc()
    display('Processing Imag denoised image - MRONLM')
    MRIdenoised_imag = MRIDenoisingMRONLM(avg_imag,hfinal_imag,beta,...
                                           patacharea, searcharea, rician, 0);toc()
end

MRIdenoised = complex(MRIdenoised_real,MRIdenoised_imag);





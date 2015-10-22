function [pha, swi_n, swi_p, mag] = phaserecon2D(kimg,kimgsos,a,intpl,thr)
% Usage: calculate phase images from complex inputs
% Inputs:
%        kimg -- reconstructed kspace signals
%        kimgsos -- sos of square on complex images
%        a    -- gaussian filter size 0 to 20, 10 is good balance
%        intpl -- interplation factor
%        thr -- thresholding: 0 no thresholding
%                             0.05 a good choice
% Outputs
%        pha -- phase images
%        mag -- magnatitue images
%        swi -- phase weighted magnatitue images
%
% Zhaolin Chen @ Howard Florey Institute
% log            12-05-09 SWI part fixed
%
% Amendents:
% Michael Eager Monash Biomedical Imaging
% 09-2015 2D implementation of homodyne filter
% -------------------------------------------


[Nfe,Npe,Npe2] = size(kimg);

if isempty(kimgsos)
    kimgsos=kimg;
end


% create a Gaussian LPF
win = gausswin(Nfe,a)*gausswin(Npe,a)';
%win3=gausswin(Npe2,a);
%win=zeros(Nfe,Npe,Npe2);
%for k=1:Npe2
%    win(:,:,k)=win2*win3(k);
%end
win=win./(max(win(:)));

%win=fspecial3('gaussian',[Nfe,Npe,Npe2],a);
% creat a rectangular LPF
% win = zeros(Nfe,Npe);
% L = 32;W = 32;
% win(round(Nfe/2-L/2):round(Nfe/2+L/2-1),round(Npe/2-W/2):round(Npe/2+W/2-1)) = ones(L,W);


% creat a Kaiser LPF
%win = hann(Nfe)*hann(Npe)';
%win = kaiser(Nfe,a)*kaiser(Npe,a)';
if numel(size(kimg)==5)
    kimg2=kimg(:,:,:,1,1);
    [xx,yy,zz] = find(kimg2 == max(kimg2(:)));
else
   [xx,yy,zz] = find(kimg == max(kimg(:)));
end
win = circshift(win,[xx,yy]-[floor(Nfe/2), floor(Npe/2)]);

img = fftshift(fftn(kimg));%,[intpl*Nfe,intpl*Npe,intpl*Npe2]));
imgsos = fftshift(fftn(kimgsos));%,[intpl*Nfe,intpl*Npe,intpl*Npe2]));

if numel(size(kimg)==5)
    for echo=1:size(kimg,5)
        for channel=1:size(kimg,4)
            for slice=1:size(kimg,3)
                img_lpf(:,:,slice,channel,echo) = fftshift(fft2(kimg(:,:,slice,channel,echo).*win));
            end
        end
    end
else
    for slice=1:size(kimg,3)
        img_lpf(:,:,slice) = fftshift(fft2(kimg(:,:,slice).*win));
    end
end

img_hpf = img ./ img_lpf;

mag = abs(img);

thd = (thr/sqrt(sqrt(intpl)))*max(abs(imgsos(:)));


pha = angle(img_hpf);


if (thr ~= 0)
    %    for i = 1:intpl*Nfe
    %        for j = 1:intpl*Npe
    %            if (abs(imgsos(i,j)) <= thd)
    %               pha(i,j) = pha(i,j)/10;
    %            end
    %        end
    %    end
    x = find(abs(imgsos(:)) <= thd);
    pha(x)=pha(x)/10;
end

%mag = sqrt(mag);


phasemask_n = (pha+pi)./pi;

[x] = find(pha>=0);
phasemask_n(x) = 1;

phasemask_p = (pi-pha)./pi;

[x] = find(pha<=0);
phasemask_p(x) = 1;

swi_n = phasemask_n.^4.*mag;
swi_p = phasemask_p.^4.*mag;


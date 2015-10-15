function call_pipeline2(in1,in2,out,NLfilter,hfinal,hfactor,searcharea,patcharea,rician)
% Calling non-local means pipeline 2
%
% - (C) 2015 Michael Eager (michael.eager@monash.edu)
% - Monash Biomedical Imaging

[a,b,c] = fileparts(mfilename('fullpath')) ;
[a,b,c] = fileparts(a) ;
root_path=a;
warning off MATLAB:dispatcher:nameConflict
addpath(fullfile(root_path,'../matlab'))
addpath(fullfile(root_path,'../matlab/NIFTI'))
addpath(fullfile(root_path, '../matlab/Agilent/'))
addpath(fullfile(root_path, '../matlab/NLmeans'))
% add recursive directories in MRI denoise package
addpath(genpath(fullfile(root_path, '../matlab/NLmeans/MRIDenoisingPackage')))
addpath(fullfile(root_path, ['../matlab/NLmeans/' ...
    'MRIDenoisingModified']))
run (fullfile(root_path, '../matlab/NLmeans/vlfeat/toolbox/vl_setup.m'))

display('Calling non-local means filter pipeline 1')
if nargin < 4
    NLfilter=0;
end
if nargin < 5
    hfinal=[];hfactor=[];rician=[];
    searcharea=[];patcharea=[];
end
%% Clean input strings
in1 = regexprep(in1,'["\[\]]','');
out = regexprep(out,'["\[\]]','');
in2 = regexprep(in2,'["\[\]]','');


if exist(in1,'file')==2 && ~isempty(strfind(in1,'.nii'))
    nii1_in=load_nii(in1);
    img1=nii1_in.img;
    voxelsize1 = nii1_in.hdr.dime.pixdim(2:4);
elseif ~isempty(strfind(in1,'.img')) && isdir(in1)
    [img1 hdr] =readfdf(in1);
    %    voxelsize=hdr.FOVcm/size(img)*10;
    %     voxelsize1 = hdr.roi*10/hdr.matrix;
    voxelsize1 = hdr.voxelsize*10;
elseif ~isempty(strfind(in1,'.fid')) && isdir(in1)
    [img1, hdr, ksp, RE, IM] = readfid(in1);
    %    voxelsize1=hdr.FOVcm*10/size(img);
    voxelsize1=hdr.voxelmm;
else
    display(['Cannot find ' in1])
    return
end

img1=squeeze(img1);



if exist(in2,'file') && ~isempty(strfind(in2,'.nii'))
    nii2_in=load_nii(in2);
    img2=nii2_in.img;
    voxelsize2 = nii2_in.hdr.dime.pixdim(2:4);
elseif ~isempty(strfind(in2,'.img')) && isdir(in2)
    [img2 hdr] =readfdf(in2);
    %    voxelsize=hdr.FOVcm/size(img)*10;
    voxelsize2 = hdr.voxelsize*10;
elseif ~isempty(strfind(in2,'.fid')) && isdir(in2)
    [img2, hdr, ksp, RE, IM] = readfid(in2);
    %    voxelsize2=hdr.FOVcm*10/size(img);
    voxelsize2 = hdr.voxelmm;
else
    display(['Cannot find ' in2])
    return
end

img2=squeeze(img2);


if sum(voxelsize1) ~= sum(voxelsize2)
    display(['Voxelsizes don''t match: ' str2num(voxelsize1) ' ' ...
        str2num(voxelsize2)])
    return
end
voxelsize=voxelsize1;

% Normalise images

[img1,Range] = NormaliseImage2(squeeze(img1));
img1 = img1*256.0;
img2 = (img2 - Range(1))/(Range(2)-Range(1));
min_img2=min(img2(:));
if min_img2 < 0
  img2=img2-min_img2;
  img1=img1-min_img2;
end


if length(size(img)) == 3
    
    display('Calling pipeline 2 on 3D volume')
    tic(),[MRIdenoised2,sigma,filtername] = pipeline2(img1, img2, NLfilter,...
        hfinal, hfactor, searcharea, patcharea, ...
        rician);toc()
    
elseif length(size(img)) == 4
    
    display 'Processing 4D image'
    %    img=img(:,:,:,1);
    for vol=1:size(img,4)
        display(['Calling pipeline 2 on 4D volume ' num2str(vol)])
        tic()
        [MRIdenoised2(:,:,:,vol), sigma, filtername] = pipeline2(img1(:,:,:,vol), img2(:,:,:,vol), NLfilter,...
            hfinal, hfactor, searcharea, patcharea, ...
            rician);
        toc()
    end
    
elseif length(size(img)) == 5
    display 'Processing 5D image'
    %    img=img(:,:,:,1,1);
    for echo=1:size(img,5)
        for vol=1:size(img,4)
            display(['Calling pipeline 2 on 5Dvolume ' num2str(vol) ' echo ' num2str(echo)])
            tic()
            [MRIdenoised2(:,:,:,vol,echo),sigma,filtername] = pipeline2(img1(:,:,:,vol,echo), img2(:,:,:,vol,echo), NLfilter,...
                hfinal, hfactor, searcharea, patcharea, rician);
            toc()
        end
    end
    
elseif length(size(img)) > 5
    display 'Unable to process images greater than 5D'
    return
end




%% Save image to nifti

filter_line = ['filter' filtername '_sigma' num2str(sigma) ];


if exist(out,'file') == 2
    delete(out)
    denoised_file = [ out(1:end-8) '_' filter_line out(end-7:end)];
    [a,b,c] = fileparts(out) ;
    raw_file = [ a, '/raw_average.nii.gz'];
else
    if isdir(out)
        denoised_file = [out, '/pipeline2_' filter_line '.nii.gz'];
        raw_file = [out, '/raw_average.nii.gz'];
    else
        mkdir(out)
        denoised_file = [out, '/pipeline2_' filter_line '.nii.gz'];
        raw_file = [out, '/raw_average.nii.gz'];
    end
end
if exist(raw_file,'file') ~= 2
    display(['Saving ' raw_file ])
    save_nii(make_nii(abs(img1+img2)/2,voxelsize,[],16),raw_file)
    % display(['Deleting ' raw_file ])
    % delete(raw_file)
end
if exist(denoised_file,'file')
    display(['Deleting old ' denoised_file ])
    delete(denoised_file)
end

display(['Saving ' denoised_file ])
save_nii(make_nii(MRIdenoised2,voxelsize,[],16),denoised_file)

function call_pipeline1(in,out,in2,NLfilter,hfinal,hfactor,searcharea,patcharea,rician)
% Calling non-local means pipeline 1
%  Use rician noise estimator to calculate std of noise in one MRI
%  magnitude image
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
% - Monash Biomedical Imaging

display([' In call_pipeline1 ' in  out])
[a,b,c] = fileparts(mfilename('fullpath'))
[a,b,c] = fileparts(a) ;
root_path=a
warning off MATLAB:dispatcher:nameConflict
addpath(fullfile(root_path,'../matlab'))
addpath(fullfile(root_path,'../matlab/NIFTI'))
addpath(fullfile(root_path, '../matlab/Agilent/'))
addpath(fullfile(root_path, '../matlab/NLmeans'))
% add recursive directories in MRI denoise package
addpath(genpath(fullfile(root_path, '../matlab/NLmeans/MRIDenoisingPackage')))
addpath(fullfile(root_path, '../matlab/NLmeans/MRIDenoisingModified'))

run (fullfile(root_path,'../matlab/NLmeans/vlfeat/toolbox/vl_setup.m'))

display('Calling non-local means filter pipeline 1')

voxelsize=[];
if nargin < 4
    hfinal=[];hfactor=[];rician=[];
    searcharea=[];patcharea=[];
end
if nargin < 3
    NLfilter=0;
end
%% Clean input strings
in = regexprep(in,'"','');
out = regexprep(out,'"','');
if nargin >= 3
    if isstr(in2)
        in2 = regexprep(in2,'"','');
    else
        in2=[];
    end
end



if exist(in,'file')==2 && ~isempty(strfind(in,'.nii'))
    nii_in=load_nii(in);
    img = nii_in.img;
    voxelsize = nii_in.hdr.dime.pixdim(2:4);
elseif ~isempty(strfind(in,'.img')) && isdir(in)
    [img hdr] =readfdf(in);
    %    voxelsize=hdr.FOVcm/size(img)*10;
    % voxelsize=hdr.roi*10/hdr.matrix;
    voxelsize = hdr.voxelsize * 10;
elseif ~isempty(strfind(in,'.fid')) && isdir(in)
    [img, hdr, ksp, RE, IM] = readfid(in);
    % voxelsize = hdr.FOVcm * 10 / size(img)
    voxelsize = hdr.voxelmm;
    img=abs(img);
else
    display(['Cannot find ' in])
    return
end

if nargin >= 3 && ~isempty(in2)
    display('Pipeline 1 with real and imag')
    if exist(in2,'file') && ~isempty(strfind(in2,'.nii'))
        nii2_in=load_nii(in2);
        img2=nii2_in.img;
        voxelsize2 = nii2_in.dime.pixdim(2:4);
    elseif ~isempty(strfind(in2,'.img')) && isdir(in2)
        [img2 hdr] =readfdf(in2);
        %    voxelsize=hdr.FOVcm/size(img)*10;
        voxelsize2 = hdr.voxelsize * 10;
    elseif ~isempty(strfind(in2,'.fid')) && isdir(in2)
        [img2, hdr, ksp, RE, IM] = readfid(in2);
        %    voxelsize2=hdr.FOVcm*10/size(img);
        voxelsize2 = hdr.voxelmm;
    else
        display(['Cannot find ' in2])
        return
    end
    img = abs(complex(img,img2))
end

% Squeeze and normalise image
img = NormaliseImage2(squeeze(img))*256.0;


if length(size(img)) == 3
    display 'Processing 3D image'
    display(['Calling pipeline 1'])
    tic()
    [MRIdenoised1,sigma,filtername] = pipeline1(img, NLfilter,...
        hfinal, hfactor, searcharea, patcharea, rician);
    toc()


elseif length(size(img)) == 4
    display 'Processing 4D image'
    for vol=1:size(img,4)
        display(['Calling pipeline 1 on volume ' num2str(vol)])
        tic()
        [MRIdenoised1(:,:,:,vol),sigma,filtername] = pipeline1(img(:,:,:,vol), NLfilter,...
            hfinal, hfactor, searcharea, patcharea, rician);
        toc()
    end
    
    
elseif length(size(img)) == 5
        display 'Processing 5D image'
        for echo=1:size(img,5)
            for vol=1:size(img,4)
                display(['Calling pipeline 1 on volume ' num2str(vol) ' echo ' num2str(echo)])
                tic()
                [MRIdenoised1(:,:,:,vol,echo),sigma,filtername] = pipeline1(NormaliseImage2(img(:,:,:,vol,echo))*256, NLfilter,...
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

if exist(out,'file')==2
    delete(out)
    denoised_file = [ out(1:end-8) '_' filter_line out(end-7:end)];
    [a,b,c] = fileparts(out) ;
    raw_file = [ a, '/raw_average.nii.gz'];
else
    if ~isdir(out)
        mkdir(out)
    end
    denoised_file = [out '/pipeline1_' filter_line '.nii.gz'];
    raw_file = [out, '/raw_average.nii.gz'];
end
% Save raw input image average if not already saved
if exist(raw_file,'file') ~= 2
    save_nii(make_nii(img,voxelsize,[],16),raw_file)
end
% Save the denoised image
save_nii(make_nii(MRIdenoised1,voxelsize,[],16),denoised_file)

function call_pipeline3(in1,in2,out,flags,NLfilter)
% Calling non-local means pipeline 3
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
run  (fullfile(root_path, '../matlab/NLmeans/vlfeat/toolbox/vl_setup.m'))

display('Calling non-local means filter pipeline 3')
saveRI=0;savePhase=0;
if nargin <4
    flags=0
end
if flags>=4
    savePhase=1;
end
if mod(flags,2)==1
    saveRI=1;
end
if nargin < 5
    NLfilter=0;
end

%% Clean input strings
in1 = regexprep(in1,'"','');
in2 = regexprep(in2,'"','');
out = regexprep(out,'"','');

voxelsize=[];
shifted=0;
if exist(in1,'file')==2 && ~isempty(strfind(in1,'.nii')) 
    nii1_in=load_nii(in1);
    img1=nii1_in.img;
    ksp1=fftshift(fftn(img1));shifted=1;
    voxelsize1=nii1_in.hdr.dime.pixdim(2:4);
elseif ~isempty(strfind(in1,'.img')) && isdir(in1)
    [img1 hdr] =readfdf(in1);
    %    voxelsize=hdr.FOVcm/size(img)*10;
    ksp1=fftshift(fftn(img1));shifted=1;
    %    voxelsize1=hdr.roi*10/hdr.matrix;
    voxelsize1 = hdr.voxelsize*10;
elseif ~isempty(strfind(in1,'.fid')) && isdir(in1)
    [img1, hdr, ksp1, RE, IM] = readfid(in1);
    voxelsize1=hdr.voxelmm;
    %    voxelsize=hdr.FOVcm*10/size(img);
else
    display(['Cannot find ' in1])
    return
end


ksp1=squeeze(ksp1);img1=squeeze(img1);

if length(size(ksp1)) == 4
    display 'Reducing 4D image down to 3D'
  ksp1=ksp1(:,:,:,1);img1=img1(:,:,:,1);
elseif length(size(ksp1)) == 5
    display 'Reducing 5D image down to 3D'
    ksp1=ksp1(:,:,:,1,1);img1=img1(:,:,:,1,1);
elseif length(size(ksp1)) > 5
    display 'Unable to process images greater than 5D'
    return
end


if exist(in2,'file')==2 && ~isempty(strfind(in2,'.nii')) 
    nii2_in=load_nii(in2);
    img2 = nii2_in.img;
    ksp2 = fftshift(fftn(img2));shifted=1;
    voxelsize2=nii2_in.hdr.dime.pixdim(2:4);
elseif ~isempty(strfind(in2,'.img')) && isdir(in2)
    [img2 hdr] =readfdf(in2);
    %    voxelsize=hdr.FOVcm/size(img)*10;
    ksp2=fftshift(fftn(img2));shifted=1;
    %voxelsize2=hdr.roi*10/hdr.matrix;
    voxelsize2 = hdr.voxelsize*10;
elseif ~isempty(strfind(in2,'.fid')) && isdir(in2)
    [img2, hdr, ksp2, RE, IM] = readfid(in2);
    %    voxelsize2=hdr.FOVcm*10/size(img);
    voxelsize2=hdr.voxelmm;
else
    display(['Cannot find ' in2])
    return
end

ksp2=squeeze(ksp2);img2=squeeze(img2);

if length(size(ksp2)) == 4
    display 'Reducing 4D image down to 3D'
    ksp2=ksp2(:,:,:,1); img2=img2(:,:,:,1);
elseif length(size(ksp2)) == 5
    display 'Reducing 5D image down to 3D'
    ksp2=ksp2(:,:,:,1,1); img2=img2(:,:,:,1,1);
elseif length(size(ksp2)) > 5
    display 'Unable to process images greater than 5D'
    return
end



if sum(voxelsize1) ~= sum(voxelsize2)
  display(['Voxelsizes don''t match: ' str2num(voxelsize1) ' ' ...
           str2num(voxelsize2)])
  return
end
voxelsize=voxelsize1;

display 'Calling pipeline 3'
tic(),MRIdenoised3 = pipeline3(ksp1,ksp2,NLfilter);toc()
if shifted
    MRIdenoised3 = ifftshift(MRIdenoised3);
end

if exist(out,'file') == 2
    delete(out)
    denoised_file=out
    [a,b,c] = fileparts(out) ;
    raw_file = [ a, '/raw_average.nii.gz'];						   
else
    if ~isdir(out)
    % if not a file or a dir, create dir
    mkdir (out)
    end
    denoised_file=[out '/pipeline3_magn.nii.gz'];
    raw_file = [out, '/raw_average_magn.nii.gz'];

end

if exist(raw_file,'file')
    display(['Deleting ' raw_file ])
    delete(raw_file)
end
if exist(denoised_file,'file')
    display(['Deleting ' denoised_file ])
    delete(denoised_file)
end
display(['Saving ' raw_file ])
save_nii(make_nii(abs(img1+img2)/2,voxelsize,[],16),raw_file)
display(['Saving ' denoised_file ])
save_nii(make_nii(abs(MRIdenoised3),voxelsize,[],16),denoised_file)

if savePhase
    denoised_phase_file = regexprep(denoised_file,'magn','real');
    if exist(denoised_phase_file,'file')
      display(['Deleting ' denoised_phase_file ])
      delete(denoised_phase_file)
    end
    display(['Saving ' denoised_phase_file ])
    save_nii(make_nii(angle(MRIdenoised3),voxelsize,[],16),denoised_phase_file)
end
if saveRI
    denoised_real_file = regexprep(denoised_file,'magn','real');
    denoised_imag_file = regexprep(denoised_file,'magn','imag')
    if exist(denoised_real_file,'file')
       display(['Deleting ' denoised_real_file ])
       delete(denoised_real_file)
    end
    if exist(denoised_imag_file,'file')
       display(['Deleting ' denoised_imag_file ])
       delete(denoised_imag_file)
    end
    display(['Saving ' denoised_real_file ])
    save_nii(make_nii(real(MRIdenoised3),voxelsize,[],16),denoised_real_file)
    display(['Saving ' denoised_imag_file ])
    save_nii(make_nii(imag(MRIdenoised3),voxelsize,[],16),denoised_imag_file)
end    

function call_pipeline2(in1,in2,out)
% Calling non-local means pipeline 2
%
% - (C) 2015 Michael Eager (michael.eager@monash.edu)
% - Monash Biomedical Imaging

[a,b,c] = fileparts(mfilename('fullpath')) ;
[a,b,c] = fileparts(a) ;
root_path=a;
addpath(fullfile(root_path,'matlab'))
addpath(fullfile(root_path,'matlab/NIFTI'))
addpath(fullfile(root_path, 'matlab/Agilent/'))
addpath(fullfile(root_path, 'matlab/NLmeans'))
% add recursive directories in MRI denoise package
addpath(genpath(fullfile(root_path, 'matlab/NLmeans/MRIDenoisingPackage')))
run  matlab/NLmeans/vlfeat-0.9.17/toolbox/vl_setup.m

display('Calling non-local means filter pipeline 1')

if exist(in1,'file')==2 && ~isempty(strfind(in1,'.nii'))
    nii1_in=load_nii(in1);
    img1=nii1_in.img;
    voxelsize1 = nii1_in.dime.pixdim(2:4);
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

if exist(in2,'file') && ~isempty(strfind(in2,'.nii')) 
    nii2_in=load_nii(in2);
    img2=nii2_in.img;
    voxelsize2 = nii2_in.dime.pixdim(2:4);
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

if sum(voxelsize1) ~= sum(voxelsize2)
  display(['Voxelsizes don''t match: ' str2num(voxelsize1) ' ' ...
           str2num(voxelsize2)])
  return
end
voxelsize=voxelsize1;

display('Calling pipeline 2')
tic(),MRIdenoised2=pipeline2(img1,img2);toc()


if exist(out,'file')==2
    delete(out)
else
    if isdir(out)
        out=[out, '/pipeline2.nii.gz'];
    else
        mkdir(out)
        out=[out, '/pipeline2.nii.gz'];
    end
end
save_nii(make_nii(MRIdenoised2,voxelsize,[],16),out)

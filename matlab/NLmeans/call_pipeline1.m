function call_pipeline1(in,out)
% Calling non-local means pipeline 1
%
% - (C) Michael Eager 2015 (michael.eager@monash.edu)
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

voxelsize=[];

if exist(in,'file')==2 && ~isempty(strfind(in,'.nii'))
    nii_in=load_nii(in);
    img=nii_in.img;
    voxelsize=nii_in.dime.pixdim(2:4);
elseif ~isempty(strfind(in,'.img')) && isdir(in)
    [img hdr] =readfdf(in);
    %    voxelsize=hdr.FOVcm/size(img)*10;
    voxelsize=hdr.roi*10/hdr.matrix;
    voxelsize2 = hdr.voxelsize*10;
elseif ~isempty(strfind(in,'.fid')) && isdir(in)
    [img, hdr, ksp, RE, IM] = readfid(in);
    % voxelsize=hdr.FOVcm*10/size(img)
    voxelsize=hdr.voxelmm;
    img=abs(img);
else
    display(['Cannot find ' in])
    return
end

display 'Calling pipeline 1'
tic(),MRIdenoised1=pipeline1(NormaliseImage2(img)*256);toc()

if exist(out,'file')==2
    delete(out)
else
if isdir(out)
    out=[out '/pipeline1.nii.gz'];
else 
    mkdir(out)
    out=[out '/pipeline1.nii.gz'];
end
end
save_nii(make_nii(MRIdenoised1,voxelsize,[],16),out)

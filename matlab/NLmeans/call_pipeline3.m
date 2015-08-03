function call_pipeline3(in1,in2,out)
% Calling non-local means pipeline 3
%
% - (C) 2015 Michael Eager (michael.eager@monash.edu)
% - Monash Biomedical Imaging

[a,b,c] = fileparts(mfilename('fullpath')) ;
[a,b,c] = fileparts(a) ;
root_path=a;
addpath(fullfile(root_path,'../matlab'))
addpath(fullfile(root_path,'../matlab/NIFTI'))
addpath(fullfile(root_path, '../matlab/Agilent/'))
addpath(fullfile(root_path, '../matlab/NLmeans'))
% add recursive directories in MRI denoise package
addpath(genpath(fullfile(root_path, '../matlab/NLmeans/MRIDenoisingPackage')))
run  (fullfile(root_path, '../matlab/NLmeans/vlfeat/toolbox/vl_setup.m'))

display('Calling non-local means filter pipeline 3')

voxelsize=[];

if exist(in1,'file')==2 && ~isempty(strfind(in1,'.nii')) 
    nii1_in=load_nii(in1);
    img=nii1_in.img;
    ksp1=fftn(img);
    voxelsize1=nii1_in.dime.pixdim(2:4);
elseif ~isempty(strfind(in1,'.img')) && isdir(in1)
    [img hdr] =readfdf(in1);
    %    voxelsize=hdr.FOVcm/size(img)*10;
    ksp1=fftn(img);
    %    voxelsize1=hdr.roi*10/hdr.matrix;
voxelsize1 = hdr.voxelsize*10;
elseif ~isempty(strfind(in1,'.fid')) && isdir(in1)
    [img, hdr, ksp1, RE, IM] = readfid(in1);
    voxelsize1=hdr.voxelmm;
    %    voxelsize=hdr.FOVcm*10/size(img);
else
    display(['Cannot find ' in1])
    return
end

if exist(in2,'file')==2 && ~isempty(strfind(in2,'.nii')) 
    nii2_in=load_nii(in2);
    img=nii2_in.img;
    ksp2=fftn(img);
    voxelsize2=nii2_in.dime.pixdim(2:4);
elseif ~isempty(strfind(in2,'.img')) && isdir(in2)
    [img hdr] =readfdf(in2);
    %    voxelsize=hdr.FOVcm/size(img)*10;
    ksp2=fftn(img);
    %voxelsize2=hdr.roi*10/hdr.matrix;
    voxelsize2 = hdr.voxelsize*10;
elseif ~isempty(strfind(in2,'.fid')) && isdir(in2)
    [img, hdr, ksp2, RE, IM] = readfid(in2);
    %    voxelsize2=hdr.FOVcm*10/size(img);
    voxelsize2=hdr.voxelmm;
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

display 'Calling pipeline 3'
tic(),MRIdenoised3 = pipeline3(ksp1,ksp2);toc()

if exist(out,'file')~=2 && ~isdir(out)
    %if not a file or a dir, create dir
    mkdir (out)
end
if isdir(out)
    out=[out '/pipeline3.nii.gz'];
end
if exist(out,'file')
    delete(out)
end
save_nii(make_nii(MRIdenoised3,voxelsize,[],16),out)

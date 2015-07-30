function call_swi(in1,in2,out,posflag)
% Calling susceptibility weighted imaging filter
%
% - (C) 2015 Michael Eager (michael.eager@monash.edu)
% - Monash Biomedical Imaging

[a,b,c] = fileparts(mfilename('fullpath')) ;
[a,b,c] = fileparts(a) ;
root_path=a;
addpath(fullfile(root_path,'matlab'))
addpath(fullfile(root_path,'matlab/NIFTI'))
addpath(fullfile(root_path, 'matlab/Agilent/'))

display('Calling SWI')

if nargin == 3
    posflag=0;
end
    

voxelsize=[];
ksp1=[];ksp2=[];

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

if ~isempty(in2)
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
    %    return
end
end


if ~isempty(in2) && sum(voxelsize1) ~= sum(voxelsize2)
  display(['Voxelsizes don''t match: ' str2num(voxelsize1) ' ' ...
           str2num(voxelsize2)])
  return
end
voxelsize=voxelsize1;


[pha1, swi_n1, swi_p1, mag1] = phaserecon_v1(ksp1,ksp1,0.4,1,0.05);
% Necessary translations to match FDF images
swi_n1=flipdim(flipdim(flipdim(swi_n1,1),2),3);
swi_n1=circshift(swi_n1,[1,1,1]);

if ~isempty(ksp2)
    [pha2, swi_n2, swi_p2, mag2] = phaserecon_v1(ksp2,ksp2,0.4,1, ...
                                                 0.05);
    swi_n2=flipdim(flipdim(flipdim(swi_n2,1),2),3);
    swi_n2=circshift(swi_n2,[1,1,1]);
end

if exist(out,'file')~=2 && ~isdir(out)
    %if not a file or a dir, create dir
    mkdir (out)
end
if isdir(out)
    out=[out '/swi_neg.nii.gz'];
end
if exist(out,'file')
    delete(out)
end

if isempty(ksp2)
    save_nii(make_nii(swi_n1,voxelsize,[],16),out)
else
    save_nii(make_nii((swi_n1+swi_n2)/2,voxelsize,[],16),out)
end

if posflag == 1
   swi_p1=flipdim(flipdim(flipdim(swi_p1,1),2),3);
   swi_p1=circshift(swi_p1,[1,1,1]);
 
    if isempty(ksp2)
        save_nii(make_nii(swi_p1,voxelsize,[],16),out)
    else
        swi_p2=flipdim(flipdim(flipdim(swi_p2,1),2),3);
        swi_p2=circshift(swi_p2,[1,1,1]);
        save_nii(make_nii((swi_p1+swi_p2)/2,voxelsize,[],16),out)
    end
end

function call_mci(in1,in2,out)
% Calling MCI - max contrast imaging
%
% - (C) 2015 Michael Eager (michael.eager@monash.edu)
% - Monash Biomedical Imaging

[a,b,c] = fileparts(mfilename('fullpath')) ;
[a,b,c] = fileparts(a) ;
root_path=a;
addpath(fullfile(root_path,'matlab'))
addpath(fullfile(root_path,'matlab/NIFTI'))
addpath(fullfile(root_path, 'matlab/Agilent/'))

display('Calling MCI')
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
mag1=flipdim(flipdim(flipdim(mag1,1),2),3);
mag1=circshift(mag1,[1,1,1]);
pha1=flipdim(flipdim(flipdim(pha1,1),2),3);
pha1=circshift(pha1,[1,1,1]);
if ~isempty(ksp2)
    [pha2, swi_n2, swi_p2, mag2] = phaserecon_v1(ksp2,ksp2,0.4,1, ...
                                                 0.05);
    mag2=flipdim(flipdim(flipdim(mag2,1),2),3);
    mag2=circshift(mag2,[1,1,1]);
    pha2=flipdim(flipdim(flipdim(pha2,1),2),3);
    pha3=circshift(pha2,[1,1,1]);
end

stdmask=stdfilt(mag1);
posmask=find(pha1>0);
stdmask=stdfilt(mag1);
stdphmask=stdfilt(pha1);
gbmask =getbiggestobject(stdphmask<0.013 | stdmask>2000);
mask=(gbmask.*(pha1>0));

r_mag = mean(mag1(mask==1));
r_pha = mean(pha1(mask==1));

    mci_mag1 = mci(mag1,pha1,[r_mag, r_pha]);




if exist(out,'file')~=2 && ~isdir(out)
    %if not a file or a dir, create dir
    mkdir (out)
end
if isdir(out)
    out=[out '/mci.nii.gz'];
end
if exist(out,'file')
    delete(out)
end

if isempty(ksp2)
    save_nii(make_nii(mci_mag1,voxelsize,[],16),out)
else
    mci_mag2 = mci(mag2,pha2,[r_mag,rpha]);
    save_nii(make_nii((mci_mag1+mci_mag2)/2,voxelsize,[],16),out)
end




function output = mci(magnitude,phase,reference)
% MCI create Maximum Contrast Image
%
% OUTPUT = MCI(MAGNITUDE,PHASE,REFERENCE)
%
% MAGNITUDE = magnitude image
% PHASE = phase image
% REFERENCE = reference point
%
% Created by Zhaolin Chen 
% Adapted by Amanda Ng on 11 March 2009
% Updated by Michael Eager July 2015
    
    phasereg=phase.*0;
    phasereg(phase>0)=phase(phase>0)-reference(2);
    phasereg(phase<0)=phase(phase<0)+reference(2);
    
    output = sqrt((magnitude-reference(1)).^2 + (phasereg).^2);

function call_swi(in1,in2,out,order,preprocess,saveRI,swineg,swipos)
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
%% Clean input strings
in1 = regexprep(in1,'["\[\]]','');  
if ~isempty(in2)
 if isstr(in2)
        in2 = regexprep(in2,'["\[\]]','');
    else
        in2=[];
    end
end
out = regexprep(out,'["\[\]]',''); %"


display('Calling SWI')
display(in1)
display (in2)
display (out)

if nargin < 8
    swipos=0;
end
if nargin < 7
    swipneg=0;
end
if nargin < 6
    saveRI=0;
end
if nargin < 5
    preprocess=0;    
end
if nargin < 4
    order=0;    
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
    display(['Cannot find second image ' in2])
    %    return
end
end


if ~isempty(in2) && sum(voxelsize1) ~= sum(voxelsize2)
  display(['Voxelsizes don''t match: ' str2num(voxelsize1) ' ' ...
           str2num(voxelsize2)])
  return
end
voxelsize=voxelsize1;


[pha, swi_n, swi_p, mag] = phaserecon_v1(ksp1,ksp1,0.4,1,0.05);
% Necessary translations to match FDF images
swi_n=flipdim(flipdim(flipdim(swi_n,1),2),3);
swi_n=circshift(swi_n,[1,1,1]);
img=mag.*exp(1i*pha);

if ~isempty(ksp2)
    [pha2, swi_n2, swi_p2, mag2] = phaserecon_v1(ksp2,ksp2,0.4,1, ...
                                                 0.05);
    swi_n2=flipdim(flipdim(flipdim(swi_n2,1),2),3);
    swi_n2=circshift(swi_n2,[1,1,1]);
swi_n = (swi_n+swi_n2)/2;
    swi_p = (swi_p+swi_p2)/2;
img = (img+ (mag2.*exp(1i*pha2)))/2;

		
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


    save_nii(make_nii(swi_n,voxelsize,[],16),out)


if swipos
   swi_p=flipdim(flipdim(flipdim(swi_p,1),2),3);
   swi_p=circshift(swi_p,[1,1,1]);
   outp = regexprep(out,'neg','pos');

        save_nii(make_nii(swi_p,voxelsize,[],16),outp)

end

if saveRI
        save_nii(make_nii(real(img),voxelsize,[],16),regexprep(out,'neg','real'))
        save_nii(make_nii(imag(img),voxelsize,[],16),regexprep(out,'neg','imag'))

end
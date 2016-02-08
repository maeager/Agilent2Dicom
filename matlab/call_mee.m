function call_mee(in1,in2,out,porder,preprocess,saveRI,useswi)
% Calling MEE - multi-echo enhancement
%
% - (C) 2015 Michael Eager (michael.eager@monash.edu)
% - Monash Biomedical Imaging

[a,b,c] = fileparts(mfilename('fullpath')) ;
[a,b,c] = fileparts(a) ;
root_path=a;
addpath(fullfile(root_path,'./matlab'))
addpath(fullfile(root_path,'./matlab/NIFTI'))
addpath(fullfile(root_path, './matlab/Agilent/'))

display('Calling MCI')
if nargin < 4 
    porder=3;
end
if nargin < 5
    preprocess=0;
end
if nargin < 6
    saveRI=0;
end
if nargin < 7
    useswi=0;
end
%% Clean input strings
in1 = regexprep(in1,'["\[\]]','');
if ~isempty(in2)
 if isstr(in2)
        in2 = regexprep(in2,'["\[\]]','');
    else
        in2=[];
    end
end
out = regexprep(out,'["\[\]]','');


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

%% If using preprocessing
% add content

%% Homodyne filter

[pha1, swi_n1, swi_p1, mag1] = phaserecon_v1(ksp1,ksp1,0.4,1,0.05);
% Necessary translations to match FDF images



if ~useswi
    mag1=flipdim(flipdim(flipdim(mag1,1),2),3);
    mag1=circshift(mag1,[1,1,1]);
    pha1=flipdim(flipdim(flipdim(pha1,1),2),3);
    pha1=circshift(pha1,[1,1,1]);
    if ~isempty(ksp2)
        [pha2, swi_n2, swi_p2, mag2] = phaserecon_v1(ksp2,ksp2,0.4,1, 0.05);
        mag2=flipdim(flipdim(flipdim(mag2,1),2),3);
        mag2=circshift(mag2,[1,1,1]);
        pha2=flipdim(flipdim(flipdim(pha2,1),2),3);
        pha2=circshift(pha2,[1,1,1]);
    end

    mee_mag1 = mee((mag1.*exp(1i*pha1)),porder);
else
    swi1=flipdim(flipdim(flipdim(swi_n1,1),2),3);
    swi1=circshift(swi1,[1,1,1]);
    mee_mag1 = mee(swi1,porder);
end


if exist(out,'file')~=2 && ~isdir(out)
    %if not a file or a dir, create dir
    mkdir (out)
end
if isdir(out)
    out=[out '/mee_magn.nii.gz'];
end
if exist(out,'file')
    delete(out)
end

if isempty(ksp2)
    save_nii(make_nii(abs(mee_mag1),voxelsize,[],16),out)
else
    if ~useswi
        mee_mag2 = mee((mag2.*exp(1i*pha2)),porder);
    else
        swi2=flipdim(flipdim(flipdim(swi_n2,1),2),3);
        swi2=circshift(swi2,[1,1,1]);
        mee_mag2 = mee(swi2,porder);
    end
    save_nii(make_nii(abs((mee_mag1+mee_mag2)/2.0),voxelsize,[],16),out)
end


if saveRI && ~useswi
    if isempty(ksp2)
        save_nii(make_nii(real(mee_mag1),voxelsize,[],16),regexprep(out,'magn','real'))
        save_nii(make_nii(imag(mee_mag1),voxelsize,[],16),regexprep(out,'magn','imag'))
    else
        save_nii(make_nii(real((mee_mag1+mee_mag2)/2.0),voxelsize, ...
                          [],16),regexprep(out,'magn','real'))
        save_nii(make_nii(imag((mee_mag1+mee_mag2)/2.0),voxelsize, ...
                          [],16),regexprep(out,'magn','imag'))
    end
end




function output = mee(img,order)
% Multi-echo enhancement
%
% img = 5D magnitude image
% 
% Created by Michael Eager July 2015

    sz=size(img);output=[];
    if len(sz)==4
        sz(5)=sz(4);sz(4)=1;
        img = reshape(img,sz);
    elseif len(sz)~=5
        return
    end
    if sz(4)!=1 && sz(5) == 1
        output=img;
        return
    end

    img = (abs(img));

    if nargin==1
        p=3.0;
    else
        p=order;
    end

    mee = (sum(img.^-p,5)./sz).^(-1/p);

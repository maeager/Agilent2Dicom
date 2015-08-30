function call_pipeline1(in,out,in2,NLfilter,hfinal)
% Calling non-local means pipeline 1
%  Use rician noise estimator to calculate std of noise in one MRI
% magnitude image
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
run (fullfile(root_path,'../matlab/NLmeans/vlfeat/toolbox/vl_setup.m'))

display('Calling non-local means filter pipeline 1')

voxelsize=[];
hfinal=[];
NLfilter=[];

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

if nargin == 3 && ~isempty(in2)
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

img=squeeze(img);

if length(size(img)) == 4
    display 'Reducing 4D image down to 3D'
    img=img(:,:,:,1);
elseif length(size(img)) == 5
    display 'Reducing 5D image down to 3D'
    img=img(:,:,:,1,1);
elseif length(size(img)) > 5
    display 'Unable to process images greater than 5D'
    return
end

display 'Calling pipeline 1'

tic(),MRIdenoised1=pipeline1(NormaliseImage2(img)*256,NLfilter, ...
                             hfinal);toc()

if exist(out,'file')==2
    delete(out)
    denoised_file = out
    [a,b,c] = fileparts(out) ;
    raw_file = [ a, '/raw_average.nii.gz'];	
else
   if ~isdir(out) 
      mkdir(out)
   end
   denoised_file = [out '/pipeline1.nii.gz'];
   raw_file = [out, '/raw_average.nii.gz'];
end
% Save raw input image average if not already saved
if exist(raw_file,'file') ~= 2
    save_nii(make_nii(img,voxelsize,[],16),raw_file)
end
%Save the denoised image
save_nii(make_nii(MRIdenoised1,voxelsize,[],16),denoised_file)

function make_diffusion(diffusion_path)



addpath ~/Code/Agilent
addpath ~/Monash016/eagerm/mrtrix-0.2.12/matlab
% cd ~/Monash016/eagerm/Agilent2Dicom

if nargin < 1
   diffusion_path='../example_data/s_2014061202/epi-dir30_01.img'
end


%%
[im,hdr] = readfdf(diffusion_path);
[pp, acq] = readpp(fullfile(diffusion_path,'procpar'), true);
image.data=single(squeeze(im));



%% 
% DW=zeros(31,4);
% jj=2
% for ii = 5:2:64
% DW(jj,1:3) = abs(eig([pp.bvalpp(ii) pp.bvalrp(ii) pp.bvalsp(ii); pp.bvalrs(ii) pp.bvalrr(ii) pp.bvalrs(ii); pp.bvalsp(ii) pp.bvalrp(ii) pp.bvalss(ii)])');
% DW(jj,4)=pp.bvalSave; %bval(ii);
% jj=jj+1;
% end

DW=zeros(31,4);
jj=2;
for ii = 5:2:64
DW(jj,1:3) = [pp.dpe(ii)  pp.dro(ii) pp.dsl(ii)];
DW(jj,4)=pp.bvalSave; %bval(ii);
jj=jj+1;
end


%% Diffusion weighted image
 image.DW_scheme=DW; % a NDWx4 matrix of gradient directions [optional]
 write_mrtrix(image,'../output_mif/diff.mif')
 image.vox= hdr.voxelsize;  %        N-vector of voxel sizes (in mm) (default: { 2 }) [optional]
    % image.comments= 'Diffusion weighted image created by Agilent2Dicom in MATLAB.' %   a cell array of strings [optional]
    % image.datatype='float32';  %   the datatype specifier (default: float32) [optional]
 image.transform = [hdr.orientation(1:3) 0; hdr.orientation(4:6) 0; hdr.orientation(7:9) 0; 0 0 0 1];   %  a 4x4 matrix [optional]
 write_mrtrix(image,'../output_mif/diff2.mif')


%% Write Mask

aveDif= mean(image.data(:,:,:,2:end),4);
composition = imdilate(single(aveDif>250),strel('ball',3,3));
aveDwi.data = single(composition > 3.05);
write_mrtrix(aveDwi,'../output_mif/mask_DWI1.mif')
aveDWi.vox= hdr.voxelsize; 
aveDwi.transform = [hdr.orientation(1:3) 0; hdr.orientation(4:6) 0; hdr.orientation(7:9) 0; 0 0 0 1]; 
write_mrtrix(aveDwi,'../output_mif/mask_DWI2.mif')

composition2 = imerode(single(aveDif > 250),strel('ball',3,3));
comp2temp = imfill(imdilate(single(composition2 > -2.5),strel('ball',3,3)),'holes');
aveDwi2.data=single(comp2temp>3.05)
write_mrtrix(aveDwi2,'../output_mif/mask_DWI.mif')
aveDWi2.vox= hdr.voxelsize; 
aveDwi2.transform = [hdr.orientation(1:3) 0; hdr.orientation(4:6) 0; hdr.orientation(7:9) 0; 0 0 0 1]; 
write_mrtrix(aveDwi2,'../output_mif/mask2_DWI.mif')


%% ADC

aveDif= mean(image.data(:,:,:,2:end),4);
adc.data = single(aveDif./(image.data(:,:,:,1)));
write_mrtrix(adc,'../output_mif/adc_DWI.mif')
adc.vox= hdr.voxelsize; 
adc.transform = [hdr.orientation(1:3) 0; hdr.orientation(4:6) 0; hdr.orientation(7:9) 0; 0 0 0 1]; 
write_mrtrix(adc,'../output_mif/adc_DWI2.mif')







function make_diffusion(diffusion_path)
%% make_diffusion quick calculation of diffusion images
%
%
% Copyright (C) 2014 Michael Eager  (michael.eager@monash.edu)

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.



addpath ~/Code/Agilent
addpath /usr/local/mrtrix/0.3.12/matlab
% cd ~/Monash016/eagerm/Agilent2Dicom

if nargin < 1
   diffusion_path='../example_data/s_2014061202/epi-dir30_01.img';
end


%%
[im,hdr] = readfdf(diffusion_path);
[pp, ~] = readpp(fullfile(diffusion_path,'procpar'), true);
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







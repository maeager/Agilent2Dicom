#!/usr/bin/env bash 
#
# -(c) Michael Eager 2014 (michael.eager@monash.edu)

set -u # nounset
set -e # errexit
set -o # pipefail
# exec 2>> $(dirname $0)/test_diffusion_error.log
# set -x

outpath='../output_mif'

module load  mrtrix/0.2.12 matlab/r2013a

[ ! -d ${outpath} ] &&  mkdir ${outpath}

rm -f ${outpath}/*.mif
rm -f ${outpath}/*.tck

[ ! -d ../example_data/s_2014061202/epi-dir30_01.img ] && (echo "No heart tissue FDF directory."; exit 1)

## Create diffusion weighted image from FDF

# Standard method
# ./fdf2dcm.sh -i ../example_data/s_2014061202/epi-dir30_01.img -o ../output_data/hearttissue.dcm
# mrconvert ../output_data/hearttissue.dcm/ ../output_mif/diff.mif
# Masking
# average ${outpath}/diff.mif -axis 3 - | threshold - -| median3D - - | median3D - ${outpath}/mask_DWI.mif


# Quick and dirty MATLAB method with masking
matlab -nosplash -nodesktop -r "addpath "$(dirname $0)";make_diffusion;quit"

mrinfo ${outpath}/diff.mif | grep 'DW scheme'
[ $? -ne 0 ] && (echo 'Diffusion MIF file does not contain DW scheme'; exit 1)


## Diffusion tensor imaging
dwi2tensor ${outpath}/diff.mif ${outpath}/dt_DWI.mif

## FA map
 tensor2FA ${outpath}/dt_DWI.mif  - | mrmult - ${outpath}/mask_DWI.mif ${outpath}/fa_DWI.mif

## eigenvector map
 tensor2vector ${outpath}/dt_DWI.mif  - | mrmult - ${outpath}/fa_DWI.mif ${outpath}/ev_DWI.mif

## Constrained spherical deconvolution (CSD), masking single voxels
 erode ${outpath}/mask_DWI.mif -npass 3 - | mrmult ${outpath}/fa_DWI.mif - - | threshold - -abs 0.7 ${outpath}/sf_DWI.mif


## Fibre tracking
streamtrack DT_STREAM ${outpath}/diff.mif -seed ${outpath}/mask_DWI.mif -mask ${outpath}/mask_DWI.mif -num 100000 ${outpath}/DWI_wholetissue.tck

tracks2prob ${outpath}/DWI_wholetissue.tck -colour -vox 0.5 ${outpath}/tdi.mif


## Response function coefficient
 estimate_response ${outpath}/diff.mif ${outpath}/sf_DWI.mif ${outpath}/response.txt

## display response
 disp_profile -response ${outpath}/response.txt

##
estimate_response ${outpath}/diff.mif ${outpath}/sf_DWI.mif -lmax 6 ${outpath}/response-6max.txt

## CSD computation:
csdeconv ${outpath}/diff.mif ${outpath}/response-6max.txt -lmax 8 -mask ${outpath}/mask_DWI.mif ${outpath}/CSD8.mif


streamtrack SD_PROB ${outpath}/CSD8.mif -seed ${outpath}/mask_DWI.mif -mask ${outpath}/mask_DWI.mif ${outpath}/sd_prob.tck -num 50000
 tracks2prob ${outpath}/sd_prob.tck -colour -vox 0.5 ${outpath}/tdi-sd.mif
## Track density
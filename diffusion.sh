#!/usr/bin/env bash 
## diffusion.sh - Correction of DCMULTI conversion of dicoms
#
# Copyright (C) 2014 Michael Eager  (michael.eager@monash.edu)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

set -u # nounset
set -e # errexit
set -o # pipefail
# exec 2>> $(dirname $0)/test_diffusion_error.log
# set -x

inputdicom=$1
if [ $2 -eq 1 ];then
    outpath
else
    outpath=${inputdicom}/../mrtrix
fi
MAKENIFTI=0
if [ $# -eq 3 ];then
    MAKENIFTI=1
    niioutpath=${outpath}/../nifti
fi
STOPATMASK=0
if [ $# -eq 4 ];then
	STOPATMASK=1
fi


if test ${MASSIVEUSERNAME+defined}; then
    module load  mrtrix matlab
fi

[ ! -d ${outpath} ] &&  mkdir ${outpath}

rm -f ${outpath}/*.mif
rm -f ${outpath}/*.tck

[ ! -d "${inputdicom}" ] && (echo "No heart tissue FDF directory."; exit 1)

## Create diffusion weighted image from FDF

# Standard method
# ./fdf2dcm.sh -i ../example_data/s_2014061202/epi-dir30_01.img -o ../output_data/hearttissue.dcm
mrconvert ${inputdicom} ${outpath}/diff.mif -datatype float32

# Masking
average ${outpath}/diff.mif -axis 3 - | threshold - -| median3D - - | median3D - ${outpath}/mask_DWI.mif


# Quick and dirty MATLAB method with masking
# matlab -nosplash -nodesktop -r "addpath "$(dirname $0)";make_diffusion;quit"

mrinfo ${outpath}/diff.mif | grep 'DW scheme'
[ $? -ne 0 ] && (echo 'Diffusion MIF file does not contain DW scheme'; exit 1)

if [ $STOPATMASK -eq 1 ]; then 
    exit 0;
fi

## Diffusion tensor image
dwi2tensor ${outpath}/diff.mif ${outpath}/dt_DWI.mif

## FA map
tensor2ADC ${outpath}/dt_DWI.mif  - | mrmult - ${outpath}/mask_DWI.mif ${outpath}/adc_DWI.mif

## FA map
tensor2FA ${outpath}/dt_DWI.mif  - | mrmult - ${outpath}/mask_DWI.mif ${outpath}/fa_DWI.mif

## eigenvector map
tensor2vector ${outpath}/dt_DWI.mif  - | mrmult - ${outpath}/fa_DWI.mif ${outpath}/ev_DWI.mif

## Constrained spherical deconvolution (CSD), masking single voxels
erode ${outpath}/mask_DWI.mif -npass 3 - | mrmult ${outpath}/fa_DWI.mif - - | threshold - -abs 0.7 ${outpath}/sf_DWI.mif


## Fibre tracking
streamtrack DT_STREAM ${outpath}/diff.mif -seed ${outpath}/mask_DWI.mif -mask ${outpath}/mask_DWI.mif -num 100000 ${outpath}/DWI_wholetissue.tck
## Track density
tracks2prob ${outpath}/DWI_wholetissue.tck -colour -vox 0.5 ${outpath}/tdi.mif


## Response function coefficient
estimate_response ${outpath}/diff.mif ${outpath}/sf_DWI.mif ${outpath}/response.txt

## display response
disp_profile -response ${outpath}/response.txt

##
estimate_response ${outpath}/diff.mif ${outpath}/sf_DWI.mif -lmax 6 ${outpath}/response-6max.txt

## CSD computation:
csdeconv ${outpath}/diff.mif ${outpath}/response-6max.txt -lmax 8 -mask ${outpath}/mask_DWI.mif ${outpath}/CSD8.mif
## Create tracks from CSD
streamtrack SD_PROB ${outpath}/CSD8.mif -seed ${outpath}/mask_DWI.mif -mask ${outpath}/mask_DWI.mif ${outpath}/sd_prob.tck -num 5000
## Track density with CSD
tracks2prob ${outpath}/sd_prob.tck -colour -vox 0.5 ${outpath}/tdi-sd.mif




if [ $MAKENIFTI -eq 1 ]; then

    [ ! -d $niioutpath ] && mkdir -p $niioutpath
    ## fahdshariff.blogspot.com.au/2012/10/shell-scripting-best-practices.html
    ## don't use for i in $(ls ...
    find "${outpath}" -type f -name "*.mif" -maxdepth 1 -print0 | while IFS= read -r -d $'\0' miffile; do
	[ ! -e "$miffile" ] || continue
	mrconvert "$miffile" "$niioutpath"/$(basename "$miffile" .mif).nii
    done

fi


#  mkdir dti_rep1_01.mif
#  cd dti_rep1_01.mif
#  mrconvert ../dti_rep1_01.dcm dwi_01.mif -datatype float32
#  dwi2tensor dwi_01.mif dt_DWI_01.mif
#  average dwi_01.mif -axis 3 - | threshold - -| median3D - -| median3D - mask_DWI_01.mif
#  tensor2FA dt_DWI_01.mif - |mrmult - mask_DWI_01.mif fa_DWI_01.mif
#  tensor2vector dt_DWI_01.mif - |mrmult - fa_DWI_01.mif ev_DWI_01.mif
#  erode mask_DWI_01.mif -npass 3 - | mrmult fa_DWI_01.mif - -| threshold - -abs 0.7 sf_SWI_01.mif
#  estimate_response dwi_01.mif sf_DWI_01.mif  -lmax 6 response.txt
#  csdeconv dwi_01.mif  response.txt -lmax 8 -mask mask_DWI_01.mif CSD8.mif
#  tensor2ADC dt_DWI_01.mif - | mrmult - mask_DWI_01.mif adc_DWI_01.mif

# # Stream tracking

#  cd ..
#  mkdir dti_rep1_01.nii
#  for miffile in $(ls dti_rep1_01.mif/*.mif) ;do
#      mrconvert -scale 1000 -info $miffile dti_rep1_01.nii/$(basename $miffile .mif).nii  -datatype float32
#  done

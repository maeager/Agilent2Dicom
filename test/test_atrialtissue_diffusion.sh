#!/usr/bin/env bash 
#
# - Michael Eager  (michael.eager@monash.edu)

  # Copyright (C) 2014 Michael Eager  (michael.eager@monash.edu)

  # This program is free software: you can redistribute it and/or modify
  # it under the terms of the GNU General Public License as published by
  # the Free Software Foundation, either version 3 of the License, or
  # (at your option) any later version.

  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.

  # You should have received a copy of the GNU General Public License
  # along with this program.  If not, see <http://www.gnu.org/licenses/>.

set -x
set -v
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

Diffusion_FDF=${1:-../example_data/s_2014061202/epi-dir30_01.img} # use heart tissue as default
[ ! -d $Diffusion_FDF ] && (echo "No diffusion FDF directory."; exit 1)
Diffusion_DCM=${2:-../output_data/hearttissue.dcm} # use heart tissue as default
[ ! -d $Diffusion_DCM ] && (echo "No diffusion DCM directory."; exit 1)

## Create diffusion weighted image from FDF

# Standard method
if [ -f ${output_path}/diff.mif ]; then
    ./fdf2dcm.sh -i ${Diffusion_FDF} -i ${Diffusion_DCM}
    mrconvert -info -datatype float32 $Diffusion_DCM ${output_path}/diff.mif
fi



# Masking
if [ -f ../output_mif/mask_DWI.mif ]; then
    average ${outpath}/diff.mif -axis 3 - | threshold - -| median3D - - | median3D - ${outpath}/mask_DWI.mif
fi

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
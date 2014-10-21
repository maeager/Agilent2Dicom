#!/usr/bin/env bash 


EXAMPLEDATA=../..
fdffolders=$(find $EXAMPLEDATA -type d -name "*.img")

fdfarray=($fdffolders)

for fdf in "${fdfarray[@]}"
do
    output_path=$(echo ${filen} | sed -e 's/img/dcm/' -e 's/example/output/')
    [ $fdf == $output_path ] && continue
    ../fdf2dcm.sh -v -i "$fdf" -o "${output_path}"
    mrinfo "${output_path}"

done


module load mrtrix

if [ -d ../output_nii ]; 
then
    rm -f ..output_nii/*.nii
else
    mkdir ../output_nii
fi
## standard 2d

./Agilent2Dicom/fdf2dcm.sh -v -i ~/Monash016/amanda/ExampleAgilentData/standard2d/ -o ./output_data/standard2d
mrconvert -info  output_data/standard2d/ output_nii/standard2d.nii
fslview output_nii/standard2d.nii


## Multi-echo 3D

./Agilent2Dicom/fdf2dcm.sh -v -i  ~/Monash016/amanda/ExampleAgilentData/multiecho3d_magonly/ -o ./output_data/multiecho3d_magonly
mrconvert -info  ./output_data/multiecho3d_magonly/ output_nii/ME3d.nii
fslview output_nii/ME3d.nii


## Multi-echo 2D Mag and phase

./Agilent2Dicom/fdf2dcm.sh -v -i ~/Monash016/amanda/ExampleAgilentData/multiecho2d_magandphase/ -o  ./output_data/multiecho2d_magandphase
mrconvert -info  output_data/multiecho2d_magandphase/magnitude.dcm output_nii/ME2d_mag.nii
mrconvert -info  output_data/multiecho2d_magandphase/phase.dcm output_nii/ME2d_phase.nii
fslview output_nii/ME2d_mag.nii output_nii/ME2d_phase.nii


## CINE

./Agilent2Dicom/fdf2dcm.sh -v -i  ~/Monash016/amanda/ExampleAgilentData/cine -o ./output_data/cine
mrconvert -info  output_data/cine/ output_nii/CINE.nii
fslview output_nii/CINE.nii


## ASL

./Agilent2Dicom/fdf2dcm.sh -v -i ./example_data/1008.2.40.4.1.1/ASL_se_06.img -o ./output_data/ASL_se_06.dcm
mrconvert -info  output_data/ASL_se_06.dcm/ output_nii/ASL.nii
fslview output_nii/ASL.nii

#Diffusion

./Agilent2Dicom/fdf2dcm.sh -v -i ~/Monash016/amanda/ExampleAgilentData/diffusion/ -o ./output_data/diffusion
mrconvert -info  output_data/diffusion/ output_nii/Diffusion.nii
fslview output_nii/Diffusion.nii







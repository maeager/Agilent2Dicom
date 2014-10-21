#!/usr/bin/env bash
#
# Testing dicom outputs with mricron Nifti converter
#
# - Monash Biomedical Imaging
# - Michael Eager (michael.eager@monash.edu)

#   Copyright (C) 2014 Michael Eager  (michael.eager@monash.edu)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Example ~/.dcm2nii/dcm2nii.ini:
# [ BOOL]
# DebugMode=0
# UntestedFeatures=0
# UINT16toFLOAT32=1
# Verbose=0
# Anonymize=1
# AnonymizeSourceDICOM=0
# AppendAcqSeries=1
# AppendDate=1
# AppendFilename=0
# AppendPatientName=0
# AppendProtocolName=1
# AutoCrop=0
# CollapseFolders=1
# createoutputfolder=0
# CustomRename=0
# enablereorient=1
# OrthoFlipXDim=0
# EveryFile=1
# fourD=1
# Gzip=1
# ManualNIfTIConv=1
# PhilipsPrecise=0
# RecursiveUseNameAppend=0
# SingleNIIFile=1
# SPM2=0
# Stack3DImagesWithSameAcqNum=0
# Swizzle4D=1
# UseGE_0021_104F=0
#
# [ INT]
# MaxReorientMatrix=1023
# MinReorientMatrix=200
# RecursiveFolderDepth=5
# OutDirMode=0
# SiemensDTIUse0019If00181020atleast=15
# SiemensDTINoAngulationCorrectionIf00181020atleast=1000
# SiemensDTIStackIf00181020atleast=15
#
# [ STR]
# OutDir=/home/eagerm/

if test ${MASSIVEUSERNAME+defined}; then 
    module load mricron fsl
fi

for i in ../output_dcm2nii/*.nii;
do
    rm -f $i
done

        echo "MRICRON convert to Nifti"
        dcm2nii -n -o ../output_dcm2nii/ ../output_data/standard2d/
        dcm2nii -n -o ../output_dcm2nii/  ../output_data/multiecho3d_magonly/ 
        dcm2nii -n -o ../output_dcm2nii/  ../output_data/multiecho2d_magandphase/magnitude.dcm/ 
        dcm2nii -n -o ../output_dcm2nii/  ../output_data/multiecho2d_magandphase/phase.dcm/ 
        dcm2nii -n -o ../output_dcm2nii/  ../output_data/cine/ 
        dcm2nii -n -o ../output_dcm2nii/  ../output_data/ASL_se_06.dcm/ 
        dcm2nii -n -o ../output_dcm2nii/  ../output_data/diffusion/ 
        dcm2nii -n -o ../output_dcm2nii/  ../output_data/dti/ 
        dcm2nii -n -o ../output_dcm2nii/  ../output_data/epip/ 
        dcm2nii -n -o ../output_dcm2nii/  ../output_data/j6/ 


 # call: a2d_check <inputdir> <fullniftifile> 
function a2d_check(){
echo "Agilent2Dicom Check DCM and NII"; 
[ -f $2 ] && fslhd $2 | grep -e '^dim[ 1234]' | \
awk 'BEGIN{i=0} {a[ i]=$$2; ++i} END{printf("NIFTI  Dimensions:        "); for (i=0;i<4;++i){printf("%d ",a[ i]);if (i!=3){printf("x ")}}; printf("\n")}' ;
printf "DICOM";
mrinfo $1 2> /dev/null | grep 'Dimensions'
}



a2d_check  ../output_data/standard2d/                                   ../output_dcm2nii/standard2d.nii 
a2d_check  ../output_data/multiecho3d_magonly/                   	../output_dcm2nii/multiecho3d_magonly.nii 
a2d_check  ../output_data/multiecho2d_magandphase/magnitude.dcm/	../output_dcm2nii/multiecho2d_magandphase/magnitude.nii
a2d_check  ../output_data/multiecho2d_magandphase/phase.dcm/     	../output_dcm2nii/multiecho2d_magandphase/phase.nii 
a2d_check  ../output_data/cine/						../output_dcm2nii/cine.nii
a2d_check  ../output_data/ASL_se_06.dcm/                         	../output_dcm2nii/ASL_se_06.nii
a2d_check  ../output_data/diffusion/					../output_dcm2nii/diffusion.nii
a2d_check  ../output_data/dti/                                   	../output_dcm2nii/dti.nii
a2d_check  ../output_data/epip/						../output_dcm2nii/epip.nii
a2d_check  ../output_data/j6/                                    	../output_dcm2nii/j6.nii
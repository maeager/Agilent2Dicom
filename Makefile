
REQUIREMENTS= ./dpush mrconvert ./fdf2dcm.sh ./agilent2dicom storescu dcmodify

check:
	for req in $(REQUIREMENTS); do \
	   which $$req > /dev/null || echo "Missing dependency $$req"; \
	done

run_kidney:
	./fdf2dcm.sh -v -i ~/Monash016/RatKidney/Agilent/20120522/kidney512iso_01.img -o ../output_data/kidney512iso.dcm

test_kidney:
	dciodvfy   ../output_data/kidney512iso.dcm/0001.dcm 2>&1 >/dev/null | grep Err 

setup:
	mkdir ../output_data
	mkdir ../output_nii

cleanup:
	rm -f ../output_nii/*.nii

run_standard2d:
	./fdf2dcm.sh -v -i ~/Monash016/amanda/ExampleAgilentData/standard2d/ -o ../output_data/standard2d
# [ -f ../output_nii/standard2d.nii ] 
	rm -f ../output_nii/standard2d.nii
	echo "MRTIX convert to Nifti"
	mrconvert -info  ../output_data/standard2d/ ../output_nii/standard2d.nii
	fslview ../output_nii/standard2d.nii

test_standard2d:
	-dciodvfy   ../output_data/standard2d/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy   ../output_data/standard2d/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'

## Multi-echo 3D
run_me3d:
	./fdf2dcm.sh -v -i  ~/Monash016/amanda/ExampleAgilentData/multiecho3d_magonly/ -o ../output_data/multiecho3d_magonly
	-rm -f ../output_nii/ME3d.nii
	mrconvert -info  ../output_data/multiecho3d_magonly/ ../output_nii/ME3d.nii
	fslview ../output_nii/ME3d.nii

test_me3d:
	-dciodvfy   ../output_data/multiecho3d_magonly/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy  ../output_data/multiecho3d_magonly/0001.dcm 2>&1 >/dev/null | grep -e '^Warn'


## Multi-echo 2D Mag and phase
run_me2d:
	./fdf2dcm.sh -v -i ~/Monash016/amanda/ExampleAgilentData/multiecho2d_magandphase/ -o  ../output_data/multiecho2d_magandphase
	-rm -f ../output_nii/ME2d_mag.nii ../output_nii/ME2d_phase.nii
	mrconvert -info  ../output_data/multiecho2d_magandphase/magnitude.dcm/ ../output_nii/ME2d_mag.nii
	mrconvert -info  ../output_data/multiecho2d_magandphase/phase.dcm/ ../output_nii/ME2d_phase.nii
	fslview ../output_nii/ME2d_mag.nii ../output_nii/ME2d_phase.nii


test_me2d:
	-dciodvfy   ../output_data/multiecho2d_magandphase/magnitude.dcm/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy   ../output_data/multiecho2d_magandphase/magnitude.dcm/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'
	-dciodvfy   ../output_data/multiecho2d_magandphase/phase.dcm/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy   ../output_data/multiecho2d_magandphase/phase.dcm/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'



## CINE
run_cine:
	./fdf2dcm.sh -v -d -i ~/Monash016/amanda/ExampleAgilentData/cine -o ../output_data/cine
	-rm -f ../output_nii/CINE.nii
	mrconvert -info  ../output_data/cine/ ../output_nii/CINE.nii
	fslview ../output_nii/CINE.nii

test_cine:
	-dciodvfy   ../output_data/cine/0001.dcm 2>&1 >/dev/null | grep -e '^Err'
	-dciodvfy   ../output_data/cine/0001.dcm 2>&1 >/dev/null | grep -e '^Warn'


## ASL
run_asl:
	./fdf2dcm.sh -v -i ../example_data/1008.2.40.4.1.1/ASL_se_06.img -o ../output_data/ASL_se_06.dcm
	-rm -f ../output_nii/ASL.nii
	mrconvert -info  ../output_data/ASL_se_06.dcm/ ../output_nii/ASL.nii
	fslview ../output_nii/ASL.nii

test_asl:
	-dciodvfy   ../output_data/ASL_se_06.dcm/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy   ../output_data/ASL_se_06.dcm/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'


## Diffusion
run_diffusion:
	./fdf2dcm.sh -v -i ~/Monash016/amanda/ExampleAgilentData/diffusion/ -o ../output_data/diffusion
	-rm -f ../output_nii/Diffusion.nii
	mrconvert -info  ../output_data/diffusion/ ../output_nii/Diffusion.nii
	fslview ../output_nii/Diffusion.nii

test_diffusion:
	-dciodvfy  ../output_data/diffusion/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy  ../output_data/diffusion/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'



test_all: test_standard2d test_me3d test_me2d test_cine test_asl test_diffusion
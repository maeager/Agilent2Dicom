
REQUIREMENTS= ./dpush mrconvert ./fdf2dcm.sh ./agilent2dicom storescu dcmodify
MRVIEW="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/Monash016/eagerm/mrtrix-0.2.12/lib vglrun ~/Monash016/eagerm/mrtrix-0.2.12/bin/mrview"

FSLVIEW=fslhd

check:
	for req in $(REQUIREMENTS); do \
	   which $$req > /dev/null || echo "Missing dependency $$req"; \
	done

run_kidney:
	./fdf2dcm.sh -v -i ~/Monash016/RatKidney/Agilent/20120522/kidney512iso_01.img -o ../output_data/kidney512iso.dcm


test_kidney:
	-dciodvfy   ../output_data/kidney512iso.dcm/0001.dcm 2>&1 >/dev/null 

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
	$(FSLVIEW) ../output_nii/standard2d.nii

test_standard2d:
	-dciodvfy   ../output_data/standard2d/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy   ../output_data/standard2d/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'

## Multi-echo 3D
run_me3d:
	./fdf2dcm.sh -v -i  ~/Monash016/amanda/ExampleAgilentData/multiecho3d_magonly/ -o ../output_data/multiecho3d_magonly
	-rm -f ../output_nii/ME3d.nii
	mrconvert -info  ../output_data/multiecho3d_magonly/ ../output_nii/ME3d.nii -datatype float32
	$(FSLVIEW) ../output_nii/ME3d.nii

test_me3d:
	-dciodvfy   ../output_data/multiecho3d_magonly/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy  ../output_data/multiecho3d_magonly/0001.dcm 2>&1 >/dev/null | grep -e '^Warn'


## Multi-echo 2D Mag and phase
run_me2d:
	./fdf2dcm.sh -v -i ~/Monash016/amanda/ExampleAgilentData/multiecho2d_magandphase/ -o  ../output_data/multiecho2d_magandphase
	-rm -f ../output_nii/ME2d_mag.nii ../output_nii/ME2d_phase.nii
	mrconvert -info  ../output_data/multiecho2d_magandphase/magnitude.dcm/ ../output_nii/ME2d_mag.nii -datatype float32
	mrconvert -info  ../output_data/multiecho2d_magandphase/phase.dcm/ ../output_nii/ME2d_phase.nii -datatype float32
	$(FSLVIEW) ../output_nii/ME2d_mag.nii ../output_nii/ME2d_phase.nii


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
	$(FSLVIEW) ../output_nii/CINE.nii

test_cine:
	-dciodvfy   ../output_data/cine/0001.dcm 2>&1 >/dev/null | grep -e '^Err'
	-dciodvfy   ../output_data/cine/0001.dcm 2>&1 >/dev/null | grep -e '^Warn'


## ASL
run_asl:
	./fdf2dcm.sh -v -i ../example_data/1008.2.40.4.1.1/ASL_se_06.img -o ../output_data/ASL_se_06.dcm
	-rm -f ../output_nii/ASL.nii
	mrconvert -info  ../output_data/ASL_se_06.dcm/ ../output_nii/ASL.nii
	$(FSLVIEW) ../output_nii/ASL.nii

test_asl:
	-dciodvfy   ../output_data/ASL_se_06.dcm/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy   ../output_data/ASL_se_06.dcm/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'


## Diffusion
run_diffusion:
	./fdf2dcm.sh -v -i ~/Monash016/amanda/ExampleAgilentData/diffusion/ -o ../output_data/diffusion
	-rm -f ../output_nii/Diffusion.nii
	mrconvert -info  ../output_data/diffusion/ ../output_nii/Diffusion.nii
	$(FSLVIEW) ../output_nii/Diffusion.nii

test_diffusion:
	-dciodvfy  ../output_data/diffusion/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy  ../output_data/diffusion/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'



test_all: test_standard2d test_me3d test_me2d test_cine test_asl test_diffusion


## Diffusion DTI
run_dti:
	./fdf2dcm.sh -v -i ../example_data/s_2014051208/DTI_EPIP_J30_01.img/ -o ../output_data/dti/
	-rm -f ../output_nii/DTI.nii
	mrconvert -info  ../output_data/dti/ ../output_nii/DTI.nii
	$(FSLVIEW) ../output_nii/DTI.nii

test_dti:
	-dciodvfy  ../output_data/dti/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy  ../output_data/dti/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'


## Diffusion DTI
run_epip:
	./fdf2dcm.sh -v -i ../example_data/s_2014051208/epip_axi_300TR_01.img/ -o ../output_data/epip/
	-rm -f ../output_nii/EPI.nii
	mrconvert -info  ../output_data/epip/ ../output_nii/EPI.nii
	$(FSLVIEW) ../output_nii/EPI.nii

test_epip:
	-dciodvfy  ../output_data/epip/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy  ../output_data/epip/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'


## Diffusion DTI
run_J6:
	./fdf2dcm.sh -v -i ../example_data/s_2014051901/J6-500_01.img/ -o ../output_data/j6/
	-rm -f ../output_nii/J6.nii
	mrconvert -info  ../output_data/j6/ ../output_nii/J6.nii
	$(FSLVIEW) ../output_nii/J6.nii

test_J6:
	-dciodvfy  ../output_data/j6/0001.dcm 2>&1 >/dev/null | grep -e '^Error'
	-dciodvfy  ../output_data/j6/0001.dcm 2>&1 >/dev/null | grep -e '^Warning'


.PHONY: all run_kidney run_standard2d run_me3d run_me2d run_cine run_asl run_diffusion run_dti run_epip run_J6
all: run_kidney run_standard2d run_me3d run_me2d run_cine run_asl run_diffusion run_dti run_epip run_J6

.PHONY: test_all test_kidney test_standard2d test_me3d test_me2d test_cine test_asl test_diffusion test_dti test_epip test_J6
test_all: test_kidney test_standard2d test_me3d test_me2d test_cine test_asl test_diffusion test_dti test_epip test_J6

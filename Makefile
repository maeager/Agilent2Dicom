
REQUIREMENTS= python mrconvert storescu dcmodify dcmulti dciodvfy fslhd fslview
GL=vglrun # set to appropriate virtualgl 
MRVIEW=$(GL) mrview
MRCONVERT=mrconvert -info -datatype float32 
FSLHD=fslhd
FSLVIEW=$(GL) fslview 
DCM3TOOLS=$(shell /bin/ls -d ../dicom3tools_*/bin/*)
RM=/bin/rm -f
FDF2DCM=time ./fdf2dcm.sh -v 
PATH+=:$(DCM3TOOLS)

displenv:
	echo $(DCM3TOOLS)

.PHONY: check
check: 
	for req in $(REQUIREMENTS); do \
	   which $$req > /dev/null || echo "Missing dependency $$req"; \
	done
#	test -x $$req || (echo "$$req not found"; exit 1)
	python check.py


.PHONY: setup
setup:
	-mkdir ../output_data
	-mkdir ../output_nii

cleanup:
	-$(RM) ../output_nii/*.nii

# $(call a2d_convert,<inputdir>,<fullniftifile>)
define a2d_convert
echo "Convert Dicom to Nifti (mrconvert)";\
if [ ! -d $1/tmp ] ;then \
([ -f $2 ] && $(RM) $2; $(MRCONVERT) $1 $2 > /dev/null) \
else echo "Cannot use mrconvert with tmp files."; touch $2 fi
endef

# $(call a2d_check,<inputdir>,<fullniftifile>)
define a2d_check
echo "Agilent2Dicom Check DCM and NII"; \
[ -f $2 ] && fslhd $2 | grep -e '^dim[1234]' | \
awk 'BEGIN{i=0} {a[i]=$$2; ++i} END{printf("NIFTI  Dimensions:        "); for (i=0;i<4;++i){printf("%d ",a[i]);if (i!=3){printf("x ")}}; printf("\n")}' \
printf "DICOM";mrinfo ../output_data/standard2d 2> /dev/null | grep 'Dimensions'
endef

# $(call a2d_verify,<dcmfile>)
define a2d_verify
echo "Verify Dicom (dciodvfy)";\
dciodvfy $1 2>&1 >/dev/null | grep -e '^Error' -e '^Warning'
endef


run_standard2d:
	$(FDF2DCM) -i ~/Monash016/amanda/ExampleAgilentData/standard2d/ -o ../output_data/standard2d
	-$(call a2d_convert,../output_data/standard2d/,../output_nii/standard2d.nii)


check_standard2d:
	-$(call a2d_convert,../output_data/standard2d/,../output_nii/standard2d.nii)


test_standard2d:
	$(call a2d_verify,../output_data/standard2d/0001.dcm)

view_standard2d: 
	$(MRVIEW) ../output_data/standard2d &

viewn_standard2d: 
	$(FSLVIEW) ../output_nii/standard2d.nii &

.PHONY: standard2d test_standard2d run_standard2d check_standard2d view_standard2d
standard2d: run_standard2d check_standard2d test_standard2d


## Multi-echo 3D
run_me3d:
	$(FDF2DCM) -v -i  ~/Monash016/amanda/ExampleAgilentData/multiecho3d_magonly/ -o ../output_data/multiecho3d_magonly
	-$(call a2d_convert,../output_data/multiecho3d_magonly,../output_nii/ME3d.nii)
# [ ! -d ../output_data/multiecho3d_magonly/tmp ] && $(MRCONVERT)   ../output_data/multiecho3d_magonly/ ../output_nii/ME3d.nii > /dev/null 
#	$(FSLVIEW) ../output_nii/ME3d.nii

check_me3d:
	-$(call a2d_check,../output_data/multiecho3d_magonly,../output_nii/ME3d.nii)

test_me3d:
	-$(call a2d_verify,../output_data/multiecho3d_magonly/0001.dcm 2>&1)

view_me3d: 
	$(MRVIEW) ../output_data/multiecho3d_magonly &

viewn_me3d: 
	$(FSLVIEW) ../output_nii/ME3d.nii &

.PHONY: me3d run_me3d check_me3d view_me3d test_me3d
me3d: run_me3d check_me3d test_me3d

## Multi-echo 2D Mag and phase
run_me2d:
	$(FDF2DCM) -v -i ~/Monash016/amanda/ExampleAgilentData/multiecho2d_magandphase/ -o  ../output_data/multiecho2d_magandphase
	-$(call a2d_convert,../output_data/multiecho2d_magandphase/magnitude.dcm,../output_nii/ME2d_mag.nii)
	-$(call a2d_convert,../output_data/multiecho2d_magandphase/phase.dcm,../output_nii/ME2d_phase.nii)

check_me2d:
	-$(call a2d_check,../output_data/multiecho2d_magandphase/magnitude.dcm,../output_nii/ME2d_mag.nii)
	-$(call a2d_check,../output_data/multiecho2d_magandphase/phase.dcm,../output_nii/ME2d_phase.nii)


test_me2d:
	-$(call a2d_verify,../output_data/multiecho2d_magandphase/magnitude.dcm/0001.dcm)
	-$(call a2d_verify,../output_data/multiecho2d_magandphase/phase.dcm/0001.dcm)

view_me2d: 
	$(MRVIEW) ../output_data/multiecho2d_magandphase/magnitude.dcm &
	$(MRVIEW) ../output_data/multiecho2d_magandphase/phase.dcm &

viewn_me2d: 
	$(FSLVIEW) ../output_nii/ME2d_mag.nii ../output_nii/ME2d_phase.nii 



.PHONY: me2d run_me2d check_me2d test_me2d view_me2d
me2d: run_me2d check_me2d test_me2d


## CINE
run_cine:
	$(FDF2DCM) -v -i ~/Monash016/amanda/ExampleAgilentData/cine -o ../output_data/cine
	$(call a2d_convert,../output_data/cine,../output_nii/CINE.nii)

check_cine:
	$(call a2d_check,../output_data/cine,../output_nii/CINE.nii)


test_cine:
	-$(call a2d_verify,../output_data/cine/0001.dcm)

view_cine: 
	$(MRVIEW) ../output_data/cine &

viewn_cine: 
	$(FSLVIEW) ../output_nii/CINE.nii &


.PHONY: cine run_cine check_cine test_cine view_cine
cine: run_cine check_cine test_cine


## ASL
run_asl:
	$(FDF2DCM) -v -i ../example_data/1008.2.40.4.1.1/ASL_se_06.img -o ../output_data/ASL_se_06.dcm
	$(call a2d_convert,../output_data/ASL_se_06.dcm,../output_nii/ASL.nii)

check_asl:
	$(call a2d_check,../output_data/ASL_se_06.dcm,../output_nii/ASL.nii)


test_asl:
	-$(call a2d_verify,../output_data/ASL_se_06.dcm/0001.dcm)

view_asl: 
	$(MRVIEW) ../output_data/ASL_se_06.dcm &

viewn_asl: 
	$(FSLVIEW) ../output_nii/ASL.nii &

.PHONY: asl run_asl check_asl test_asl view_asl
asl: run_asl check_asl test_asl



## Diffusion
run_diffusion:
	$(FDF2DCM) -v -i ~/Monash016/amanda/ExampleAgilentData/diffusion/ -o ../output_data/diffusion
	$(call a2d_convert,../output_data/diffusion,../output_nii/Diffusion.nii)


check_diffusion:
	$(call a2d_check,../output_data/diffusion,../output_nii/Diffusion.nii)

test_diffusion:
	-$(call a2d_verify,../output_data/diffusion/0001.dcm)
	-$(call a2d_verify,../output_data/diffusion/0002.dcm)

view_diffusion: 
	$(MRVIEW) ../output_data/diffusion &

viewn_diffusion: 
	$(FSLVIEW) ../output_nii/Diffusion.nii &


.PHONY: diffusion run_diffusion check_diffusion test_diffusion view_diffusion
diffusion: run_diffusion check_diffusion test_diffusion

## Kidney 3D MGEMs

run_kidney: check
	$(FDF2DCM) -v -i ~/Monash016/RatKidney/Agilent/20120522/kidney512iso_01.img -o ../output_data/kidney512iso.dcm
	$(call a2d_convert,../output_data/kidney512iso.dcm/,../output_nii/Kidney512iso.nii)


check_kidney:
	$(call a2d_check,../output_data/kidney512iso.dcm/,../output_nii/Kidney512iso.nii)

test_kidney:
	-$(call a2d_verify,../output_data/kidney512iso.dcm/0001.dcm)

view_kidney: 
	$(MRVIEW) ../output_data/kidney512iso &

viewn_kidney: 
	$(FSLVIEW) ../output_nii/Kidney512iso.nii &


.PHONY: kidney run_kidney check_kidney test_kidney view_kidney
kidney: run_kidney check_kidney test_kidney


## Diffusion DTI
run_dti:
	$(FDF2DCM) -v -i ../example_data/s_2014051208/DTI_EPIP_J30_01.img/ -o ../output_data/dti/
	$(call a2d_convert,../output_data/dti/,../output_nii/DTI.nii)

check_dti:
	$(call a2d_check,../output_data/dti/,../output_nii/DTI.nii)

test_dti:
	$(call a2d_verify,../output_data/dti/0001.dcm)

view_dti: 
	$(MRVIEW) ../output_data/dti &

viewn_dti: 
	$(FSLVIEW) ../output_nii/DTI.nii &


.PHONY: dti run_dti check_dti view_dti test_dti
dti: run_dti check_dti test_dti





## Diffusion EPIP 
run_epip:
	$(FDF2DCM) -v -i ../example_data/s_2014051208/epip_axi_300TR_01.img/ -o ../output_data/epip/
	$(call a2d_convert,../output_data/epip/,../output_nii/EPI.nii)


check_epip:
	$(call a2d_check,../output_data/epip/,../output_nii/EPI.nii)


test_epip:
	$(call a2d_verify,../output_data/epip/0001.dcm)

view_epip: 
	$(MRVIEW) ../output_data/epip &


viewn_epip: 
	$(FSLVIEW) ../output_nii/EPI.nii &

.PHONY: epip run_epip check_epip view_epip test_epip
epip: run_epip check_epip test_epip


## Diffusion DTI
run_J6:
	$(FDF2DCM) -v -i ../example_data/s_2014051901/J6-500_01.img/ -o ../output_data/j6/
	$(call a2d_convert,../output_data/j6/,../output_nii/J6.nii)


check_J6: 
	$(call a2d_check,../output_data/j6/,../output_nii/J6.nii)

test_J6:
	$(call a2d_verify,../output_data/j6/0001.dcm)
	$(call a2d_verify,../output_data/j6/0002.dcm)

view_J6: run_J6
	$(MRVIEW) ../output_data/j6 &

viewn_J6: run_J6
	$(FSLVIEW) ../output_nii/J6.nii &

.PHONY: J6 run_J6 check_J6 view_J6 test_J6
J6: run_J6 check_J6 test_J6



## Diffusion EPI heart tissue - external recon
heart:
	$(FDF2DCM) -i ../example_data/s_2014061202/epi-dir30_01.img/ -o ../atrialtissue/dcm/
	$(call a2d_convert,../atrialtissue/dcm/,../atrialtissue/nifti/heart.nii)
	-$(RM) ../atrialtissue/mrtrix/diff.mif
	$(MRCONVERT) ../atrialtissue/dcm/ ../atrialtissue/mrtrix/diff.mif
	-$(RM) ../atrialtissue/mrtrix/dt_DWI.mif
	dwi2tensor ../atrialtissue/mrtrix/diff.mif ../atrialtissue/mrtrix/dt_DWI.mif

check_heart:
	$(call a2d_check,../atrialtissue/dcm/,../atrialtissue/nifti/heart.nii)

test_heart:
	$(call a2d_verify,../atrialtissue/dcm/0001.dcm)
	$(call a2d_verify,../atrialtissue/dcm/0002.dcm)

view_heart:
	$(MRVIEW) ../atrialtissue/dcm &

viewn_heart:
	$(FSLVIEW) ../atrialtissue/nifti/heart.nii &

.PHONY: heart run_J6 check_heart test_heart view_heart
heart: run_J6 check_heart test_heart

.PHONY: all 
all: standard2d me3d me2d cine asl diffusion dti epip J6 kidney heart

.PHONY: test_all 
test_all: test_kidney test_standard2d test_me3d test_me2d test_cine test_asl test_diffusion test_dti test_epip test_J6

.PHONY:  check_all
check_all: check_standard2d check_me3d check_me2d check_cine check_asl check_diffusion check_kidney check_epip

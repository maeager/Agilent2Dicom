#!/usr/bin/env bash
# Fixes to the MR Diffusion macro from dcmulti's conversion of 2D slices to enhanced MR format
# - Michael Eager (michael.eager@monash.edu)
# - Monash Biomedical Imaging 
# - (C) 2014 Michael Eager

# Check DCMTK on MASSIVE or Agilent console
if test ${MASSIVEUSERNAME+defined}; then
    if [ ! -x `which dcmodify` ];then
	module load dcmtk
    fi
else
    DCMTK="/home/vnmr1/src/dcmtk-3.6.0/bin"
    export PATH=${PATH}:${DCMTK}

fi

if [ ! -x `which dcmodify` ];then
    echo "ERROR: dcmodify not found (fix-dicoms.sh)"; 
    exit 1
fi


output_dir=$1
MODIFY=1
##COMMON FIXES to enhanced DICOMs
DCMODIFY="dcmodify --no-backup  " # --ignore-errors" 
files=$(find ${output_dir} -type f -name "*.dcm" | grep -v tmp)

nbdirs=$(ls -1 output_data/hearttissue.dcm/tmp/slice* | sed 's/.*image0\(.*\)echo.*/\1/' | tail -1)

for idx in $(seq 1 ${nbdirs}); do
    # Grab parameters from first slice in image
    dcmdump ${output_dir}/tmp/slice001image$(printf '%03d' $idx)echo001.dcm 2>&1  | grep 'Diff' |grep -v Sequence| awk '{print $3}' | tr -d '[]' > ${output_dir}/diffusion.tmp
# dcmdump output_data/hearttissue.dcm/tmp/slice001image002echo001.dcm 2>&1 | grep Diff | awk '{print $1,$3,$NF}'
    if [ -f  ${output_dir}/diffusion.tmp ]; then
	declare -a DiffusionSequence
	let i=0;while IFS=$'\r\n' read -r line_data; do 
	    DiffusionSequence[i]="${line_data}"; ((++i)); 
	done < ${output_dir}/diffusion.tmp
    	
	if [ ${#DiffusionSequence[*]} -eq 3 ]; then

	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9075)=NONE" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9087)=0" ${output_dir}/$(printf '%04d' $idx).dcm
	    
	else
	    echo "Diffusion sequence params: ",  ${DiffusionSequence[*]}
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9075)=${DiffusionSequence[0]}" ${output_dir}/$(printf '%04d' $idx).dcm  ## Directionality
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9076)[0].(0018,9089)=${DiffusionSequence[1]}" ${output_dir}/$(printf '%04d' $idx).dcm  ## Gradient Orientation
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9087)=${DiffusionSequence[2]}" ${output_dir}/$(printf '%04d' $idx).dcm  ## BValue
## B MATRIX
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9076)[0].(0018,9602)=${DiffusionSequence[3]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9076)[0].(0018,9603)=${DiffusionSequence[4]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9076)[0].(0018,9604)=${DiffusionSequence[5]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9076)[0].(0018,9605)=${DiffusionSequence[6]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9076)[0].(0018,9606)=${DiffusionSequence[7]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(0018,9117)[0].(0018,9076)[0].(0018,9607)=${DiffusionSequence[8]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    
	fi
    fi
    ((++index))
done

rm -f ${output_dir}/diffusion.tmp


# done

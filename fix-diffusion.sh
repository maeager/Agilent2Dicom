#!/usr/bin/env bash
## Fixes bug in the MR Diffusion macro from dcmulti's conversion 
# of 2D slices to enhanced MR format
#
# - Michael Eager (michael.eager@monash.edu)
# - Monash Biomedical Imaging 
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


# Check DCMTK on MASSIVE or Agilent console
if test ${MASSIVEUSERNAME+defined}; then
    if [ ! -x $(which dcmodify) ];then
	module load dcmtk
    fi
else
    DCMTK="/home/vnmr1/src/dcmtk-3.6.0/bin"
    export PATH=${PATH}:${DCMTK}
fi

if [ ! -x $(which dcmodify) ];then
    error_exit("dcmodify not found (fix-dicoms.sh)"); 
    exit 1
fi


output_dir=$1
MODIFY=1
## COMMON FIXES to enhanced DICOMs
DCMODIFY="dcmodify --no-backup  " # --ignore-errors" 
files=$(find ${output_dir} -type f -name "*.dcm" | grep -v tmp)

nbdirs=$(ls -1 ${output_dir}/tmp/slice* | sed 's/.*image0\(.*\)echo.*/\1/' | tail -1)

for idx in $(seq 1 ${nbdirs}); do
    # Grab parameters from first slice in image
    dcmdump ${output_dir}/tmp/slice001image$(printf '%03d' $idx)echo001.dcm 2>&1  | grep 'Diff' |grep -v Sequence| awk '{print $3}' | tr -d '[]' > ${output_dir}/diffusion.tmp
# dcmdump output_data/hearttissue.dcm/tmp/slice001image002echo001.dcm 2>&1 | grep Diff | awk '{print $1,$3,$NF}'
    if [ -f  ${output_dir}/diffusion.tmp ]; then
	declare -a DiffusionSequence
	let i=0;while IFS=$'\r\n' read -r line_data; do 
	    DiffusionSequence[i]="${line_data}"; ((++i)); 
	done < ${output_dir}/diffusion.tmp
    	echo "Diffusion sequence params: ",  ${DiffusionSequence[*]}
	if [ ${#DiffusionSequence[*]} -lt 10 ]; then ## BVal index differs -o ${DiffusionSequence[2]} -lt 10
	    ## (Shared sequence)[0].(MR Diffusion Sequence)[0]
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9075)=NONE" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9087)=0" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9147)=${DiffusionSequence[2]}" ${output_dir}/$(printf '%04d' $idx).dcm  
	else
	    
	    ## (Shared sequence)[0].(MR Diffusion Sequence)[0](Directionality)
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9075)=${DiffusionSequence[0]}" ${output_dir}/$(printf '%04d' $idx).dcm  ## 
	    ## (Shared sequence)[0].(MR Diffusion Sequence)[0](Diff Gradient Orientation)
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9076)[0].(0018,9089)=${DiffusionSequence[1]}" ${output_dir}/$(printf '%04d' $idx).dcm 
            ## (Shared sequence)[0].(MR Diffusion Sequence)[0](BValue)
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9087)=${DiffusionSequence[2]}" ${output_dir}/$(printf '%04d' $idx).dcm  
            ## (Shared sequence)[0].(MR Diffusion Sequence)[0](Anisotrophy)
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9147)=${DiffusionSequence[3]}" ${output_dir}/$(printf '%04d' $idx).dcm  

## (Shared sequence)(MR Diffusion Sequence)(B MATRIX sequence)( B val re[XYZ][XYZ])
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9076)[0].(0018,9601)[0].(0018,9602)=${DiffusionSequence[4]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9076)[0].(0018,9601)[0].(0018,9603)=${DiffusionSequence[5]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9076)[0].(0018,9601)[0].(0018,9604)=${DiffusionSequence[6]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9076)[0].(0018,9601)[0].(0018,9605)=${DiffusionSequence[7]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9076)[0].(0018,9601)[0].(0018,9606)=${DiffusionSequence[8]}" ${output_dir}/$(printf '%04d' $idx).dcm
	    ${DCMODIFY} -i "(5200,9229)[0].(0018,9117)[0].(0018,9076)[0].(0018,9601)[0].(0018,9607)=${DiffusionSequence[9]}" ${output_dir}/$(printf '%04d' $idx).dcm
	fi
    fi
done

rm -f ${output_dir}/diffusion.tmp


# done

#!/usr/bin/env  bash 
## FID/FDF DICOM converter functions
#
# - Michael Eager (michael.eager@monash.edu.au)
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
#


function get_output_dirs(){
    output_root=$(dirname $output_dir)
    output_base=$(basename $output_dir .dcm)
    dirs=$(find $output_root -maxdepth 1 -type d  -name  "$output_base*.dcm")
    echo "fid2dicom.py completed successfully. Dicom paths generated were: "
    echo $dirs
}


function yesno(){
    read -r -p "$@" response
    response=$(echo $response | awk '{print tolower($0)}')
# response=${response,,} # tolower Bash 4.0
    if [[ $response =~ ^(yes|y| ) ]]; then
	return 0
    fi
    return 1
}

function error_exit(){
    echo "${PROGNAME}: ${1:-Unknown error}" 1>&2
    exit 1
}



function move_dcm_to_tmp_fid(){    
    echo "Moving dicom images"
    for dcmdir in ${dirs}; do
	mkdir ${dcmdir}/tmp
	mv "${dcmdir}"/*.dcm "${dcmdir}"/tmp/
    done
}

function enhance_dicoms_fdf(){
    echo "Convert dicom images to single enhanced MR dicom format image"
    if [ -f ${output_dir}/MULTIECHO ]
    then
	echo "Contents of MULTIECHO"; cat ${output_dir}/MULTIECHO; echo '\n'
	nechos=$(cat ${output_dir}/MULTIECHO)
	nechos=$(printf "%1.0f" $nechos)
	echo "Multi echo sequence, $nechos echos"
	for ((iecho=1;iecho<=nechos;++iecho)); do
     	    echoext=$(printf '%03d' $iecho)
     	    echo "Converting echo ${iecho} using dcmulti"
     	    ${DCMULTI} "${output_dir}/0${echoext}.dcm" $(ls -1 ${output_dir}/tmp/*echo${echoext}.dcm | \
		sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | \
		sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
	done

# DCMULTI="dcmulti -v -makestack -sortby EchoTime -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
#-makestack -sortby ImagePositionPatient  -sortby AcquisitionNumber
#  ${DCMULTI} ${output_dir}/0001.dcm $(ls -1 ${output_dir}/tmp/*.dcm  | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')

	${RM} ${output_dir}/MULTIECHO
	echo "Multi echo sequence completed."
	
    elif  [ -f ${output_dir}/DIFFUSION ]; then

	echo "Contents of DIFFUSION"; cat ${output_dir}/DIFFUSION; echo '\n'

    # nbdirs=$(cat ${output_dir}/DIFFUSION)
    # ((++nbdirs)) # increment by one for B0
	nbdirs=$(ls -1 ${output_dir}/tmp/slice* | sed 's/.*image0\(.*\)echo.*/\1/' | tail -1)

	echo "Diffusion sequence, $nbdirs B-directions"
	for ((ibdir=1;ibdir<=nbdirs;ibdir++)); do
     	    bdirext=$(printf '%03d' $ibdir)

     	    echo "Converting bdir ${ibdir} using dcmulti"

	## Input files are sorted by image number and slice number. 
     	    ${DCMULTI} "${output_dir}/0${bdirext}.dcm" $(ls -1 ${output_dir}/tmp/*image${bdirext}*.dcm | \
		sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | \
		sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')

	done
	echo "Diffusion files compacted."

    elif  [ -f ${output_dir}/ASL ]; then

	echo "Contents of ASL"; cat ${output_dir}/ASL; echo '\n'

    # nbdirs=$(cat ${output_dir}/ASL)
    # ((++nbdirs)) # increment by one for B0
	asltags=$(ls -1 ${output_dir}/tmp/slice* | sed 's/.*image0\(.*\)echo.*/\1/' | tail -1)

	echo "ASL sequence"
	for ((iasl=1;iasl<=2;iasl++)); do
     	    aslext=$(printf '%03d' $iasl)

     	    echo "Converting ASL tag ${iasl} using dcmulti"

	## Input files are sorted by image number and slice number. 
     	    ${DCMULTI} "${output_dir}/0${aslext}.dcm" $(ls -1 ${output_dir}/tmp/*echo${aslext}.dcm | \
		sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | \
		sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')

	done
	echo "ASL files converted."


    else

    ## Dcmulti config is dependent on order of files.  The 2D standard
    ## dicoms are sorted by echo time, image number then slice
    ## number. The second argument reorders the list of 2D dicom files
    ## based on echo time, then image number, then slice number.
    ## Only one output file is required, 0001.dcm. 
	${DCMULTI} ${output_dir}/0001.dcm $(ls -1 ${output_dir}/tmp/*.dcm  | \
	    sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | \
	    sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
	[ $? -ne 0 ] && error_exit "$LINENO: dcmulti failed"

    fi
    echo "DCMULTI complete. Fixing inconsistencies."

}

function enhance_dicoms_fid(){
    
    echo "Convert dicom images to single enhanced MR dicom format image"
    if [ -f ${output_dir}/MULTIECHO ]
    then
	echo "Contents of MULTIECHO"; cat ${output_dir}/MULTIECHO; echo '\n'
	nechos=$(cat ${output_dir}/MULTIECHO)
	nechos=$(printf "%1.0f" $nechos)
	echo "Multi echo sequence, $nechos echos"
	for ((iecho=1;iecho<=nechos;++iecho)); do
     	    echoext=$(printf '%03d' $iecho)
     	    echo "Converting echo ${iecho} using dcmulti"
     	    for dcmdir in $dirs; do
		${DCMULTI} "${dcmdir}/0${echoext}.dcm" $(ls -1 ${dcmdir}/tmp/*echo${echoext}.dcm | \
		    sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | \
		    sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
	    done
	done

# DCMULTI="dcmulti -v -makestack -sortby EchoTime -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
#-makestack -sortby ImagePositionPatient  -sortby AcquisitionNumber
#  ${DCMULTI} ${output_dir}/0001.dcm $(ls -1 ${output_dir}/tmp/*.dcm  | \
#     sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | \
#     sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')

	${RM} ${output_dir}/MULTIECHO
	echo "Multi echo sequence completed."
	
    elif  [ -f ${output_dir}/DIFFUSION ]; then

	echo "Contents of DIFFUSION"; cat ${output_dir}/DIFFUSION; echo '\n'

    # nbdirs=$(cat ${output_dir}/DIFFUSION)
    # ((++nbdirs)) # increment by one for B0
	nbdirs=$(ls -1 ${output_dir}/tmp/slice* | sed 's/.*image0\(.*\)echo.*/\1/' | tail -1)

	echo "Diffusion sequence, $nbdirs B-directions"
	for ((ibdir=1;ibdir<=nbdirs;ibdir++)); do
     	    bdirext=$(printf '%03d' $ibdir)

     	    echo "Converting bdir ${ibdir} using dcmulti"
	    for dcmdir in $dirs; do
	## Input files are sorted by image number and slice number. 
     		${DCMULTI} "${dcmdir}/0${bdirext}.dcm" $(ls -1 ${dcmdir}/tmp/*image${bdirext}*.dcm | \
		    sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | \
		    sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
	    done
	done
	echo "Diffusion files compacted."

    elif  [ -f ${output_dir}/ASL ]; then

	echo "Contents of ASL"; cat ${output_dir}/ASL; echo '\n'

    # nbdirs=$(cat ${output_dir}/ASL)
    # ((++nbdirs)) # increment by one for B0
	asltags=$(ls -1 ${output_dir}/tmp/slice* | sed 's/.*image0\(.*\)echo.*/\1/' | tail -1)

	echo "ASL sequence"
	for ((iasl=1;iasl<=2;iasl++)); do
     	    aslext=$(printf '%03d' $iasl)

     	    echo "Converting ASL tag ${iasl} using dcmulti"
	    for dcmdir in $dirs; do
	## Input files are sorted by image number and slice number. 
     		${DCMULTI} "${dcmdir}/0${aslext}.dcm" $(ls -1 ${dcmdir}/tmp/*echo${aslext}.dcm | \
		    sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | \
		    sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
	    done
	done
	echo "ASL files converted."


    else

    ## Dcmulti config is dependent on order of files.  The 2D standard
    ## dicoms are sorted by echo time, image number then slice
    ## number. The second argument reorders the list of 2D dicom files
    ## based on echo time, then image number, then slice number.
    ## Only one output file is required, 0001.dcm. 
	for dcmdir in $dirs; do
	    ${DCMULTI} ${dcmdir}/0001.dcm $(ls -1 ${dcmdir}/tmp/*.dcm  | \
		sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | \
		sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
	    [ $? -ne 0 ] && error_exit "$LINENO: dcmulti failed"
	done
    fi
    echo "DCMULTI complete. Fixing inconsistencies."
}

function fix_enhanced_dicoms_fid(){
## Corrections to dcmulti conversion
    for dcmdir in $dirs; do
	${FID2DCMPATH}/fix-dicoms.sh "${dcmdir}"
	echo "Fixing dicom dir ${dcmdir}"

    ## Additional corrections to diffusion files
	if [ -f ${output_dir}/DIFFUSION ];then
	    ${FID2DCMPATH}/fix-diffusion.sh "${dcmdir}"
	    echo "Fixed diffusion module parameters."
	fi
    ## Additional corrections to ASL files
	if [ -f ${output_dir}/ASL ];then
	    ${FID2DCMPATH}/fix_asl.sh "${dcmdir}"
	    echo "Fixed ASL module parameters."
	fi
    done
    
    [ -f ${output_dir}/DIFFUSION ] && ${RM} ${output_dir}/DIFFUSION
    [ -f ${output_dir}/ASL ] && ${RM} ${output_dir}/ASL
}

function iodvfy(){
    if (( verbosity > 0 )); then
	echo "Verifying dicom compliance using dciodvfy."
	if [ -f "${output_dir}/0001.dcm" ]; then
	    set +e
	## Send dciodvfy stderr and stdout to log file
	    dciodvfy "${output_dir}/0001.dcm" &> $(dirname ${output_dir})/$(basename ${output_dir} .dcm).log
	    set -e  
	else
	    error_exit "$LINENO: could not find ${output_dir}/0001.dcm for verification"
	fi
    fi
}

function clean_tmps_fid(){
## Cleaning up temporary directories
    (( verbosity > 0 )) && echo "Cleaning up tmp directories."
# Check the raw output only, then delete all tmp 
    if [ -d "${output_dir}/tmp" ]
    then
	if (( verbosity > 0 ))
	then
	    if yesno "Remove existing tmp output directories, y or n (default y)?"
	    then
		echo "Removing existing tmp output directory"
		for dcmdir in $dirs; do
		    ${RMDIR} "${dcmdir}/tmp"    
		done
	    else
		echo "fid2dcm completed. Temporary dicoms still remain."
		exit 0
	    fi
	else
	    echo "Removing existing tmp output directory"
	    for dcmdir in $dirs; do
		${RMDIR} "${dcmdir}/tmp"    
	    done
	fi
	[ -d "${output_dir}/tmp" ] && error_exit "$LINENO: temporary dicom directory could not be deleted."
    fi
}

#!/usr/bin/env  bash 
## FDF to DICOM converter
#   Front-end to agilent2dicom and dcmulti
#
# - Michael Eager (michael.eager@monash.edu.au)
# - Monash Biomedical Imaging 
# - (C) 2014 Michael Eager
#
#
#  "$Id:"
#  Version 0.0: Wrapper for agilent2dicom
#  Version 1.0: Full support of FDF formats
###############################################

set -o nounset  # shortform: -u
set -o errexit  # -e
# set -o pipefail
# touch $(dirname $0)/error.log
#exec 2>> $(dirname $0)/error.log
#set -x  # show debugging output


## Set config variables
FDF2DCMVERSION=1.1
PROGNAME=$(basename $0)
FDF2DCMPATH=$(dirname $0)
KERNEL_RELEASE=$(uname -r | awk -F'.' '{printf("%d.%d.%d\n", $1,$2,$3)}')
DCM3TOOLS="${FDF2DCMPATH}/../dicom3tools_1.00.snapshot.20140306142442/bin/1.${KERNEL_RELEASE}.x8664/"
## Set dcmulti and its arguments
DCMULTI="dcmulti -v -makestack -sortby AcquisitionNumber -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
DCMULTI_DTI="dcmulti -v -makestack -sortby DiffusionBValue -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
#-makestack -sortby ImagePositionPatient  -sortby AcquisitionNumber
# Check dcmtk applications on MASSIVE or Agilent console
if test ${MASSIVEUSERNAME+defined}; then
    test -x dcmodify || module load dcmtk 
else
    DCMTK="/home/vnmr1/src/dcmtk-3.6.0/bin"
    export PATH=${PATH}:${DCMTK}
fi
if [ ! -d ${DCM3TOOLS} ]; then
    echo "${DCM3TOOLS} path not found"
    exit 1
elif [ ! -x ${DCM3TOOLS}/dcmulti ]; then
    echo "Unable to find Dicom3Tool's executable dcmulti"
    exit 1
fi 
export PATH=${PATH}:${DCM3TOOLS}
declare -i verbosity=0
declare -i do_modify=1
declare input_dir=""
declare output_dir=""
declare python_args=""

E_BADARGS=65

#source ${FDF2DCMPATH}/yesno.sh
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


# Print usage information and exit
print_usage(){
    echo -e "\n" \
	"usage: ./fdf2dcm.sh -i inputdir [-o outputdir] [-v] [-m] [-p]\n" \
	"\n" \
	"-i <inputdir>  FDF source directory\n" \
	"-o <outputdir> Optional destination DICOM directory. Default is input_dir/.dcm. \n" \
	"-d             Disable dcmodify fixes to DICOMs." \
	"-m,-p          Enable magnitude and phase subdirectory conversion.  These flags are passed to agilent2dicom and should only be used from within fdf2dcm or with knowledge of input fdf data. \n" \
	"-v             verbose output. \n" \
	"-h             this help\n" \
	"\n" 
    # && exit 1
}


## Check for number of args
if [ $# -eq 0 ]; then
    echo "fdfdcm.sh must have one argument: -i, --input [directory of FDF images]"
    print_usage
    exit $E_BADARGS
fi


## Parse arguments
while getopts ":i:o:hmpdv" opt; do
    case $opt in
	i)
	    echo "Input dir:  $OPTARG" >&2
	    input_dir="$OPTARG"
	    ;;
	o)
	    echo "Output dir: $OPTARG" >&2
	    output_dir="$OPTARG"
	    ;;
	h)
	    print_usage
	    ${FDF2DCMPATH}/agilent2dicom.py -h
	    exit 0
	    ;;
	m)
	    echo "Implementing magnitude component of FDF to DICOM conversion."
	    python_args="$python_args -m"
	    ;;
	p)
	    echo "Implementing phase component of FDF to DICOM conversion."
	    python_args="$python_args -p"
	    ;;
	d)
	    do_modify=0
	    echo " Disable dcmodify correction."
	    ;;
	v)
	    verbosity=1
	    echo "Setting verbose to $verbosity."
	    python_args="$python_args -v"
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    print_usage
	    exit $E_BADARGS
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    print_usage
	    exit $E_BADARGS
	    ;;
    esac
done


# Clean up input args
if [ ! -d "$input_dir" ]; then
    echo "fdfdcm.sh must have a valid input directory of FDF images."
    exit $E_BADARGS
fi
## Set output_dir if not in args, default is INPUT/.dcm
if [ "$output_dir" == "" ]; then
    output_dir="$(dirname ${input_dir})/$(basename ${input_dir} .img).dcm"
    echo "Output dir set to: " $output_dir
fi

## Check for MAG/PHASE component directories
if [ -d  ${input_dir}/magnitude.img ] || [ -d ${input_dir}/phase.img ]; then
    echo "Input directory has 'magnitude.img' and 'phase.img' "
    verb=''; if [ "$verbosity" -eq 1 ]; then verb=' -v '; fi 
    $0 $verb -m -i ${input_dir}/magnitude.img/ -o ${output_dir}/magnitude.dcm/
    if [ "$verbosity" -eq 1 ] && ! yesno "Magnitude complete. Continue converting phase?"; then
	echo "fdf2dcm completed without phase." && exit 0
    fi	
    $0 $verb -p -i ${input_dir}/phase.img/ -o ${output_dir}/phase.dcm/
    exit 0
fi

## Check output directory
JumpToDCmulti=0
if [ -d "${output_dir}" ]; then
    if [ -d "${output_dir}/tmp" -a "$verbosity" -eq 1 ];then
	if yesno "Remove existing output directory, y or n (default y)?"; then
	    echo "Removing existing output directory"
	    rm -rf ${output_dir}
	else
	    JumpToDCmulti=1
	fi
    else
	echo "Removing existing output directory"
	rm -rf ${output_dir}
    fi	
fi

if (( "${JumpToDCmulti}" == 0 ))
then
    shopt -s nullglob  # Carefull with nullglob on csh  
    found=0
    for i in ${input_dir}/*.fdf; do
	((++found))
    done
    shopt -u nullglob
    if [ $found -eq 0 ]; then  #-o "$fdffiles" == "" 
	error_exit "$LINENO: Input directory has no FDF images"
    else
	echo $found, " FDF files were found"
    fi

    if [ ! -f ${input_dir}/procpar ]; then
	error_exit "$LINENO: Input directory has no procpar file"
    fi

# set -o errexit  # -e
# set -o pipefail


## Crux of script - conversion of FDF images to standard DICOM images
    echo  "Calling agilent2dicom"
    echo " Arguments: ", ${python_args} -i "${input_dir}" -o "${output_dir}"
    ${FDF2DCMPATH}/agilent2dicom.py ${python_args} -i "${input_dir}" -o "${output_dir}"

    [ $? -ne 0 ] && error_exit "$LINENO: agilent2dicom failed"
    
    [ ! -d ${output_dir} ] && error_exit "$LINENO: Output dir not created by agilent2dicom."

    # dcmfiles=$(ls ${output_dir}/*.dcm)  ## Bad code - use glob
    #if[ $? -ne 0 ]
    test -e ${output_dir}/0001.dcm  &&	error_exit "$LINENO: Output directory of agilent2dicom has no DICOM images."
    

    echo "Moving dicom images"
    mkdir ${output_dir}/tmp
    mv ${output_dir}/*.dcm ${output_dir}/tmp/

fi ## JumpToDCMulti

echo "Convert dicom images to single enhanced MR dicom format image"
if [ -f ${output_dir}/MULTIECHO ]; then
    echo "Contents of MULTIECHO"; cat ${output_dir}/MULTIECHO; echo '\n'
    nechos=$(cat ${output_dir}/MULTIECHO)
    echo "Multi echo sequence, $nechos echos"
    for iecho in $(seq 1 ${nechos}); do
     	echoext=$(printf '%03d' $iecho)
     	echo "Converting echo ${iecho} using dcmulti"
     	${DCMULTI} "${output_dir}/0${echoext}.dcm" $(ls -1 ${output_dir}/tmp/*echo${echoext}.dcm | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
    done

# DCMULTI="dcmulti -v -makestack -sortby EchoTime -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
#-makestack -sortby ImagePositionPatient  -sortby AcquisitionNumber
#  ${DCMULTI} ${output_dir}/0001.dcm $(ls -1 ${output_dir}/tmp/*.dcm  | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')

    rm -f ${output_dir}/MULTIECHO
elif  [ -f ${output_dir}/DIFFUSION ]; then
    echo "Contents of DIFFUSION"; cat ${output_dir}/DIFFUSION; echo '\n'
    nbdirs=$(cat ${output_dir}/DIFFUSION)
    ((++nbdirs)) # increment by one for B0
    nbdirs=$(ls -1 ${output_dir}/tmp/slice* | sed 's/.*image0\(.*\)echo.*/\1/' | tail -1)
    echo "Diffusion sequence, $nbdirs B-directions"
    for ibdir in $(seq 1 ${nbdirs}); do
     	bdirext=$(printf '%03d' $ibdir)
     	echo "Converting bdir ${ibdir} using dcmulti"
     	${DCMULTI} "${output_dir}/0${bdirext}.dcm" $(ls -1 ${output_dir}/tmp/*image${bdirext}*.dcm | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
    done
    echo "Diffusion files compacted."

else

    # Dcmulti config is dependent on order of files.  The 2D standard dicoms are sorted by echo time, image number then slice number. 
    ${DCMULTI} ${output_dir}/0001.dcm $(ls -1 ${output_dir}/tmp/*.dcm  | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
fi
echo "DCMULTI complete. Fixing inconsistencies."


## Corrections to dcmulti conversion
if [ ${do_modify} -eq 1 ]
then
    ${FDF2DCMPATH}/fix-dicoms.sh "${output_dir}"
    echo "Fixing dicoms complete."
## Additional corrections to diffusion files
    if [ -f ${output_dir}/DIFFUSION ];then
	${FDF2DCMPATH}/fix-diffusion.sh "${output_dir}"
	echo "Fixed diffusion module parameters."
	rm -f ${output_dir}/DIFFUSION
    fi
fi
[ -f ${output_dir}/DIFFUSION ] && rm -f ${output_dir}/DIFFUSION

if [ "$verbosity" -eq 1 ]; then
    echo "Verifying dicom compliance using dciodvfy."
    if [ -f "${output_dir}"/0001.dcm ]; then
	## Send dciodvfy stderr and stdout to log file
	dciodvfy ${output_dir}/0001.dcm &> $(dirname ${output_dir})/$(basename ${output_dir} .dcm).log  
    else
	error_exit "$LINENO: could not find ${output_dir}/0001.dcm for verification"
    fi
fi

## Cleaning up temporary directories
if [ -d ${output_dir}/tmp ]; then
    if [ "$verbosity" -eq 1 ]; then
	if yesno "Remove existing tmp output directory, y or n (default y)?"; then
	    echo "Removing existing tmp output directory"
	    rm -rf ${output_dir}/tmp    
	else
	    echo "fdf2dcm completed. Temporary dicoms still remain."
	    exit 0
	fi
    else
	echo "Removing existing tmp output directory"
	rm -rf ${output_dir}/tmp    
    fi
    [ -d ${output_dir}/tmp ] && error_exit "$LINENO: temporary dicom directory could not be deleted."
fi


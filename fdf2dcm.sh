#!/usr/bin/env  bash 
## FDF to DICOM converter
#   Front-end to agilent2dicom and dcmulti
#
# - Michael Eager (michael.eager@monash.edu.au)
# - Monash Biomedical Imaging 
# - (C) 2014 Michael Eager


set -o nounset  # shortform: -u
# set -o errexit  # -e
# set -o pipefail
#exec 2>> $(dirname $0)/error.log
# set -x


# Set config variables
VERBOSE=0
MODIFY=1
FDF2DCMPATH=$(dirname $0)
KERNEL_RELEASE=$(uname -r | awk -F'.' '{printf("%d.%d.%d\n", $1,$2,$3)}')
DCM3TOOLS="${FDF2DCMPATH}/../dicom3tools_1.00.snapshot.20140306142442/bin/1.${KERNEL_RELEASE}.x8664/"
export PATH=${PATH}:${DCM3TOOLS}
input_dir=""
output_dir=""
python_args=""
# Check DCMTK on MASSIVE or Agilent console
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

E_BADARGS=65
source ${FDF2DCMPATH}/yesno.sh

DCMULTI="dcmulti -v -makestack -sortby AcquisitionNumber -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
#-makestack -sortby ImagePositionPatient  -sortby AcquisitionNumber



# Print usage information and exit
print_usage(){
    echo -e "\n" \
	"usage: ./fdf2dcm.sh -i inputdir [-o outputdir] [-v] [-m] [-p]\n" \
	"\n" \
	"-i <inputdir>  FDF source directory\n" \
	"-o <outputdir> Optional destination DICOM directory. Default is input_dir/.dcm. \n" \
	"-v             verbose output. \n" \
	"-m,-p          Enable magnitude and phase subdirectory conversion.  These flags are passed to agilent2dicom and should only be used from within fdf2dcm or with knowledge of input fdf data. \n" \
	"-h             this help\n" \
	"\n" 
    # && exit 1
}


# Check for number of args
if [ $# -eq 0 ]; then
    echo "fdfdcm.sh must have one argument: -i, --input [directory of FDF images]"
    print_usage
    exit $E_BADARGS
fi




# Parge arguments
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
	    MODIFY=0
	    echo " Disable dcmodify correction."
	    ;;
	v)
	    echo "Setting verbose to 1."
	    VERBOSE=1
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
# Set output_dir if not in args, default is INPUT/.dcm
if [ "$output_dir" == "" ]; then
    output_dir="$(dirname ${input_dir})/$(basename ${input_dir} .img).dcm"
    echo "Output dir set to: " $output_dir
fi
if [ -d "$output_dir" ]; then
    if [ "$VERBOSE" -eq 1 ];then
	if yesno "Remove existing output directory, y or n (default y)?"; then
	    echo "Removing existing output directory"
	    rm -rf $output_dir
	fi
    else
	rm -rf $output_dir
    fi	
fi

if [ -d  ${input_dir}/magnitude.img ] || [ -d ${input_dir}/phase.img ]; then
    echo "Input directory has 'magnitude.img' and 'phase.img' "
    verb=''; if [ "$VERBOSE" -eq 1 ]; then verb=' -v '; fi 
    $0 $verb -m -i ${input_dir}/magnitude.img/ -o ${output_dir}/magnitude.dcm/
    if [ "$VERBOSE" -eq 1 ] && ! yesno "Magnitude complete. Continue converting phase?"; then
	echo "fdf2dcm completed without phase." && exit 0
    fi	
    $0 $verb -p -i ${input_dir}/phase.img/ -o ${output_dir}/phase.dcm/
    exit 0
fi

# Carefull with nullglob on csh  
shopt -s nullglob
found=0
for i in ${input_dir}/*.fdf; do
    ((++found))
done
shopt -u nullglob
if [ $found -eq 0 ]; then  #-o "$fdffiles" == "" 
    echo "fdf2dcm Error: Input directory has no FDF images"
    exit 1
else
    echo $found, " FDF files were found"
fi

if [ ! -f ${input_dir}/procpar ]; then
    echo "fdf2dcm Error: Input directory has no procpar file"
    exit 1
fi

# set -o errexit  # -e
# set -o pipefail


## Crux of script
echo  "Calling agilent2dicom"
echo ${FDF2DCMPATH}/agilent2dicom.py $python_args -i "${input_dir}" -o "${output_dir}"
${FDF2DCMPATH}/agilent2dicom.py $python_args -i "${input_dir}" -o "${output_dir}"

if [ $? -ne 0 ]; then 
    echo "Error: agilent2dicom failed"
    exit 1
fi
[ ! -d $output_dir ] && (echo "Output dir not created by agilent2dicom. Shutting down fdf2dcm."; exit 1)
dcmfiles=`ls ${output_dir}/*.dcm`
if [ $? -ne 0 ]; then  #-o "$fdffiles" == "" 
    echo "Error: Output directory of agilent2dicom has no DICOM images. Shutting down fdf2dcm."
    exit 1
fi

echo "Moving dicom images"
mkdir ${output_dir}/tmp
mv ${output_dir}/*.dcm ${output_dir}/tmp/


echo "Convert dicom images to single enhanced MR dicom format image"
if [ -f ${output_dir}/MULTIECHO ]; then
    echo "Contents of MULTIECHO"; cat ${output_dir}/MULTIECHO; echo '\n'
    nechos=$(cat ${output_dir}/MULTIECHO)
    echo "Multi echo sequence, $nechos echos"
    # for iecho in $(seq 1 ${nechos}); do
    # 	echoext=$(printf '%03d' $iecho)
    # 	echo "Converting echo ${iecho} using dcmulti"
    # 	${DCMULTI} "${output_dir}/0${echoext}.dcm" $(ls -1 ${output_dir}/tmp/*echo${echoext}.dcm | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
    # done

DCMULTI="dcmulti -v -makestack -sortby EchoTime -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
#-makestack -sortby ImagePositionPatient  -sortby AcquisitionNumber
 ${DCMULTI} ${output_dir}/0001.dcm $(ls -1 ${output_dir}/tmp/*.dcm  | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')

    rm -f ${output_dir}/MULTIECHO
else

    # Dcmulti config is dependent on order of files.  The 2D standard dicoms are sorted by echo time, image number then slice number. 
    ${DCMULTI} ${output_dir}/0001.dcm $(ls -1 ${output_dir}/tmp/*.dcm  | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
fi


## Corrections to dcmulti conversion
${FDF2DCMPATH}/fix-dicoms.sh "${output_dir}"

if [ "$VERBOSE" -eq 1 ]; then
    echo "Verifying dicom compliance using dciodvfy."
    if [ -f "$output_dir"/0001.dcm ]; then
	dciodvfy $output_dir/0001.dcm &> $(dirname $output_dir)/$(basename $output_dir .dcm).log
    else
	echo "ERROR: could not find $output_dir/0001.dcm for verification"
    fi
fi

## Cleaning up
if [ -d ${output_dir}/tmp ]; then

    if [ "$VERBOSE" -eq 1 ] && ! (yesno "Remove existing tmp output directory, y or n (default y)?"); then
	echo "fdf2dcm completed. Temporary dicoms still remain." && exit 0
    fi

    echo "Removing existing tmp output directory"
    rm -rf ${output_dir}/tmp    

fi

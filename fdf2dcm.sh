#!/usr/bin/env  bash 
## FDF to DICOM converter
#   Front-end to agilent2dicom and dcmulti
#
# - Michael Eager (michael.eager@monash.edu.au)
# - (c) 2014

set -o nounset  # shortform: -u
# set -o errexit  # -e
# set -o pipefail

# Set config variables
VERBOSE=0
FDF2DCMPATH=$(dirname $0)
DCM3TOOLS="${FDF2DCMPATH}/../dicom3tools_1.00.snapshot.20140306142442/bin/1.2.6.35.x8664/"
if [ ! -d ${DCM3TOOLS} ]; then
   echo "${DCM3TOOLS} not found"
   exit 1
elif [ ! -f ${DCM3TOOLS}/dcmulti ]; then
    echo "Unable to find dcmulti"
    exit 1
fi 
export PATH=${PATH}:${DCM3TOOLS}
E_BADARGS=65
source ${FDF2DCMPATH}/yesno.sh

DCMULTI="dcmulti -v  -makestack -sortby ImagePositionPatient  -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "


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
while getopts ":i:o:hmpv" opt; do
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
      ;;
    p)
      echo "Implementing phase component of FDF to DICOM conversion."
      ;;
    v)
      echo "Setting verbose to 1."
      VERBOSE=1
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
    output_dir="$(dirname ${input_dir})/$(basename ${input_dir})/.dcm"
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

magphflag=`ls ${input_dir}/magnitude.img ${input_dir}/phase.img` 2> /dev/null
if [ $? -eq 0 ]; then  
    echo "Input directory has 'magnitude.img' and 'phase.img' "
    $0 -m -i ${input_dir}/magnitude.img/ -o ${output_dir}/magnitude.dcm
    if [ "$VERBOSE" -eq 1 ] && ! yesno "Magnitude complete. Continue converting phase?"; then
	echo "fdf2dcm completed without phase." && exit 0
    fi	
    $0 -p -i ${input_dir}/phase.img -o ${output_dir}/phase.dcm
    exit 0
fi

fdffiles=`ls ${input_dir}/*.fdf`
if [ $? -ne 0 ]; then  #-o "$fdffiles" == "" 
    echo "Error: Input directory has no FDF images"
    exit 1
fi
procfiles=`ls ${input_dir}/procpar`
if [ $? -ne 0 ]; then  #-o "$procfiles" == "" 
    echo "Error: Input directory has no procpar file"
    exit 1
fi

set -o errexit  # -e
set -o pipefail


## Crux of script
echo  "Calling agilent2dicom"
${FDF2DCMPATH}/agilent2dicom $@

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
    for iecho in $(seq 1 ${nechos}); do
	echoext=$(printf '%03d' $iecho)
	echo "Converting echo ${iecho} using dcmulti"
	${DCMULTI} "${output_dir}/0${echoext}.dcm" $(ls -1 ${output_dir}/tmp/*echo${echoext}.dcm)
    done
    rm -f ${output_dir}/MULTIECHO
else
    ${DCMULTI} ${output_dir}/0001.dcm $(ls -1 ${output_dir}/tmp/*.dcm)
fi

#if 0
#then

DCMODIFY=`which dcmodify`  #"$(dirname $0)/dmodify"
files=$(find ${output_dir} -name "*.dcm")
 # ${DCMODIFY} -m "(0020,0060)=" $files  # Laterality  # fixed in agilent2dicom
 ${DCMODIFY} -i "(5200,9229)[0].(0018,9125)[0].(0018,1312)=ROW" $files # In-plane phase encoding direction

firsttmpdcm=$(ls -1 ${output_dir}/tmp/*.dcm| head -1)
declare -a FrameAnatomySequence=(`dcdump ${firsttmpdcm} 2>&1 >/dev/null | grep '(0x0008,0x010' | awk -F'>' '/</ {print $4}'| tr -d '<'`)
echo "Frame Anatomy Seq: size: " ${#FrameAnatomySequence[*]}
if [[ ${#FrameAnatomySequence[*]} -ne 8 ]];then
    echo "DCM modification error. Not enough Frame Anatomy Sequence parameters."
    echo " Ignoring Anatomy modifications."
else
echo "FrameAnt 7: " ${FrameAnatomySequence[7]}
${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2218)[0].(0008,0100)=${FrameAnatomySequence[0]}" $files
${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2218)[0].(0008,0104)=${FrameAnatomySequence[1]}" $files   #CodeMeaning=
${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2218)[0].(0008,0102)=SRT" $files
# ${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2220)[0].(0008,0100)=${FrameAnatomySequence[2]}" $files   #CodeValue=
# ${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2220)[0].(0008,0104)=${FrameAnatomySequence[3]}" $files
# ${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2220)[0].(0008,0102)=SRT" $files
${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2228)[0].(0008,0100)=${FrameAnatomySequence[4]}" $files
${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2228)[0].(0008,0104)=${FrameAnatomySequence[5]}" $files
${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2228)[0].(0008,0102)=SRT" $files
# ${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2230)[0].(0008,0100)=${FrameAnatomySequence[6]}" $files
# ${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2230)[0].(0008,0104)=${FrameAnatomySequence[7]}" $files
# ${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2230)[0].(0008,0102)=SRT" $files
${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0020,9072)=R" $files
#${DCMODIFY} -i "(0018,9125)[0].(0018,1312)=ROW" $files
fi # array check
echo "DCModify done"
# fi #debugging



## Cleaning up
if [ -d ${output_dir}/tmp ]; then

    if [ "$VERBOSE" -eq 1 ] && ! (yesno "Remove existing tmp output directory, y or n (default y)?"); then
	echo "fdf2dcm completed. Temporary dicoms still remain." && exit 0
    fi

    echo "Removing existing tmp output directory"
    rm -rf ${output_dir}/tmp    

fi

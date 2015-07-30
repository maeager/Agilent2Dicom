#!/usr/bin/env  bash 
## Enhanced DICOM converter
#   Front-end to dcmulti
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
###############################################

function error_exit(){                                                               
     echo "${PROGNAME}: ${1:-Unknown error}" 1>&2                                     
     # if [ -x `which mutt` && "$USER" == "vnmr1" ]; then                             
     #     logfiles=$(find ${FID2DCMPATH} ${PWD} -name *.log -print0)                 
     #   EMAIL="$USER@$HOST"                                                          
     #   echo "Error occured `date`" | mutt -s "${PROGNAME}: ${1:-Unknown error}" -a $logfiles  michael.eager@monash.edu                                                   
      # fi                                                                         

     exit 1                                                                       
}

## Set config variables
FID2DCMPATH="$(dirname "$0")"
source "${FID2DCMPATH}/agilent2dicom_globalvars.py"
PROGNAME=$(basename "$0")

if test -x dciodvfy;then
    DCM3TOOLS=$(dirname $(which dciodvfy))
else
# KERNEL_RELEASE=$(uname -r | awk -F'.' '{printf("%d.%d.%d\n", $1,$2,$3)}')
# DCM3TOOLS="${FID2DCMPATH}/../dicom3tools_1.00.snapshot.20140306142442/bin/1.${KERNEL_RELEASE}.x8664/"
    DCM3TOOLS=$(/bin/ls -d "${FID2DCMPATH}"/../dicom3tools_*/bin/*)
#DCM3TOOLS="${FID2DCMPATH}/../dicom3tools_1.00.snapshot.20140306142442/bin/"
#DCM3TOOLS=$(echo "${DCM3TOOLS}"$(ls "${DCM3TOOLS}")"/")
fi
if [ ! -d "${DCM3TOOLS}" ]; then
    error_exit "${DCM3TOOLS} path not found"
elif [ ! -x "${DCM3TOOLS}/dcmulti" ]; then
    error_exit  "Unable to find Dicom3Tool's executable dcmulti"
fi   
export PATH=${PATH}:${DCM3TOOLS}
RM='/bin/rm -f'
RMDIR='/bin/rm -rf'
## Set dcmulti's arguments
DCMULTI="dcmulti -v -makestack -sortby AcquisitionNumber -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
#DCMULTI_DTI="dcmulti -v -makestack -sortby DiffusionBValue -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
#-makestack -sortby ImagePositionPatient  -sortby AcquisitionNumber

## Check dcmtk applications on MASSIVE or Agilent console
if test ${MASSIVEUSERNAME+defined}; then
    test -x dcmodify || module load dcmtk 
else
    DCMTK="/home/vnmr1/src/dcmtk-3.6.0/bin"
    export PATH=${PATH}:${DCMTK}
fi


declare -i verbosity=0
declare stddcmdir=""
declare output_dir=""
declare -i echoes=1
declare -i bdirs=0

### Debugging variables and config
set -o nounset  # shortform: -u
set -o errexit  # -e
# set -o pipefail

# touch $(dirname $0)/error.log
# exec 2>> $(dirname $0)/error.log
## set -x  # show debugging output

E_BADARGS=65


function yesno(){
    read -r -p "$@" response
    response=$(echo "$response" | awk '{print tolower($0)}')
# response=${response,,} # tolower Bash 4.0
    if [[ "$response" =~ ^(yes|y| ) ]]; then
	return 0
    fi
    return 1
}


# Print usage information and exit
print_usage(){
    echo -e "\n" \
	"usage: ./enh_mr.sh -i inputdir -o <outputdir> [-v] \n

  Convert standard MR files to enhanced MR format.\n" \
	"\n" \
        "-i <inputdir>  Standard DICOM directory. \n" \
	"-o <outputdir> Directory for output enhanced Dicom files. \n" \
        "-e <echos>     Number of echoes. \n" \
        "-d <bdirs>     Diffusion directions. \n" \
	"-x             Debug mode. \n" \
	"-v             Verbose output. \n" \
	"-h             this help\n" \
	"\n" 
    # && exit 1
#	"-s             Sequence type (one of MULTIECHO,DIFFUSION,ASL. \n" \
}


## Check for number of args
if [ $# -eq 0 ]; then
	echo "enh_mr.sh must have one argument: -i, --input [directory of standard DICOM images]"
	print_usage
	exit $E_BADARGS
fi


## Parse arguments
while getopts ":i:o:e:d:hxv" opt; do
    case $opt in
	i)
	    echo "Input dir of standard DICOMs:  $OPTARG" >&2
	    stddcmdir="$OPTARG"
	    ;;
	o)
	    echo "Output dir: $OPTARG" >&2
	    output_dir="$OPTARG"
	    ;;
	e)
	    echo "Multi echo: $OPTARG" >&2
	    echoes="$OPTARG"
	    ;;
	d)
	    echo "Diffusion: $OPTARG" >&2
	    bdirs="$OPTARG"
	    ;;       	
	v)
	    ((++verbosity))
	    echo "Setting verbose to $verbosity."
	    ((verbosity==1)) && python_args="$python_args -v"
	    ;;
	x)
	    set -x  ## print all commands
	    log_file="$(dirname $0)/error.log"
	    exec &> >(tee -a "$log_file")
	    # exec 2> $(dirname $0)/error.log
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
if [ ! -d "$stddcmdir" ]; then
    echo "enh_mr.sh must have a valid input directory of standard DICOM images."
    exit $E_BADARGS
fi
## Set output_dir if not in args, default is parent or stddcmdir
if [ -z "$output_dir" ]
then #test for empty string
    echo "Moving dicom images"
    mkdir "${stddcmdir}/tmp"
    mv "${stddcmdir}"/*.dcm "${stddcmdir}"/tmp/
    output_dir = "${stddcmdir}"
    stddcmdir="${output_dir}"/tmp
    echo "Output dir set to: ${output_dir}"
fi

## Check existing output directory for enhanced dicoms
dcmfiles="${output_dir}"/0*.dcm
for dcmfile in $dcmfiles; do
    if [ -f "${dcmfile}" ];then
	echo "Deleting " "${dcmfile}"
	${RM} "${dcmfile}"
    fi
done

    shopt -s nullglob  
    found=0
    for i in "${stddcmdir}"/*.dcm; do
	if [ -e "$i" ];then 
	    (( ++found ))
	else
	    error_exit "$LINENO: dicom file does not exist $i"
	fi
    done
    shopt -u nullglob
    if [ $found -eq 0 ]; then  
	error_exit "$LINENO: Input directory has no Dicom images"
    else
	echo $found, " Dicom files were found"
    fi

##############################

    # output_root="$(dirname "${output_dir}")"
    # output_base=$(basename "${output_dir}" .dcm)
    # dirs=$(find "${output_root}" -maxdepth 1 -type d  -name  "$output_base*.dcm")
    # echo "fid2dicom.py completed successfully. Dicom paths generated were: "
    # echo "$dirs"

    # dcmfiles=$(ls ${output_dir}/*.dcm)  ## Bad code - use glob
    #if[ $? -ne 0 ]
    # test -e "${output_dir}"/0001.dcm && error_exit "$LINENO: Output directory of fid2dicom has no DICOM images."
    


echo "Convert dicom images to single enhanced MR dicom format image"
if [ -f "${stddcmdir}/MULTIECHO" -o echoes > 1 ]
then
    if [ -f "${stddcmdir}/MULTIECHO" ]
    then 
	echo "Contents of MULTIECHO"; cat "${stddcmdir}/MULTIECHO"; printf '\n'
	nechos=$(cat "${stddcmdir}/MULTIECHO")
	nechos=$(printf "%1.0f" "$nechos")
    else 
	nechos=$echoes
    fi
    echo "Multi echo sequence, $nechos echos"
    for ((iecho=1;iecho<=nechos;++iecho)); do
     	echoext=$(printf '%03d' $iecho)
     	echo "Converting echo ${iecho} using dcmulti"
     	${DCMULTI} "${output_dir}"/0"${echoext}".dcm $(ls -1 "${stddcmdir}"/*echo"${echoext}".dcm | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
	[ $? -ne 0 ] && error_exit "$LINENO: dcmulti MULTI ECHO failed"
    done


#    ${RM} "${output_dir}/MULTIECHO"
    echo "Multi echo sequence completed."
    
elif  [ -f "${stddcmdir}/DIFFUSION"  ]; then
    #if  [ -f "${stddcmdir}/DIFFUSION" ]; then
#	echo "Contents of DIFFUSION"; cat "${stddcmdir}/DIFFUSION"; printf '\n'

    # nbdirs=$(cat ${output_dir}/DIFFUSION)
    # ((++nbdirs)) # increment by one for B0
    # Get b directions from maximum image number
	nbdirs=$(ls -1 "${stddcmdir}/slice*" | sed 's/.*image0\(.*\)echo.*/\1/' | tail -1)
 #   else
#	nbdirs=$bdirs
#    fi
    echo "Diffusion sequence, $nbdirs B-directions"
    for ((ibdir=1;ibdir<=nbdirs;ibdir++)); do
     	bdirext=$(printf '%03d' $ibdir)
     	echo "Converting bdir ${ibdir} using dcmulti"
	#for dcmdir in $dirs; do
	## Input files are sorted by image number and slice number. 
     	    ${DCMULTI} "${output_dir}/0${bdirext}.dcm" $(ls -1 "${stddcmdir}"/*image"${bdirext}"*.dcm | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
	#done
	    [ $? -ne 0 ] && error_exit "$LINENO: dcmulti DIFFUSION failed"
    done
    echo "Diffusion files compacted."

elif  [ -f "${output_dir}/ASL" ]; then

    echo "Contents of ASL"; cat "${output_dir}/ASL"; printf '\n'

    # nbdirs=$(cat ${output_dir}/ASL)
    # ((++nbdirs)) # increment by one for B0
    # asltags=$(ls -1 "${output_dir}/tmp/slice*" | sed 's/.*image0\(.*\)echo.*/\1/' | tail -1)
    echo "ASL sequence"
    for ((iasl=1;iasl<=2;iasl++)); do
     	aslext=$(printf '%03d' $iasl)
     	echo "Converting ASL tag ${iasl} using dcmulti"
	
	## Input files are sorted by image number and slice number. 
     	${DCMULTI} "${output_dir}/0${aslext}.dcm" $(ls -1 "${stddcmdir}"/*echo"${aslext}".dcm | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
	
	[ $? -ne 0 ] && error_exit "$LINENO: dcmulti ASL failed"
    done
    echo "ASL files converted."


else

    ## Dcmulti config is dependent on order of files.  The 2D standard
    ## dicoms are sorted by echo time, image number then slice
    ## number. The second argument reorders the list of 2D dicom files
    ## based on echo time, then image number, then slice number.
    ## Only one output file is required, 0001.dcm. 
    
    ${DCMULTI} "${output_dir}/0001.dcm" $(ls -1 "${stddcmdir}"/*.dcm  | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')
    [ $? -ne 0 ] && error_exit "$LINENO: dcmulti failed"
    
fi
echo "DCMULTI complete. Fixing inconsistencies."

## Corrections to dcmulti conversion

#    for dcmdir in $dirs; do
	"${FID2DCMPATH}/fix-dicoms.sh" "${output_dir}"
	echo "Fixing dicom dir ${output_dir}"

    ## Additional corrections to diffusion files
	if [ -f "${output_dir}/DIFFUSION" -o bdirs > 1 ];then
	    "${FID2DCMPATH}"/fix-diffusion.sh "${output_dir}"
	    echo "Fixed diffusion module parameters."
	fi
    ## Additional corrections to ASL files
	if [ -f "${output_dir}/ASL" ];then
	    "${FID2DCMPATH}"/fix_asl.sh "${output_dir}"
	    echo "Fixed ASL module parameters."
	fi
#    done


if (( verbosity > 0 )); then
    echo "Verifying dicom compliance using dciodvfy."
    if [ -f "${output_dir}/0001.dcm" ]; then
	set +e
	## Send dciodvfy stderr and stdout to log file
	dciodvfy "${output_dir}/0001.dcm" &> "$(dirname "${output_dir}")/$(basename "${output_dir}" .dcm).log"
	set -e  
    else
	error_exit "$LINENO: could not find ${output_dir}/0001.dcm for verification"
    fi
fi


exit 0

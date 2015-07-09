#!/usr/bin/env  bash 
## FID to DICOM converter
#   Front-end to fid2dicom and dcmulti
#
# - Michael Eager (michael.eager@monash.edu.au)
# - Monash Biomedical Imaging 
#
#  "$Id: fid2dcm.sh,v 8c9c3b886d01 2015/03/10 05:23:17 michael $"
#  Version 0.1: FID2DCM based on FDF2DCM with fid2dicom core
#  Version 0.5: Major update to input args
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



## Set config variables
FID2DCMPATH="$(dirname "$0")"
source "${FID2DCMPATH}/agilent2dicom_globalvars.py"
## variable collected global FID2DCMVERSION=0.5
PROGNAME=$(basename "$0")
FID2DICOM=fid2dicom.py
KERNEL_RELEASE=$(uname -r | awk -F'.' '{printf("%d.%d.%d\n", $1,$2,$3)}')
DCM3TOOLS="${FID2DCMPATH}/../dicom3tools_1.00.snapshot.20140306142442/bin/1.${KERNEL_RELEASE}.x8664/"
DCM3TOOLS=$(/bin/ls -d "${FID2DCMPATH}"/../dicom3tools_*/bin/*)
#DCM3TOOLS="${FID2DCMPATH}/../dicom3tools_1.00.snapshot.20140306142442/bin/"
#DCM3TOOLS=$(echo "${DCM3TOOLS}"$(ls "${DCM3TOOLS}")"/")
RM='/bin/rm -f'
RMDIR='/bin/rm -rf'

## Check dcmtk applications on MASSIVE or Agilent console
if test ${MASSIVEUSERNAME+defined}; then
    test -x dcmodify || module load dcmtk 
else
    DCMTK="/home/vnmr1/src/dcmtk-3.6.0/bin"
    export PATH=${PATH}:${DCMTK}
fi
if [ ! -d "${DCM3TOOLS}" ]; then
    echo "${DCM3TOOLS} path not found"
    exit 1
elif [ ! -x "${DCM3TOOLS}/dcmulti" ]; then
    echo "Unable to find Dicom3Tool's executable dcmulti"
    exit 1
fi 
export PATH=${PATH}:${DCM3TOOLS}
declare -i verbosity=0
declare -i do_modify=1
declare -i do_filter=0
declare input_dir=""
declare output_dir=""
declare kspace_dir=""
declare python_args=""

### Debugging variables and config
set -o nounset  # shortform: -u
set -o errexit  # -e
# set -o pipefail

# touch $(dirname $0)/error.log
# exec 2>> $(dirname $0)/error.log
## set -x  # show debugging output

E_BADARGS=65

#source ${FID2DCMPATH}/yesno.sh
function yesno(){
    read -r -p "$@" response
    response=$(echo "$response" | awk '{print tolower($0)}')
# response=${response,,} # tolower Bash 4.0
    if [[ "$response" =~ ^(yes|y| ) ]]; then
	return 0
    fi
    return 1
}

function error_exit(){
    echo "${PROGNAME}: ${1:-Unknown error}" 1>&2
    # if [ -x `which mutt` && "$USER" == "vnmr1" ]; then
    #     logfiles=$(find ${FID2DCMPATH} ${PWD} -name *.log -print0)
    # 	EMAIL="$USER@$HOST" 
    # 	echo "Error occured `date`" | mutt -s "${PROGNAME}: ${1:-Unknown error}" -a $logfiles  michael.eager@monash.edu 
    # fi
    exit 1
}


# Print usage information and exit
print_usage(){
    echo -e "\n" \
	"usage: ./fid2dcm.sh -i inputdir [-o outputdir] [-v] [-m] [-p] [-r]
[-k] [-N] [[-g 1.0 -j 0 -e wrap] [-l 1.0] [-n 5] [ -w 5 -z 0.001][-y
1.0 -j 0 -e wrap]] [[-D] [-G 1.0 ] [-L 1.0] [-Y 1.0 ]]\n\n

  To export recon components use magnitude (-m), phase (-p), real and
  imag (-r) or k-space (-k). Filtering is available for Gaussian filter
  (-g sigma), Laplace Gaussian filter (-l sigma), median filter (-n
  window_size), Wiener filter (-w window_size) or Epanechnikov filter
  (-y <bandwidth>). K-space Fourier filtering is avalable for Gaussian
  (-G <sigma>), Laplace of Gaussian (-L) and Epanechnikov (-Y <bandwidth>) filters.\n" \
	"\n" \
	"-i <inputdir>  FID source directory\n" \
	"-o <outputdir> Optional destination DICOM directory. Default is input_dir/.dcm. \n" \
	"-d             Disable dcmodify fixes to DICOMs.\n" \
	"-m,-p,         Save magnitude and phase components.  These flags are passed to fid2dicom 
                        and should only be used from within fid2dcm or with knowledge of input fid data. \n" \
	"-r             Save real and imaginary components of filtered image. \n" \
	"-k             Save Kspace data. \n" \
	"-N             Save filtered outputs to NIFTI.\n" \
	"-g <sigma>     Gaussian filter smoothing of reconstructed RE and
	                IM components. Sigma variable argument, default 1/srqt(2). \n" \
	"-G <sigma>     Fourier Gaussian filter smoothing of k-space RE and
	                IM components. Sigma variable argument, default size/1/srqt(2). \n" \
	"-j <order>     Gaussian filter order variable argument, default 0. \n" \
	"-e {'wrap','reflect','nearest','mirror'}  Gaussian filter mode variable 
                        argument, default=nearest. \n" \
	"-l <simga>     Gaussian Laplace filter smoothing of reconstructed RE and 
                        IM components. Sigma variable argument, default 1/srqt(2).\n" \
	"-L <sigma>     Fourier Laplace-of-Gaussian filter smoothing of k-space RE and
	                IM components. Sigma variable argument, default size/1/srqt(2). \n" \
	"-n <wsize>     Median filter smoothing of reconstructed RE and IM components. 
                        Size variable argument, default 5.\n" \
	"-s <wsize>     Local Std deviation filter of reconstructed RE and IM components. 
                        Size variable argument, default 5.\n" \
	"-w <wsize>     Wiener filter smoothing of reconstructed RE and IM components. 
                        Size variable argument, default 5.\n" \
	"-z <noise>     Wiener filter noise variable. 
                        Default 0 (or none) variance calculated in local region and 
                        can be quite computationally expensive.\n" \
	"-y <bwidth>    Epanechnikov filter smoothing of reconstructed RE and
	                IM components. Bandwidth variable argument, default sqrt(Dim+4))/srqt(2). \n" \
	"-Y <bwidth>    Fourier Epanechnikov filter. Smoothing of k-space RE and
	                IM components. Sigma variable argument, default size/sqrt(Dim+4)/srqt(2). \n" \
	"-D             Double resolution in k-space with zero-padding and recon. 
                        !!!Warning!!! this increases size on disk by a factor of 8.\n" \
	"-C             Disable k-space shift centering to maxima \n" \
	"-x             Debug mode. \n" \
	"-v             Verbose output. \n" \
	"-h             this help\n" \
	"\n" 
    # && exit 1
#	"-s             Sequence type (one of MULTIECHO,DIFFUSION,ASL. \n" \
}


## Check for number of args
if [ $# -eq 0 ]; then
	echo "fid2dcm.sh must have one argument: -i, --input [directory of FID images]"
	print_usage
	exit $E_BADARGS
fi


## Parse arguments
while getopts ":i:o:g:l:j:e:n:s:w:y:z:G:L:Y:DhmprkdNCxv" opt; do
    case $opt in
	i)
	    echo "Input dir:  $OPTARG" >&2
	    input_dir="$OPTARG"
	    ;;
	o)
	    echo "Output dir: $OPTARG" >&2
	    output_dir="$OPTARG"
	    ;;
	k)
	    echo "K-space dir:  $OPTARG" >&2
	    kspace_dir="$OPTARG"
	    python_args="$python_args --kspace $kspace_dir"
	    ;;
	g)
	    echo "Gaussian filter sigma: $OPTARG" >&2
	    gaussian_sigma="$OPTARG"
	    python_args="$python_args --gaussian_filter --sigma $gaussian_sigma --gaussian_order 0 --gaussian_mode nearest"
	    do_filter=1
	    ;;
	G)
	    echo "Fourier domain Gaussian filter sigma: $OPTARG" >&2
	    gaussian_sigma="$OPTARG"
	    python_args="$python_args --FT_gaussian_filter --sigma $gaussian_sigma --gaussian_order 0 --gaussian_mode nearest"
	    do_filter=1
	    ;;
	j)
	    [ ${do_filter} != 1 ] && (echo "Must have -g before -j"; print_usage; exit 1)
	    echo "Gaussian filter order: $OPTARG" >&2
	    gaussian_order="$OPTARG"
	    #python_args=$(echo $python_args | sed 's/--gaussian_order\s[0-3]\s/--gaussian_order '$gaussian_order' /')
	    python_args="${python_args//--gaussian_order\s[0-3]\s/--gaussian_order $gaussian_order /}"
	    ;;
	e)
	    [ ${do_filter} != 1 ] && (echo "Must have -g before -e"; print_usage; exit 1)
	    echo "Gaussian filter mode: $OPTARG" >&2
	    gaussian_mode="$OPTARG"
	    python_args=$(echo $python_args | sed 's/--gaussian_mode\s\w+\s/--gaussian_mode $gaussian_mode/')
	    ;;
	l)
	    echo "Gaussian Laplace filter sigma: $OPTARG" >&2
	    gaussian_sigma="$OPTARG"
	    python_args="$python_args --gaussian_laplace --sigma $gaussian_sigma"
	    do_filter=2
	    ;;
	L)
	    echo "Fourier Domain Gaussian Laplace filter sigma: $OPTARG" >&2
	    gaussian_sigma="$OPTARG"
	    python_args="$python_args --FT_gaussian_laplace --sigma $gaussian_sigma"
	    do_filter=2
	    ;;
	n)
	    echo "Median filter size: $OPTARG" >&2
	    median_window_size="$OPTARG"
	    python_args="$python_args --median_filter --window_size $median_window_size"
	    do_filter=3
	    ;;
	s)
	    echo "Standard deviation filter [mode/size]: $OPTARG" >&2
	    stdev_mode="${OPTARG:0:1}"
	    stdev_window_size="${OPTARG:2}"
	    [ "$stdev_mode" == "0" ] && python_args="$python_args --standard_deviation_filter"
	    [ "$stdev_mode" == "1" ] && python_args="$python_args --standard_deviation_magn"
	    [ "$stdev_mode" == "2" ] && python_args="$python_args --standard_deviation_phase"
	    python_args="$python_args --window_size $stdev_window_size"
	    do_filter=6
	    ;;	
	w)
	    echo "Wiener filter size: $OPTARG" >&2
	    wiener_window_size="$OPTARG"
	    python_args="$python_args --wiener_filter --window_size $wiener_window_size"
	    do_filter=4
	    ;;
	z)
	    [ ${do_filter} != 4 ] && (echo "Must have -w before -z"; print_usage; exit 1)
	    echo "Wiener noise: $OPTARG" >&2
	    wiener_noise="$OPTARG"
	    python_args="$python_args --wiener_noise $wiener_noise"
	    ;;
	y)
	    echo "Epanechnikov  filter size: $OPTARG" >&2
	    epan_bandwidth="$OPTARG"
	    python_args="$python_args --epanechnikov_filter --epanechnikov_bandwidth $epan_bandwidth"
	    do_filter=5
	    ;;
	Y)
	    echo "Fourier Epanechnikov  filter size: $OPTARG" >&2
	    epan_bandwidth="$OPTARG"
	    python_args="$python_args --FT_epanechnikov_filter --epanechnikov_bandwidth $epan_bandwidth"
	    do_filter=5
	    ;;
	D)
	    echo "Implementing super-resolution of FID to double the resolution of DICOM conversion."
	    python_args="$python_args --double_resolution"
	    ;;
	C)
	    echo "Disable kspace shift function to centre data to maximum."
	    python_args="$python_args --no_centre_shift"
	    ;;
	h)
	    print_usage
	    # "${FID2DCMPATH}"/fid2dicom.py -h
	    exit 0
	    ;;
	m)
	    echo "Implementing magnitude component of FID to DICOM conversion."
	    python_args="$python_args --magnitude"
	    ;;
	r)
	    echo "Save real and imaginary components of FID conversion."
	    python_args="$python_args --realimag"
	    ;;
	p)
	    echo "Implementing phase component of FID to DICOM conversion."
	    python_args="$python_args --phase"
	    ;;
	N)
	    echo "Saving filtered outputs to NIFTI."
	    python_args="$python_args --nifti"
	    ;;
	# s)
	#     echo "Sequence type: $OPTARG" >&2
	#     sequence="$OPTARG"
	#     python_args="$python_args -s $sequence"
	#     ;;
	d)
	    do_modify=0
	    echo " Disable dcmodify correction."
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
if [ ! -d "$input_dir" ]; then
    echo "fid2dcm.sh must have a valid input directory of FID images."
    exit $E_BADARGS
fi
## Set output_dir if not in args, default is INPUT/.dcm
if [ -z "$output_dir" ]
then #test for empty string
    output_dir=$(dirname "${input_dir}")/$(basename "${input_dir}" .fid).dcm
    echo "Output dir set to: ${output_dir}"
fi
## Set kspace_dir if not in args, default is INPUT.dcm
if [ "$kspace_dir" != "" ]; then 
    kspace_dir=$(dirname "${output_dir}")/$(basename "${input_dir}" .fid)_kspace.dcm
    echo "K space output dir set to: ${kspace_dir}"
    [ ! -d "$kspace_dir" ] && mkdir -p "$kspace_dir"
fi

## Check existing output directories
JumpToDCmulti=0
output_root=$(dirname "${output_dir}")
output_base=$(basename "${output_dir}" .dcm)
dirs=$(find "$output_root" -maxdepth 1 -type d  -name  "$output_base*.dcm")
echo "Output dicom paths exist already: "
echo "${dirs}"
if [ -d "${output_dir}" ]; then
    if test -d "${output_dir}/tmp" && (( verbosity > 0 ))
    then
	if yesno "Remove existing output directory, y or n (default y)?"; then
	    echo "Removing existing output directories"
	    for dcmdir in ${dirs}; do
		${RMDIR} "${dcmdir}"
	    done
	else
	    JumpToDCmulti=1
	fi
    else
	echo "Removing existing output directories"
	for dcmdir in ${dirs}; do
	    ${RMDIR} "${dcmdir}"
	done
    fi	
fi

if (( JumpToDCmulti == 0 ))
then

    shopt -s nullglob  
    found=0
    for i in "${input_dir}"/fid; do
	if [ -e "$i" ];then 
	    (( ++found ))
	else
	    error_exit "$LINENO: fid file does not exist $i"
	fi
    done
    shopt -u nullglob
    if [ $found -eq 0 ]; then  #-o "$fidfiles" == "" 
	error_exit "$LINENO: Input directory has no FID images"
    else
	echo $found, " FID files were found"
    fi

    if [ ! -f "${input_dir}/procpar" ]; then
	error_exit "$LINENO: Input directory has no procpar file"
    fi

# set -o errexit  # -e
# set -o pipefail


## Crux of script - conversion of FID images to standard DICOM images
    echo  "Calling fid2dicom"
    echo " Arguments: ${python_args} --inputdir ${input_dir} --outputdir ${output_dir}"
    "${FID2DCMPATH}/${FID2DICOM}" ${python_args} --inputdir "${input_dir}" --outputdir "${output_dir}"

    [ $? -ne 0 ] && error_exit "$LINENO: fid2dicom failed"
    
    [ ! -d "${output_dir}" ] && error_exit "$LINENO: Output dir not created by fid2dicom."

    output_root="$(dirname "${output_dir}")"
    output_base=$(basename "${output_dir}" .dcm)
    dirs=$(find "${output_root}" -maxdepth 1 -type d  -name  "$output_base*.dcm")
    echo "fid2dicom.py completed successfully. Dicom paths generated were: "
    echo "$dirs"

    # dcmfiles=$(ls ${output_dir}/*.dcm)  ## Bad code - use glob
    #if[ $? -ne 0 ]
    # test -e "${output_dir}"/0001.dcm && error_exit "$LINENO: Output directory of fid2dicom has no DICOM images."
    
    echo "Moving dicom images"
    for dcmdir in ${dirs}; do
	mkdir "${dcmdir}/tmp"
	mv "${dcmdir}"/*.dcm "${dcmdir}"/tmp/
    done
fi ## JumpToDCMulti

echo "Convert dicom images to single enhanced MR dicom format image"
if [ -f "${output_dir}/MULTIECHO" ]
then
    echo "Contents of MULTIECHO"; cat "${output_dir}/MULTIECHO"; printf '\n'
    nechos=$(cat "${output_dir}/MULTIECHO")
    nechos=$(printf "%1.0f" "$nechos")
    
    for dcmdir in $dirs; do
	"${FID2DCMPATH}"/enh_mr.sh -e $nechos -i "${dcmdir}"/tmp/ -o "${dcmdir}"
    done
    

# DCMULTI="dcmulti -v -makestack -sortby EchoTime -dimension StackID FrameContentSequence -dimension InStackPositionNumber FrameContentSequence -of "
#-makestack -sortby ImagePositionPatient  -sortby AcquisitionNumber
#  ${DCMULTI} ${output_dir}/0001.dcm $(ls -1 ${output_dir}/tmp/*.dcm  | sed 's/\(.*\)slice\([0-9]*\)image\([0-9]*\)echo\([0-9]*\).dcm/\4 \3 \2 \1/' | sort -n | awk '{printf("%sslice%simage%secho%s.dcm\n",$4,$3,$2,$1)}')


    echo "Multi echo sequence completed."
    
elif  [ -f "${output_dir}/DIFFUSION" ]; then

    echo "Contents of DIFFUSION"; cat "${output_dir}/DIFFUSION"; printf '\n'
    for dcmdir in $dirs; do
	"${FID2DCMPATH}"/enh_mr.sh -d 1 -i "${dcmdir}"/tmp/ -o "${dcmdir}"
    done

    echo "Diffusion files completed."

else
    ## ASL and standard dicoms
    for dcmdir in $dirs; do
	"${FID2DCMPATH}"/enh_mr.sh -i "${dcmdir}"/tmp/ -o "${dcmdir}"
    done
fi
[ -f "${output_dir}/MULTIECHO" ] && ${RM} "${output_dir}/MULTIECHO"
[ -f "${output_dir}/DIFFUSION" ] && ${RM} "${output_dir}/DIFFUSION"
[ -f "${output_dir}/ASL" ] && ${RM} "${output_dir}/ASL"



## Cleaning up temporary directories
echo "Cleaning up tmp directories."
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
exit 0

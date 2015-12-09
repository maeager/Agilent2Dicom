#!/usr/bin/env bash

# - Michael Eager (michael.eager@monash.edu)
# - Monash Biomedical Imaging 
# - (C) 2014 Michael Eager

VERBOSE=0

E_BADARGS=65

#Check whether we are on MASSIVE or the Agilent console (Redhat 6.5)
if test ${MASSIVE_USERNAME+defined}; then 
    FDF2DCMPATH=$(dirname $0)
    # export PATH=${PATH}:${DCMTK}
else
    DCMTKPATH="/home/vnmr1/src/dcmtk-3.6.0/bin"
    FDF2DCMPATH="/home/vnmr1/src/Agilent2Dicom.git"
    echo $DCMTKPATH

# load dcm tools
# DCMODIFY=./dcmtools/bin/DCMODIFY
    export PATH=${PATH}:${DCMTKPATH}
fi

echo $FDF2DCMPATH
#KERNEL_RELEASE=$(uname -r | awk -F'.' '{printf("%d.%d.%d\n", $1,$2,$3)}')
#DCM3TOOLS="${FDF2DCMPATH}/../dicom3tools_1.00.snapshot.20140306142442/bin/1.${KERNEL_RELEASE}.x8664/"
DCM3TOOLS=$(/bin/ls -d "${FDF2DCMPATH}"/../dicom3tools_*/bin/*)
export PATH=${PATH}:${DCM3TOOLS}

DCIODVFY=${DCM3TOOLS}/dciodvfy


# Print usage information and exit
print_usage(){
    echo -e "\n" \
	"usage: ./dcheck.sh -o outputdir\n" \
	"\n" \
	"-o <outputdir> DICOM directory.  \n" \
	"-v             verbose output. \n" \
	"-h             this help\n" \
	"\n" 
    # && exit 1
}


# Check for number of args
if [ $# -eq 0 ]; then
	echo "dcheck.sh must have one argument: -i, --input [directory of FDF images]"
	print_usage
	exit $E_BADARGS
fi



# Parge arguments
while getopts ":o:hv" opt; do
    case $opt in
	o)
	    echo "Output dir: $OPTARG" >&2
	    output_dir="$OPTARG"
	    ;;
	h)
	    print_usage
	    exit 0
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
if [ ! -d "$output_dir" ]; then
    echo "dcheck.sh must have a valid DICOM directory."
    exit $E_BADARGS
fi

if [ ! -f "$output_dir"/0001.dcm ]; then
    echo "dcheck.sh must have a valid DICOM file in the directory."
    exit $E_BADARGS
fi

if [ ! -x ${DCIODVFY} ]; then
    echo "dciodvfy not found"
    exit 1
fi

${DCIODVFY} "$output_dir"/0001.dcm

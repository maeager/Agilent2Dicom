#!/usr/bin/env bash
#
# - Michael Eager (michael.eager@monash.edu)
# - Monash Biomedical Imaging 
# - (C) 2014 Michael Eager
#
############################################



# Check DCMTK on MASSIVE or Agilent console
if test "${MASSIVEUSERNAME+defined}"; then
    if [ ! -x $(which dcmodify) ];then
	module load dcmtk
    fi
else
    DCMTK="/home/vnmr1/src/dcmtk-3.6.0/bin"
    export PATH=${PATH}:${DCMTK}

fi

if [ ! -x $(which dcmodify) ];then
    echo "fix-dicoms.sh: dcmodify not found."; 
    exit 1
fi
if [ ! -x $(which dcdump) ];then
    echo "fix-dicoms.sh: dcdump not found."; 
    exit 1
fi
if [ ! -x $(which dciodvfy) ];then
    echo "fix-dicoms.sh: dciodvfy not found."; 
    exit 1
fi


output_dir=$1
MODIFY=1
##COMMON FIXES to enhanced DICOMs
DCMODIFY="dcmodify --no-backup  " # --ignore-errors" 

# Find dcmulti converted DICOMs - do not descending into tmp directory
files=$(find ${output_dir} -maxdepth 1 -type f -name "*.dcm")



 # ${DCMODIFY} -m "(0020,0060)=" $files  # Laterality  # fixed in agilent2dicom
# In-plane phase encoding direction
${DCMODIFY} -i "(5200,9229)[0].(0018,9125)[0].(0018,1312)=ROW" $files 


# Transmit Coil Type
  #   > (0x0018,0x9051) CS Transmit Coil Type      VR=<CS>   VL=<0x0008>  <UNKNOWN >
transcoiltype=$(dcdump ${output_dir}/tmp/slice001image001echo001.dcm 2>&1 >/dev/null | grep 'Transmit Coil Type' | awk '{print $9}' | tr -d '<>' ) 
echo "Fixing Receive Coil Type :" $transcoiltype
${DCMODIFY} -m "(5200,9229)[0].(0018,9049)[0].(0018,9051)=$transcoiltype" $files


# Tranmitter Frequency
  # > (0x0018,0x9098) FD Transmitter Frequency   VR=<FD>   VL=<0x0008>  {0}
transfrq=$(dcdump ${output_dir}/tmp/slice001image001echo001.dcm 2>&1 >/dev/null | grep 'Transmitter Frequency' | sed 's/^.*{\(.*\)} *$/\1/' )
echo "Fixing Tranmitter Frequency :" $transfrq 
${DCMODIFY} -m "(5200,9229)[0].(0018,9006)[0].(0018,9098)=$transfrq" $files


 # > (0x0018,0x9042) SQ MR Receive Coil Sequence        VR=<SQ>   VL=<0xffffffff>
 #    > (0x0018,0x9043) CS Receive Coil Type       VR=<CS>   VL=<0x0008>  <UNKNOWN >
reccoiltype=$(dcdump ${output_dir}/tmp/slice001image001echo001.dcm 2>&1 >/dev/null | grep 'Receive Coil Type' | awk '{print $9}' | tr -d '<>' ) 
echo "Fixing Receive Coil Type :" $reccoiltype
${DCMODIFY} -m "(5200,9229)[0].(0018,9042)[0].(0018,9043)=$reccoiltype" $files



if [[ $output_dir = *cine* ]]  ## pattern match cine in output directory string
then
    echo "Disabling Frame anatomy modiufication in CINE";
    MODIFY=0
	${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2218)[0].(0008,0102)=SRT" $files
	${DCMODIFY} -i "(5200,9229)[0].(0020,9071)[0].(0008,2228)[0].(0008,0102)=SRT" $files

fi

firsttmpdcm=$(find ${output_dir}/tmp/ -name "*.dcm"  | head -1)

# multiple spin echo (0018,9011) - not in diffusion or asl
multspinecho=$(dcdump $firsttmpdcm 2>&1 | grep 'Multiple Spin Echo' | awk '{print $8}' | tr -d '<>')
${DCMODIFY} -i "(0018,9011)=$multspinecho" $files

if (( MODIFY == 1 )); then
#"$(dirname $0)/dmodify"

    
    dcdump "${firsttmpdcm}" 2>&1 >/dev/null | grep '(0x0008,0x010' | awk -F'>' '/</ {print $4}'| tr -d '<' > ${output_dir}/anatomy.tmp
    if [ -f  ${output_dir}/anatomy.tmp ]; then
	declare -a FrameAnatomySequence
	let i=0;while IFS=$'\r\n' read -r line_data; do 
	    FrameAnatomySequence[i]="${line_data}"; ((++i)); 
	done < ${output_dir}/anatomy.tmp
	rm ${output_dir}/anatomy.tmp
    else
	echo "Cannot find anatomy.tmp in output dir"
    fi

    echo "Frame Anatomy Seq: " ${firsttmpdcm} " size: " ${#FrameAnatomySequence[*]}
    if [ ${#FrameAnatomySequence[*]} -ne 8 ];then
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
fi #debugging modify


echo "Removing Per-frame Anatomy sequences"
index=0
total_anatseq=$(dciodvfy ${output_dir}/0001.dcm 2>&1 >/dev/null | grep -e '^Error - Functional Group Sequence already used in Shared Functional Groups Sequence - (0x0020,0x9071) Frame Anatomy Sequence - in Per-frame Functional Groups Sequence' | wc -l)
echo "Total Frame Anatomy Errors ", "$total_anatseq" 
current_anatseq=$total_anatseq
#while (( current_anatseq > 0 )); do
    for ((i=0;index<total_anatseq;++index)); do
	echo "# $index of $total_anatseq"
	${DCMODIFY} -ea "(5200,9230)[$index].(0020,9071)" "$files"
#	((++index))
    done
#    current_anatseq=$(dciodvfy ${output_dir}/0001.dcm 2>&1 >/dev/null | grep -e '^Error - Functional Group Sequence already used in Shared Functional Groups Sequence - (0x0020,0x9071) Frame Anatomy Sequence - in Per-frame Functional Groups Sequence' | wc -l)
#done


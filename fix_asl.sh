#!/usr/bin/env bash
## fix_asl - Correct DCMULTI headers for ASL sequences

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



# $ dcdump ../output_data/ASL_se_06.dcm/tmp/slice001image001echo001.dcm 2>&1| grep -A30 'MR Arterial' | sed 's/.*VL=<[^>]*>\s\(.*\)/\1/'                                
# (0x0018,0x9251) SQ MR Arterial Spin Labeling Sequence    VR=<SQ>   VL=<0xffffffff>
#   ----:
#     > (0x0018,0x9250) CS Arterial Spin Labeling Contrast         VR=<CS>   VL=<0x000a>  <CONTINUOUS>
#     > (0x0018,0x9252) LO ASL Technique Description       VR=<LO>   VL=<0x0004>  <FAIR>
#     > (0x0018,0x9257) CS ASL Context     VR=<CS>   VL=<0x0006>  <LABEL >
#     > (0x0018,0x9259) CS ASL Crusher Flag        VR=<CS>   VL=<0x0002>  <NO>
#     > (0x0018,0x925b) LO ASL Crusher Description         VR=<LO>   VL=<0x0014>  <crusher description >
#     > (0x0018,0x925c) CS ASL Bolus Cut-off Flag          VR=<CS>   VL=<0x0002>  <NO>
#     > (0x0018,0x925d) SQ ASL Bolus Cut-off Timing Sequence       VR=<SQ>   VL=<0xffffffff>
#   ----:
#     > (0x0018,0x925e) LO ASL Bolus Cut-off Technique     VR=<LO>   VL=<0x0000>  <>
#     > (0x0018,0x925f) UL ASL Bolus Cut-off Delay Time    VR=<UL>   VL=<0x0004>  [0x00000000]

#     > (0x0018,0x9260) SQ ASL Slab Sequence       VR=<SQ>   VL=<0xffffffff>
#   ----:
#     > (0x0008,0x2218) SQ Anatomic Region Sequence        VR=<SQ>   VL=<0xffffffff>
#   ----:
#     > (0x0008,0x0100) SH Code Value      VR=<SH>   VL=<0x0008>  <T-D1100 >
#     > (0x0008,0x0104) LO Code Meaning    VR=<LO>   VL=<0x0004>  <Head>
#     > (0x0008,0x2220) SQ Anatomic Region Modifier Sequence       VR=<SQ>   VL=<0xffffffff>
#   ----:
#     > (0x0008,0x0100) SH Code Value      VR=<SH>   VL=<0x0006>  <G-A138>
#     > (0x0008,0x0104) LO Code Meaning    VR=<LO>   VL=<0x0008>  <Coronal >

#     > (0x0018,0x9253) US ASL Slab Number         VR=<US>   VL=<0x0002>  [0x0001]
#     > (0x0018,0x9254) FD ASL Slab Thickness      VR=<FD>   VL=<0x0008>  {12}
#     > (0x0018,0x9255) FD ASL Slab Orientation    VR=<FD>   VL=<0x0018>  {1.38,180.661,87.387}
#     > (0x0018,0x9256) FD ASL Mid Slab Position   VR=<FD>   VL=<0x0018>  {0,0,0}
#     > (0x0018,0x9258) UL ASL Pulse Train Duration        VR=<UL>   VL=<0x0004>  [0x00000000]
                      


firsttmpdcm=$(ls -1 "${output_dir}"/tmp/*.dcm| head -1)
dcdump "${firsttmpdcm}" 2>&1| grep -A30 'MR Arterial' | sed -e 's/.*VL=<[^>]*>\s*\(.*\)/\1/' -e '/---/d' -e '/^$/d' | tr -d '<>{}[]' > "${output_dir}"/arterial.tmp

if [ -f  ${output_dir}/arterial.tmp ]; then
    declare -a ArterialSpinLabelling
    let i=0;while IFS=$'\r\n' read -r line_data; do 
	ArterialSpinLabelling[i]="${line_data}"; ((++i)); 
    done < "${output_dir}"/arterial.tmp
    rm -f "${output_dir}"/arterial.tmp
else
    echo "Cannot find arterial.tmp in output dir"
fi

echo "ASL : ${firsttmpdcm} size: " ${#ArterialSpinLabelling[*]}
if [[ ${#ArterialSpinLabelling[*]} -ne 17 ]];then
    echo "DCM modification error. Not enough ASL parameters."
    echo " Ignoring ASL modifications."
else
    echo "FrameAnt 7: " ${ArterialSpinLabelling[7]}

    ${DCMODIFY} -i "(0018,9251)[0].(0018,9250)=${ArterialSpinLabelling[0]}" $files	   #  Arterial Spin Labeling Contrast	      VR=<CS>	VL=<000a>
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9252)=${ArterialSpinLabelling[1]}" $files     # ASL Technique Description       VR=<LO>   VL=<0004>  <FAIR> 
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9257)=${ArterialSpinLabelling[2]}" $files	   #ASL Context	    VR=<CS>   VL=<0006>  <LABEL > 
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9259)=${ArterialSpinLabelling[3]}" $files     #ASL Crusher Flag        VR=<CS>   VL=<0002>  <NO>   
    ${DCMODIFY} -i "(0018,9251)[0].(0018,925b)=${ArterialSpinLabelling[4]}" $files	   #ASL Crusher Description	    VR=<LO>   VL=<0014>  <crusher description > 
    ${DCMODIFY} -i "(0018,9251)[0].(0018,925c)=${ArterialSpinLabelling[5]}" $files     #ASL Bolus Cut-off Flag          VR=<CS>   VL=<0002>  <NO>
    # ASL Bolus Cut-off Timing Sequence
    ${DCMODIFY} -i "(0018,9251)[0].(0018,925d)[0].(0018,925e)=${ArterialSpinLabelling[6]}" $files	#ASL Bolus Cut-off Technique	 VR=<LO>   VL=<0000>	<>	     
    ${DCMODIFY} -i "(0018,9251)[0].(0018,925d)[0].(0018,925f)=${ArterialSpinLabelling[7]}" $files	#ASL Bolus Cut-off Delay Time	 VR=<UL>   VL=<0004>	[00000000]
    # ASL Slab Sequence	
    #    Anatomic Region Sequence
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9260)[0].(0008,2218)[0].(0008,0100)=${ArterialSpinLabelling[8]}" $files  #Code Value	   VR=<SH>   VL=<0008>  <T-D1100 > 
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9260)[0].(0008,2218)[0].(0008,0104)=${ArterialSpinLabelling[9]}" $files  #Code Meaning    VR=<LO>   VL=<0004>  <Head>     
    #    Anatomic Region Modifier Sequence 
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9260)[0].(0008,2220)[0].(0008,0100)=${ArterialSpinLabelling[10]}" $files #Code Value	   VR=<SH>   VL=<0006>  <G-A138>   
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9260)[0].(0008,2220)[0].(0008,0104)=${ArterialSpinLabelling[11]}" $files #Code Meaning    VR=<LO>   VL=<0008>  <Coronal > 

    ${DCMODIFY} -i "(0018,9251)[0].(0018,9260)[0].(0018,9253)=${ArterialSpinLabelling[12]}" $files #ASL Slab Number		VR=<US>	  VL=<0002>  [0001]
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9260)[0].(0018,9254)=${ArterialSpinLabelling[13]}" $files #ASL Slab Thickness      VR=<FD>   VL=<0008>  {12}
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9260)[0].(0018,9255)=${ArterialSpinLabelling[14]}" $files #ASL Slab Orientation	VR=<FD>	  VL=<0018>  {1.38,180.661,87.387}
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9260)[0].(0018,9256)=${ArterialSpinLabelling[15]}" $files #ASL Mid Slab Position   VR=<FD>   VL=<0018>  {0,0,0}
    ${DCMODIFY} -i "(0018,9251)[0].(0018,9260)[0].(0018,9258)=${ArterialSpinLabelling[16]}" $files #ASL Pulse Train Duration	VR=<UL>	  VL=<0004>  [00000000]
    

fi # array check
echo "DCModify ASL done"
exit 0
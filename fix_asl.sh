


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




(0018,9251)[0]. 
  ----:               
    (0018,9251)[0].(0x0018,0x9250) 
    (0018,9251)[0].(0x0018,0x9252) 
    (0018,9251)[0].(0x0018,0x9257) 
    (0018,9251)[0].(0x0018,0x9259) 
    (0018,9251)[0].(0x0018,0x925b) 
    (0018,9251)[0].(0x0018,0x925c) 
    (0018,9251)[0].(0x0018,0x925d) 
  ----:               
    (0018,9251)[0].(0018,925d)[0].(0x0018,0x925e) 
    (0018,9251)[0].(0018,925d)[0].(0x0018,0x925f) 
                      
    (0018,9251)[0].(0x0018,0x9260) 
  ----:               
    (0018,9251)[0].(0018,9260)[0].(0x0018,0x9253) 
    (0018,9251)[0].(0018,9260)[0].(0x0018,0x9254) 
    (0018,9251)[0].(0018,9260)[0].(0x0018,0x9255) 
    (0018,9251)[0].(0018,9260)[0].(0x0018,0x9256) 
    (0018,9251)[0].(0018,9260)[0].(0x0018,0x9258) 

    (0018,9251)[0].(0018,9260)[0].(0x0008,0x2218)[0].(0x0008,0x0100) 
    (0018,9251)[0].(0018,9260)[0].(0x0008,0x2218)[0].(0x0008,0x0104) 

    (0018,9251)[0].(0018,9260)[0].(0x0008,0x2220)[0].(0x0008,0x0100) 
    (0018,9251)[0].(0018,9260)[0].(0x0008,0x2220)[0].(0x0008,0x0104) 
                      
                      



    firsttmpdcm=$(ls -1 ${output_dir}/tmp/*.dcm| head -1)
    dcdump ${firsttmpdcm} 2>&1 >/dev/null | grep '(0x0008,0x010' | awk -F'>' '/</ {print $4}'| tr -d '<' > ${output_dir}/anatomy.tmp
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
    if [[ ${#FrameAnatomySequence[*]} -ne 8 ]];then
	echo "DCM modification error. Not enough Frame Anatomy Sequence parameters."
	echo " Ignoring Anatomy modifications."
    else
	echo "FrameAnt 7: " ${FrameAnatomySequence[7]}
	${DCMODIFY} -i "(0018,9251)[0].(0020,9071)[0].(0008,2218)[0].(0008,0100)=${FrameAnatomySequence[0]}" $files
	${DCMODIFY} -i "(0018,9251)[0].(0020,9071)[0].(0008,2218)[0].(0008,0104)=${FrameAnatomySequence[1]}" $files   #CodeMeaning=
	${DCMODIFY} -i "(0018,9251)[0].(0020,9071)[0].(0008,2218)[0].(0008,0102)=SRT" $files

	${DCMODIFY} -i "(0018,9251)[0].(0020,9071)[0].(0008,2228)[0].(0008,0100)=${FrameAnatomySequence[4]}" $files
	${DCMODIFY} -i "(0018,9251)[0].(0020,9071)[0].(0008,2228)[0].(0008,0104)=${FrameAnatomySequence[5]}" $files
	${DCMODIFY} -i "(0018,9251)[0].(0020,9071)[0].(0008,2228)[0].(0008,0102)=SRT" $files

	${DCMODIFY} -i "(0018,9251)[0].(0020,9071)[0].(0020,9072)=R" $files

    fi # array check
    echo "DCModify ASL done"

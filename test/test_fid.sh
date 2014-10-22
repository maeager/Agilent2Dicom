#!/usr/bin/env bash
#
# Test examples for FID to Dicom (no filtering) 
#
# - Monash Biomedical Imaging
# - Michael Eager (michael.eager@monash.edu)

#   Copyright (C) 2014 Michael Eager  (michael.eager@monash.edu)

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

set -o errexit  # -e
set -o pipefail
set -u  # nounset

FID2DCM=../fid2dicom.py   #../fid2dcm.sh
EXAMPLEDATA=../..
fidfolders=`find $EXAMPLEDATA -type d -name "*.fid"`
echo $fidfolders
fidarray=($fidfolders)
echo "Number of FIDs " ${#fidarray[*]}
for exfid in "${fidarray[@]}"
do
    echo $exfid
    output_path=$(echo ${exfid} | sed -e 's/fid/dcm/' -e 's/example/output/')
    [ $exfid == $output_path ] && continue
    echo "Converting " $exfid " to " $output_path
    $FID2DCM -v -i "$exfid" -o "${output_path}"
    mrinfo "${output_path}"
break
done



# ../example_data/s_2013081903/35tr-2echo-20d-1nex_01.fid
# ../example_data/s_2014051208/_02.fid
# ../example_data/s_2014051208/_01.fid
# ../example_data/s_2014051208/prescan_power.fid
# ../example_data/s_2014051208/diff-trace_01.fid
# ../example_data/s_2014051208/cor-vj4_01.fid
# ../example_data/s_2014051208/sag-vj4_01.fid
# ../example_data/s_2014051208/DTI_EPIP_J30_01.fid
# ../example_data/s_2014051208/axi-vj4_01.fid
# ../example_data/s_2014051208/sag-vj4-cru3_01.fid
# ../example_data/s_2014051208/prescan_freq.fid
# ../example_data/s_2014051208/epip_axi_300TR_01.fid
# ../example_data/s_2014061202/sag-fatsat_01.fid
# ../example_data/s_2014061202/_02.fid
# ../example_data/s_2014061202/_01.fid
# ../example_data/s_2014061202/fse3d-trace-b1500_01.fid
# ../example_data/s_2014061202/fatsat_01.fid
# ../example_data/s_2014061202/cor-fatsat_01.fid
# ../example_data/s_2014061202/fse3d_01.fid
# ../example_data/s_2014061202/axi-fatsat_02.fid
# ../example_data/s_2014061202/prescan_power.fid
# ../example_data/s_2014061202/fse3d-b0_02.fid
# ../example_data/s_2014061202/fse3d-b0_01.fid
# ../example_data/s_2014061202/epi-dir30_01.fid
# ../example_data/s_2014061202/fse3d-trace-b3000_01.fid
# ../example_data/s_2014061202/fatsat-axi_01.fid
# ../example_data/s_2014061202/epi_01.fid
# ../example_data/s_2014061202/epi-dir30_02.fid
# ../example_data/s_2014061202/prescan_freq.fid
# ../example_data/s_2014061202/fse3d-dir30_01.fid
# ../example_data/s_2014061202/fse3d-dir30_02.fid
# ../example_data/s_2014061301/axi-fatsat_01.fid
# ../example_data/s_2014061301/epi-test_01.fid
# ../example_data/s_2014061301/_02.fid
# ../example_data/s_2014061301/epi-dir30_03.fid
# ../example_data/s_2014061301/_01.fid
# ../example_data/s_2014061301/axi-fatsat-longecho_01.fid
# ../example_data/s_2014061301/prescan_power.fid
# ../example_data/s_2014061301/fse3d-b0_02.fid
# ../example_data/s_2014061301/cor_01.fid
# ../example_data/s_2014061301/fse3d-30dirs_01.fid
# ../example_data/s_2014061301/fse3d-b0_06.fid
# ../example_data/s_2014061301/fse3d-b0_04.fid
# ../example_data/s_2014061301/fse3d-b0_01.fid
# ../example_data/s_2014061301/axi_02.fid
# ../example_data/s_2014061301/epi-dir30_01.fid
# ../example_data/s_2014061301/epi-dir30_06.fid
# ../example_data/s_2014061301/fse3d-b0_05.fid
# ../example_data/s_2014061301/axi_03.fid
# ../example_data/s_2014061301/epi-dir30_02.fid
# ../example_data/s_2014061301/prescan_freq.fid
# ../example_data/s_2014061301/fse3d-b0_03.fid
# ../example_data/s_2014061301/axi-fatsat-shortecho_01.fid
# ../example_data/s_2014061301/fse3d-30dirs_02.fid
# ../example_data/s_2014061301/epi-dir30_04.fid
# ../example_data/s_2014061301/sag_01.fid
# ../example_data/s_2014061301/epi-test_02.fid
# ../example_data/s_2014061301/epi-dir30_05.fid
# ../example_data/s_2014061301/axi_01.fid
# ../example_data/s_2014061301/fse3d-dw-readout_01.fid
# ../example_data/s_2014072901/T2-cor_01.fid
# ../example_data/s_2014072901/dti-rep3_01.fid
# ../example_data/s_2014072901/prescan_power.fid
# ../example_data/s_2014072901/cor_01.fid
# ../example_data/s_2014072901/dti-rep2_01.fid
# ../example_data/s_2014072901/dti-rep1_01.fid
# ../example_data/s_2014072901/prescan_freq.fid
# ../example_data/s_2014072901/sag_01.fid
# ../example_data/s_2014072901/dti-rep4_01.fid
# ../example_data/s_2014072901/dti-test_01.fid
# ../example_data/s_2014072901/axi_01.fid
# ../example_data/s_2014080502/T2-cor_01.fid
# ../example_data/s_2014080502/dti-rep3_01.fid
# ../example_data/s_2014080502/prescan_power.fid
# ../example_data/s_2014080502/cor_01.fid
# ../example_data/s_2014080502/dti-rep2_01.fid
# ../example_data/s_2014080502/dti-rep1_01.fid
# ../example_data/s_2014080502/prescan_freq.fid
# ../example_data/s_2014080502/sag_01.fid
# ../example_data/s_2014080502/dti-rep4_01.fid
# ../example_data/s_2014080502/dti-test_01.fid
# ../example_data/s_2014080502/axi_01.fid
# ../example_data/s_2014032106/axi_01.fid
# ../example_data/s_2014032106/_02.fid
# ../example_data/s_2014032106/prescan_freq.fid
# ../example_data/s_2014032106/prescan_power.fid
# ../example_data/s_2014032106/sag_01.fid
# ../example_data/s_2014032106/cor_01.fid
# ../example_data/s_2014032106/3d-part1_01.fid
# ../example_data/s_2014032106/fse-part1_01.fid
# ../example_data/s_2014032106/axi_02.fid
# ../example_data/s_2014032106/sag_02.fid
# ../example_data/s_2014032106/cor_02.fid
# ../example_data/s_2014032106/fse-part2_01.fid
# ../example_data/s_2014032106/fse-part2-hr_01.fid
# ../example_data/s_2014032106/sag_03.fid
# ../example_data/s_2014032106/cor_03.fid
# ../SheepfetusBrain/s_2014082901/mge3d-100um-2_01.fid
# ../SheepfetusBrain/s_2014082901/mge3d-100um-3_01.fid
# ../SheepfetusBrain/s_2014082901/mge3d-100um_01.fid
# ../SheepfetusBrain/s_2014082901/sag_01.fid
# ../SheepfetusBrain/s_2014072501/mge3d-100um_01.fid
# ../SheepfetusBrain/s_2014072501/mge3d-100um_02.fid
# ../SheepfetusBrain/s_2014062002/mge3d-100um_01.fid
# ../SheepfetusBrain/s_2014062002/mge3d-1.33ml_01.fid
# ../SheepfetusBrain/s_2014062002/T2-cor_01.fid
# ../SheepfetusBrain/s_2014062002/_01.fid
# ../SheepfetusBrain/s_2014062002/axi_01.fid
# ../SheepfetusBrain/s_2014062002/cor_01.fid
# ../SheepfetusBrain/s_2014062002/prescan_freq.fid
# ../SheepfetusBrain/s_2014062002/prescan_power.fid
# ../SheepfetusBrain/s_2014062002/sag_01.fid
# ../SheepfetusBrain/s_2014062002/se-axi_01.fid
# ../SheepfetusBrain/s_2014062002/se-cor-test_01.fid
# ../SheepfetusBrain/s_2014062002/se-cor_01.fid
# ../SheepfetusBrain/s_2014062002/se-sag_01.fid
# ../SheepfetusBrain/s_2014062003/_02.fid
# ../SheepfetusBrain/s_2014062003/axi_01.fid
# ../SheepfetusBrain/s_2014062003/cor_01.fid
# ../SheepfetusBrain/s_2014062003/cor_02.fid
# ../SheepfetusBrain/s_2014062003/mge3d-100um_01.fid
# ../SheepfetusBrain/s_2014062003/prescan_freq.fid
# ../SheepfetusBrain/s_2014062003/prescan_power.fid
# ../SheepfetusBrain/s_2014062003/sag_01.fid
# ../SheepfetusBrain/s_2014073102/mge3d-100um_01.fid
# ../ExampleAgilentData/kidney512iso_01.fid

#!/usr/bin/env bash

export PATH="./dicom3tools_1.00.snapshot.20140306142442/bin/1.2.6.35.x8664/:$PATH"

module load dcmtk

dcmdump ./output_data/standard2d

dcmdump  output_data/multiecho2d_magandphase/magnitude.dcm/0001.dcm

dcmdump  output_data/multiecho2d_magandphase/phase.dcm/0001.dcm

dcmdump  output_data/multiecho3d_magnitudeonly/0001.dcm

dcmdump  output_data/cine/0001.dcm

dcmdump  output_data/diffusion/0001.dcm

dcmdump  output_data/ASL_se_06.dcm/0001.dcm



# See http://www.dclunie.com/dicom3tools/dciodvfy.html for documentation

dciodvfy -dump ./output_data/standard2d/0001.dcm 2>&1 >/dev/null | grep Required | uniq


dciodvfy -dump ./output_data/multiecho3d_magonly/0001.dcm 2>&1 >/dev/null | grep Error | uniq

dciodvfy -dump ./output_data/multiecho2d_magandphase/magnitude.dcm/0001.dcm 2>&1 >/dev/null | grep Error | uniq

dciodvfy -dump ./output_data/multiecho2d_magandphase/phase.dcm/0001.dcm 2>&1 >/dev/null | grep Error | uniq 

dciodvfy -dump ./output_data/cine/0001.dcm 2>&1 >/dev/null | grep Error | uniq

dciodvfy -dump ./output_data/ASL_se_06.dcm/0001.dcm 2>&1 | grep Required | uniq

 dciodvfy -dump ./output_data/diffusion/0001.dcm 2>&1 >/dev/null | grep Required | uniq

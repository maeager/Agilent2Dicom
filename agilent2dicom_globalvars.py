#!/usr/bin/env python
# pylint: disable=bad-whitespace
# agilent2dicom is used to convert Agilent FDF and FID images to DICOM format.

  # Copyright (C) 2014 Michael Eager  (michael.eager@monash.edu)

  # This program is free software: you can redistribute it and/or modify
  # it under the terms of the GNU General Public License as published by
  # the Free Software Foundation, either version 3 of the License, or
  # (at your option) any later version.

  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.

  # You should have received a copy of the GNU General Public License
  # along with this program.  If not, see <http://www.gnu.org/licenses/>.


### Ensure no spaces in variable declaration so that shell and python comply ###
AGILENT2DICOM_VERSION="1.7.0"
AGILENT2DICOM_APP_VERSION="1.8.2"
FDF2DCMVERSION="1.2"
FID2DCMVERSION="1.6"
DVCS_STAMP="$Id: agilent2dicom_globalvars.py,v 473d4cda6b9d 2015/03/04 23:24:19 michael $"

UID_ROOT="1.3.6.1.4.1" # Agilent Root UID 1.3.6.1.4.1, default "2.25"
UID_Type_InstanceCreator="0"
UID_Type_MediaStorageSOPInstance="1"
UID_Type_StudyInstance="2"
UID_Type_SeriesInstance="3"
UID_Type_FrameOfReference="4"
UID_Type_DimensionIndex1="5"
UID_Type_DimensionIndex2="6"

# Hard coded DICOM tag values
DICOM_Tag_Manufacturer="Agilent Technologies"
DICOM_Tag_InstitutionName="Monash Biomedical Imaging"
DICOM_Tag_ManufacturerModelName="vnmrs"
DICOM_Tag_DeviceSerialNumber="unknown"
DICOM_Tag_SoftwareVersions="VnmrJ 3.2, MBI's Agilent2Dicom "

Standard_MR_SOPClassUID="1.2.840.10008.5.1.4.1.1.4" # MR Image SOP
Implementation_Class_UID="1.3.6.1.4.1.25371.1.1.2"

Derivation_Description="Dicom generated from Agilent2Dicom, an FDF/FID dataset converter to enhanced Dicom format. Monash Biomedical Imaging. Copyright 2014 Michael Eager. "
COIL_Manufacturer="Agilent Technologies"
FDF2DCM_Image_Comments="MBI's FDF2DCM converter."
FID2DCM_Image_Comments="MBI's FID2DCM converter."

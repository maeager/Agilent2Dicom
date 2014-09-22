#!/usr/bin/env python

# agilent2dicom is used to convert Agilent FDF files to DICOM format.
# (c) 2014 Michael Eager  (michael.eager@monash.edu)


VersionNumber="1.1.0"
DVCSstamp="$Id$"

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
DICOM_Tag_SoftwareVersions="VnmrJ 3.2, Agilent2Dicom 1.1"

Standard_MR_SOPClassUID="1.2.840.10008.5.1.4.1.1.4" # MR Image SOP
Implementation_Class_UID="1.3.6.1.4.1.25371.1.1.2"

Derivation_Description="Dicom generated from Agilen2Dicom, an FDF dataset converter to Dicom. Monash Biomedical Imaging, Imaging Team."
COIL_Manufacturer="Agilent Technologies"
FDF2DCM_Image_Comments="MBI's FDF2DCM converter."
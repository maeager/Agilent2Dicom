#!/usr/bin/env python

"""ParseFDF
  ParseFDF, and the accompanying methods ParseDiffusionFDF and ParseASLFDF,
  are used to parse Agilent FDF image files
  and correct errors in the Dicom dataset produced by procpar


  Copyright (C) 2014 Michael Eager  (michael.eager@monash.edu)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
import re
import math
import numpy
import argparse

from dicom.sequence import Sequence
from dicom.dataset import Dataset


import agilent2dicom_globalvars as A2D
import ProcparToDicomMap


def AssertImplementation(testval, fdffilename, comment, assumption):
    """ASSERTIMPLEMENTATION - Check FDF properties match up with interpretation
    of procpar Due to lack of documentation on the use of the procpar file, the
    interpretation implemented in this script is based on various documents and
    scripts and may have errors.  This function seeks to double check some of the
    interpretations against the fdf properties.

    """
    if testval:
        if len(fdffilename) > 0:
            fdfstr = "fdf file: " + fdffilename + "\n"
        else:
            fdfstr = ""

        print "\nImplementation check error:\n" + fdfstr + comment + '\nAssumption:' + assumption + '\n'
        # sys.exit(1)


def ParseDiffusionFDF(ds, procpar, fdf_properties, args):
    """ParseDiffusionFDF

    :param ds: Dicom dataset
    :param procpar: Procpar dictionary tag/value pairs
    :param fdf_properties: Tag/value pairs of local fdf file

    :param args: Input arguments
    :returns: Dicom struct
    """
    if args.verbose:
        print 'Processing diffusion image'

    # Get procpar diffusion parameters
    bvalue = procpar['bvalue']  # 64 element array
    bvaluesortidx = numpy.argsort(bvalue)
    # not used -- BvalSave = procpar['bvalSave']
    # if 'bvalvs' in procpar.keys():
    #    BvalVS = procpar['bvalvs']
    # excluded in external recons by vnmrj, unused here
    bvalueRS = procpar['bvalrs']  # 64
    bvalueRR = procpar['bvalrr']  # 64
    bvalueRP = procpar['bvalrp']  # 64
    bvaluePP = procpar['bvalpp']  # 64
    bvalueSP = procpar['bvalsp']  # 64
    bvalueSS = procpar['bvalss']  # 64

    if procpar['recon'] == 'external':
        diffusion_idx = 0
        while True:
            if math.fabs(bvalue[diffusion_idx] -
                         fdf_properties['bvalue']) < 0.005:
                break
            diffusion_idx += 1
        # diffusion_idx = fdf_properties['array_index'] - 1
    else:
        diffusion_idx = fdf_properties['array_index'] * 2

    if diffusion_idx > len(bvalue):
        print '''Procpar bvalue does not contain enough values
        determined by fdf_properties array_index'''

    if args.verbose:
        print 'Diffusion index ', diffusion_idx, ' arrary index ',
        fdf_properties['array_index']

    # Sort diffusion based on sorted index of bvalue instead of
    # fdf_properties['array_index']
    ds.AcquisitionNumber = bvaluesortidx[diffusion_idx]

    if math.fabs(bvalue[diffusion_idx] - fdf_properties['bvalue']) > 0.005:
        print 'Procpar and fdf B-value mismatch: procpar value ',
        bvalue[diffusion_idx], ' and  local fdf value ',
        fdf_properties['bvalue'], ' array idx ', fdf_properties['array_index']

    # MR Diffusion Sequence (0018,9117) see DiffusionMacro.txt
    # B0 scan does not need the MR Diffusion Gradient Direction Sequence macro
    # and its directionality should be set to NONE the remaining scans relate
    # to particular directions hence need the direction macro
    diffusionseq = Dataset()
    if fdf_properties['bvalue'] < 20:
        diffusionseq.DiffusionBValue = 0
        diffusionseq.DiffusionDirectionality = 'NONE'
    else:
        diffusionseq.DiffusionBValue = int(fdf_properties['bvalue'])
        # TODO  One of: DIRECTIONAL,  BMATRIX, ISOTROPIC, NONE
        diffusionseq.DiffusionDirectionality = 'BMATRIX'

    # Diffusion Gradient Direction Sequence (0018,9076)
        diffusiongraddirseq = Dataset()
        # Diffusion Gradient Orientation  (0018,9089)
        # diffusiongraddirseq.add_new((0x0018,0x9089), 'FD',[
        # fdf_properties['dro'],  fdf_properties['dpe'],
        # fdf_properties['dsl']])
        diffusiongraddirseq.DiffusionGradientOrientation = [
            fdf_properties['dro'], fdf_properties['dpe'],
            fdf_properties['dsl']]
        diffusionseq.DiffusionGradientDirectionSequence = Sequence(
            [diffusiongraddirseq])
        # diffusionseq.add_new((0x0018,0x9076), 'SQ',
        # Sequence([diffusiongraddirseq]))

    # Diffusion b-matrix Sequence (0018,9601)
        diffbmatseq = Dataset()
        diffbmatseq.DiffusionBValueXX = bvalueRR[diffusion_idx]
        diffbmatseq.DiffusionBValueXY = bvalueRP[diffusion_idx]
        diffbmatseq.DiffusionBValueXZ = bvalueRS[diffusion_idx]
        diffbmatseq.DiffusionBValueYY = bvaluePP[diffusion_idx]
        diffbmatseq.DiffusionBValueYZ = bvalueSP[diffusion_idx]
        diffbmatseq.DiffusionBValueZZ = bvalueSS[diffusion_idx]
        diffusionseq.DiffusionBMatrixSequence = Sequence([diffbmatseq])

    # TODO  One of: FRACTIONAL, RELATIVE, VOLUME_RATIO
    diffusionseq.DiffusionAnisotropyType = 'FRACTIONAL'
    ds.MRDiffusionSequence = Sequence([diffusionseq])

    MRImageFrameType = Dataset()
    MRImageFrameType.FrameType = [
        "ORIGINAL", "PRIMARY", "DIFFUSION", "NONE"]  # same as ds.ImageType
    MRImageFrameType.PixelPresentation = ["MONOCHROME"]
    MRImageFrameType.VolumetrixProperties = ["VOLUME"]
    MRImageFrameType.VolumeBasedCalculationTechnique = ["NONE"]
    MRImageFrameType.ComplexImageComponent = ["MAGNITUDE"]
    MRImageFrameType.AcquisitionContrast = ["DIFFUSION"]
    ds.MRImageFrameTypeSequence = Sequence([MRImageFrameType])

    return ds


def ParseASL(ds, procpar, fdf_properties):
    # (0018,9257)	    1C	The purpose of the Arterial Spin Labeling.
    #   Enumerated Values:
    #	       LABEL
    #              CONTROL
    #		M_ZERO_SCAN
    #	Required if Frame Type (0008,9007) is
    #	ORIGINAL. May be present otherwise.
    #	See C.8.13.5.14.1 for further
    #	explanation.
    if fdf_properties["asltag"] == 1:
        ds.MRArterialSpinLabeling[0].ASLContext = 'LABEL'
    elif fdf_properties["asltag"] == -1:
        ds.MRArterialSpinLabeling[0].ASLContext = 'CONTROL'
    else:
        ds.MRArterialSpinLabeling[0].ASLContext = 'M_ZERO_SCAN'

    # FIX ME : this could be either array_index or slice_no
    ds.MRArterialSpinLabeling[0].ASLSlabSequence[
        0].ASLSlabNumber = fdf_properties["array_index"]

    # ASL Mid slab position

    # The Image Plane Attributes, in conjunction with the Pixel Spacing
    # Attribute, describe the position and orientation of the image slices
    # relative to the patient-based coordinate system. In each image frame the
    # Image Position (Patient) (0020,0032) specifies the origin of the image
    # with respect to the patient-based coordinate system. RCS and the Image
    # Orientation (Patient) (0020,0037) attribute values specify the
    # orientation of the image frame rows and columns. The mapping of pixel
    # location i, j to the RCS is calculated as follows:
    #					 X x i Yx j 0 S x
    #				Px				  i	   i
    #					 X y i Yy j 0 S y
    #				Py				  j	   j
    #								    =M
    #					 X z i Yz j 0 S z	  0	   0
    #				Pz
    #				1	   0	   0	  01	  1	   1
    # Where:
    #	 Pxyz The coordinates of the voxel (i,j) in the frame's
    # image plane in units of mm.
    #	 Sxyz The three values of the Image Position (Patient) (0020,0032)
    # attributes. It is the
    #	       location in mm from the origin of the RCS.
    #	 Xxyz The values from the row (X) direction cosine of the Image
    # Orientation (Patient)
    #	       (0020,0037) attribute.
    #	 Yxyz The values from the column (Y) direction cosine of the Image
    # Orientation (Patient)
    #	       (0020,0037) attribute.
    #	 i     Column index to the image plane. The first column is index zero.
    #	   i Column pixel resolution of the Pixel Spacing (0028,0030)
    # attribute in units of mm.
    #	 j     Row index to the image plane. The first row index is zero.
    # j Row pixel resolution of the Pixel Spacing (0028,0030) attribute in
    # units of mm.

    # ds.MRArterialSpinLabelingSequence.ASLSlabSequence[0].ASLMidSlabPosition
    # = [str(ImagePositionPatient[0]), str(ImagePositionPatient[1]),
    # str(ImagePositionPatient[2] + (islice-1)*SliceThickness))]

    #            print ImagePositionPatient
    #           M = numpy.matrix([[PixelSpacing[0] *
    #           ImageOrientationPatient[0], PixelSpacing[1] *
    #           ImageOrientationPatient[1], SliceThinkness *
    #           ImageOrientationPatient[2] ImagePositionPatient[0]],
    #           [PixelSpacing[0] * ImageOrientationPatient[3], PixelSpacing[1]
    #           * ImageOrientationPatient[4], SliceThinkness *
    #           ImageOrientationPatient[5] ImagePositionPatient[1]],
    #           [PixelSpacing[0] * ImageOrientationPatient[6], PixelSpacing[1]
    #           * ImageOrientationPatient[8], SliceThinkness *
    #           ImageOrientationPatient[8] ImagePositionPatient[2]], [0, 0, 0,
    #           1]])
    #           pos = numpy.matrix([[ceil(ds.Rows / 2)],[ ceil(ds.Columns /
    #           2],[fdf_properties['slice_no'],[1]])
    #           Pxyz = M * pos

    # ds.MRArterialSpinLabelingSequence.ASLSlabSequence[0].ASLMidSlabPosition
    # = [str(Pxyz[0,0]),str(Pxyz[1,0]),str(Pxyz[2,0]))]
    return ds


def ParseFDF(ds, fdf_properties, procpar, args):
    """
    ParseFDF modify the dicom dataset structure based on FDF
    header information.

    Comment text copied from VNMRJ Programming.pdf

    :param ds:       Dicom dataset
    :param fdf_properties: Dict of fdf header label/value pairs
    :param procpar:  Dict of procpar label/value pairs
    :param args:     Argparse object
    :return ds:      Return updated dicom dataset struct
    :return fdfrank: Number of dimensions (rank) of fdf file
    """

    # if procpar['recon'] == 'external' and fdf_properties['rank'] == '3'
    # and procpar:
    #     fdf_tmp = fdf_properties['roi']
    #     fdf_properties['roi'][0:1] = fdf_tmp[1:2]
    #     fdf_properties['roi'][2] = fdf_tmp[0]
    #     fdf_tmp = fdf_properties['matrix']
    #     fdf_properties['matrix'][0:1] = fdf_tmp[1:2]
    #     fdf_properties['matrix'][2] = fdf_tmp[0]

    #----------------------------------------------------------
    # General implementation checks
    filename = fdf_properties['filename']

    # File dimensionality or Rank fields
    # rank is a positive integer value `(1, 2, 3, 4,...) giving the
    # number of dimensions in the data file (e.g., int rank=2;).
    fdfrank = fdf_properties['rank']
    acqndims = procpar['acqdim']
    CommentStr = '''Acquisition dimensionality (ie 2D or 3D) does not
    match between fdf and procpar'''
    AssumptionStr = '''Procpar nv2 > 0 indicates 3D acquisition and
    fdf rank property indicates dimensionality.\n''' +\
        'Using local FDF value ' + \
        str(fdfrank) + ' instead of procpar value ' + str(acqndims) + '.'
    if args.verbose:
        print 'Acqdim (type): ' + ds.MRAcquisitionType + " acqndims " + str(acqndims)

    AssertImplementation(
        acqndims != fdfrank, filename, CommentStr, AssumptionStr)

    # matrix is a set of rank integers giving the number of data
    # points in each dimension (e.g., for rank=2, float
    # matrix[]={256,256};)
    if fdfrank == 3:
        fdf_size_matrix = fdf_properties['matrix'][0:3]
    else:
        fdf_size_matrix = fdf_properties['matrix'][0:2]
    if args.verbose:
        print "FDF size matrix ", fdf_size_matrix, type(fdf_size_matrix)
    #fdf_size_matrix = numpy.array(fdf_matrix)

    # spatial_rank is a string ("none", "voxel", "1dfov", "2dfov",
    # "3dfov") for the type of data (e.g., char
    # *spatial_rank="2dfov";).
    spatial_rank = fdf_properties['spatial_rank']

    #  0018,0023 MR Acquisition Type (optional)
    # Identification of spatial data encoding scheme.
    # Defined Terms: 1D 2D 3D
    fdf_MRAcquisitionType = '2D'
    if spatial_rank == "3dfov":
        fdf_MRAcquisitionType = '3D'
    CommentStr = 'MR Acquisition type does not match between fdf and procpar'
    AssumptionStr = '''In fdf, MR Acquisition type defined by spatial_rank
    and matrix. \n
        For 2D, spatial_rank="2dfov" and matrix has two elements eg.
    {256,256}. \n
        For 3D, spatial_rank="3dfov" and matrix has three elements.\n
        In procpar, MR Acquisition type is defined by nv2 > 0 or lpe2 > 0.\n
        Using local FDF value ''' + fdf_MRAcquisitionType + \
        ' instead of procpar value ' + ds.MRAcquisitionType + '.'
    AssertImplementation(
        ds.MRAcquisitionType != fdf_MRAcquisitionType, filename, CommentStr,
        AssumptionStr)
    ds.MRAcquisitionType = fdf_MRAcquisitionType

    # Data Content Fields
    # The following entries define the data type and size.
    #  - storage is a string ("integer", "float") that defines the data
    # type (e.g., char *storage="float";).
    #  - bits is an integer (8, 16, 32, or 64) that defines the size of the
    # data (e.g., float bits=32;).
    # - type is a string ("real", "imag", "absval", "complex") that defines the
    # numerical data type (e.g., char *type="absval";).

    # roi is the size of the acquired data volume (three floating
    # point values), in centimeters, in the user's coordinate frame,
    # not the magnet frame (e.g., float roi[]={10.0,15.0,0.208};). Do
    # not confuse this roi with ROIs that might be specified inside
    # the data set.
    if fdfrank == 3:
        roi = fdf_properties['roi'][0:3]
    else:
        roi = fdf_properties['roi'][0:2]
    if args.verbose:
        print "FDF roi ", roi, type(roi)
    #roi = numpy.array(roi_text)

    # PixelSpacing - 0028,0030 Pixel Spacing (mandatory)
    PixelSpacing = map(lambda x, y: x * 10.0 / y, roi, fdf_size_matrix)
    if PixelSpacing[0] != ds.PixelSpacing[0] or \
       PixelSpacing[1] != ds.PixelSpacing[1]:
        print "Pixel spacing mismatch, procpar ", ds.PixelSpacing, " fdf spacing ", str(PixelSpacing[0]), ', ', str(PixelSpacing[1])
    if args.verbose:
        print "Pixel Spacing : Procpar   ", ds.PixelSpacing
        print "Pixel Spacing : FDF props ", PixelSpacing
    # (0028,0030) Pixel Spacing
    ds.PixelSpacing = [str(PixelSpacing[0]), str(PixelSpacing[1])]

    # FDF slice thickness
    if fdfrank == 3:
        fdfthk = fdf_properties['roi'][2] / fdf_properties['matrix'][2] * 10
    else:
        fdfthk = fdf_properties['roi'][2] * 10.0

    CommentStr = 'Slice thickness does not match between fdf and procpar'
    AssumptionStr = '''In fdf, slice thickness defined by roi[2] for 2D or
    roi[2]/matrix[2].\n
        In procpar, slice thickness defined by thk (2D) or lpe2*10/(fn2/2) or
    lpe2*10/nv2.\n
        Using local FDF value ''' + str(fdfthk) + ' instead of procpar value ' + str(ds.SliceThickness) + '.'
    if args.verbose:
        print 'fdfthk : ' + str(fdfthk)
        print 'SliceThinkness: ' + str(ds.SliceThickness)

    SliceThickness = float(ds.SliceThickness)

    # fix me Quick hack to avoid assert errors for diffusion and 3D magnitude
    # images
    # if not ('diff' in procpar.keys() and procpar["diff"] == 'y'):
    #	 if MRAcquisitionType == '3D':
    #	     print 'Not testing slicethickness in diffusion and 3D MR FDFs'
    #	else:
    AssertImplementation(
        SliceThickness != fdfthk, filename, CommentStr, AssumptionStr)

    # Slice Thickness 0018,0050 Slice Thickness (optional)
    if fdfrank == 3:
        if len(PixelSpacing) != 3:
            print "Slice thickness: 3D procpar spacing not available"
            print " fdfthk ", fdfthk
        else:
            if PixelSpacing[2] != ds.SliceThickness:
                print "Slice Thickness mismatch, procpar ", ds.SliceThickness, " fdf spacing ", PixelSpacing[2], fdfthk

    # Force slice thickness to be from fdf props
    ds.SliceThickness = str(fdfthk)
    SliceThickness = fdfthk

    #-------------------------------------------------------------------------
    # GROUP 0020: Relationship

    ds.ImageComments = A2D.FDF2DCM_Image_Comments + \
        '\n' + fdf_properties['filetext']

    orientation = numpy.array(fdf_properties['orientation']).reshape(3, 3)
    location = numpy.array(fdf_properties['location']) * 10.0
    span = numpy.array(numpy.append(fdf_properties['span'], 0) * 10.0)

    if args.verbose:
        print "FDF Span: ", span, span.shape
        print "FDF Location: ", location, location.shape

    ds, ImageTransformationMatrix = ProcparToDicomMap.CalcTransMatrix(
        ds, orientation, location, span, fdfrank, PixelSpacing, SliceThickness)

    # Nuclear Data Fields
    # Data fields may contain data generated by
    # interactions between more than one nucleus (e.g., a 2D chemical shift
    # correlation map between protons and carbon). Such data requires
    # interpreting the term ppm for the specific nucleus, if ppm to frequency
    # conversions are necessary, and properly labeling axes arising from
    # different nuclei. To properly interpret ppm and label axes, the identity
    # of the nucleus in question and the corresponding nuclear resonance
    # frequency are needed. These fields are related to the abscissa values
    # "ppm1", "ppm2", and "ppm3" in that the 1, 2, and 3 are indices into the
    # nucleus and nucfreq fields. That is, the nucleus for the axis with
    # abscissa string "ppm1" is the first entry in the nucleus field.  -
    # nucleus is one entry ("H1", "F19", same as VNMR tn parameter) for each rf
    # channel (e.g., char *nucleus[]={"H1","H1"};).  - nucfreq is the nuclear
    # frequency (floating point) used for each rf channel (e.g., float
    # nucfreq[]={200.067,200.067};).

    if fdf_properties['nucleus'][0] != ds.ImagedNucleus:
        print 'Imaged nucleus mismatch: ', fdf_properties['nucleus'], \
            ds.ImagedNucleus
    if math.fabs(fdf_properties['nucfreq'][0] -
                 float(ds.ImagingFrequency)) > 0.01:
        print 'Imaging frequency mismatch: ', fdf_properties['nucfreq'], \
            ds.ImagingFrequency

    # Change patient position and orientation in
    # if procpar['recon'] == 'external' and fdf_properties['rank'] == '3':

    #-------------------------------------------------------------------------
    # GROUP 0028: Image Presentation
    # A good short description of this section can be found here:
    # http://dicomiseasy.blogspot.com.au/2012/08/chapter-12-pixel-data.html

        # Implementation check
    CommentStr = 'Number of rows does not match between fdf and procpar'
    AssumptionStr = '''In FDF, number of rows is defined by
        matrix[1]. \n In procpar, for 3D datasets number of rows is
        either fn1/2 or nv (%s ,%s).\n For 2D datasets, number of
        rows is fn/2.0 or np (%s , %s).\n Using local FDF value %s
        instead of procpar value %s.
    ''' % (str(procpar['fn1'] / 2.0), str(procpar['nv']),
           str(procpar['fn'] / 2.0), str(procpar['np']),
           str(fdf_properties['matrix'][1]), str(ds.Rows))
    AssertImplementation(
        int(float(ds.Rows)) != int(fdf_properties['matrix'][1]),
        filename, CommentStr, AssumptionStr)
    if args.verbose:
        print 'Rows ', procpar['fn'] / 2.0, procpar['fn1'] / 2.0, \
            procpar['nv'], procpar['np'] / 2.0
        print '   Procpar: rows ', ds.Rows
        print '   FDF prop rows ', fdf_properties['matrix'][1]
    ds.Rows = fdf_properties['matrix'][1]  # (0028,0010) Rows

    # Implementation check
    CommentStr = 'Number of columns does not match between fdf and procpar'
    AssumptionStr = '''In FDF, number of columns is defined by
    matrix[0]. \n In procpar, for 3D datasets number of columns is
    either fn/2 or np (%s,%s).\n For 2D datasets, number of rows is
    fn1/2.0 or nv (%s ,%s).\n Using local FDF value %s instead of
    procpar value %s.
    ''' % (str(procpar['fn'] / 2.0), str(procpar['np']),
           str(procpar['fn1'] / 2.0), str(procpar['nv']),
           str(fdf_properties['matrix'][0]), str(ds.Columns))
    AssertImplementation(
        int(float(ds.Columns)) != int(fdf_properties['matrix'][0]),
        filename, CommentStr, AssumptionStr)
    if args.verbose:
        print 'Columns ', procpar['fn'] / 2.0, procpar['fn1'] / 2.0, \
            procpar['nv'], procpar['np'] / 2.0, fdf_properties['matrix'][0]
        print '   Procpar: Cols ', ds.Rows
        print '   FDF prop Cols ', fdf_properties['matrix'][0]
    ds.Columns = fdf_properties['matrix'][0]  # (0028,0011) Columns

    #-------------------------------------------------------------------------
    # Number of frames
    # DICOMHDR code:
    #	  elseif $tag='(0028,0008)' then	" no of frames "
    #	    $dim = 2  "default 2D"
    #	    exists('nv2','parameter'):$ex
    #	    if($ex > 0) then
    #	      if(nv2 > 0) then
    #		on('fn2'):$on		"3D data"
    #		if($on) then
    #		  $pe2 = fn2/2.0
    #		else
    #		  $pe2 = nv2
    #		endif
    #		$dim = 3
    #	      endif
    #	    endif
    #
    #	    if ($dim = 3) then
    #	      $f = $pe2	   "no of frames for 3D"
    #	    else
    #	      substr(seqcon,3,1):$spe1
    #	      if($spe1 = 's') then
    #		$f = (ns * (arraydim/nv) * ne)	 "sems type"
    #	      else
    #		$f = (ns * arraydim * ne)	"compressed gems type"
    #	      endif
    #	      if($imagesout='single') then
    #		$f = $f	"single image output: frames=(no_of_slices * \
    #    array_size * ne)"
    #	      else
    #		$f = 1				" single frame"
    #	      endif
    #	    endif
    #	   $fs=''
    #	    format($f,0,0):$fs
    #	    $value='['+$fs+']'
    #	    if $DEBUG then write('alpha','    new value = "%s"',$value) endif

    # if fdfrank == 3:
    #	 ds.NumberOfFrames = fdf_properties['matrix'][2]

    # dicom3tool uses frames to create enhanced MR
    # ds.NumberOfFrames = fdf_properties['slices']
    # ds.FrameAcquisitionNumber = fdf_properties['slice_no']

#    if 'ne' in procpar.keys() and procpar['ne'] > 1:
#	 print 'Processing multi-echo sequence image'

    if 'echo_no' in fdf_properties.keys():
        volume = fdf_properties['echo_no']

    if (len(ds.ImageType) >= 3 and ds.ImageType[2] == "MULTIECHO") and \
       ('echoes' in fdf_properties.keys() and fdf_properties['echoes'] > 1):
        print 'Multi-echo sequence'
        # TE 0018,0081 Echo Time (in ms) (optional)
        if 'TE' in fdf_properties.keys():
            if ds.EchoTime != str(fdf_properties['TE']):
                print "Echo Time mismatch: ", ds.EchoTime, fdf_properties['TE']
            ds.EchoTime = str(fdf_properties['TE'])
        # 0018,0086 Echo Number (optional)
    if 'echo_no' in fdf_properties.keys():
        if ds.EchoNumber != fdf_properties['echo_no']:
            print "Echo Number mismatch: ", ds.EchoNumber,\
                fdf_properties['echo_no']
        ds.EchoNumber = fdf_properties['echo_no']

    if len(ds.ImageType) >= 3 and ds.ImageType[2] == "ASL":
        ds = ParseASL(ds, procpar, fdf_properties)

    # if 'echoes' in fdf_properties.keys() and fdf_properties['echoes'] > 1 \
    #    and fdf_properties['array_dim'] == 1:
    #    ds.AcquisitionNumber = fdf_properties['echo_no']
    #    ds.ImagesInAcquisition = fdf_properties['echoes']
    # else:

    ds.AcquisitionNumber = fdf_properties['array_index']
    if 'array_dim' in fdf_properties.keys():
        ds.ImagesInAcquisition = fdf_properties['array_dim']
    else:
        ds.ImagesInAcquisition = 1

    if len(ds.ImageType) >= 3 and ds.ImageType[2] == 'DIFFUSION':
        ds = ParseDiffusionFDF(ds, procpar, fdf_properties, args)

    # Multi dimension Organisation and Index module
    DimOrgSeq = Dataset()
    # ds.add_new((0x0020,0x9164), 'UI', DimensionOrganizationUID)

    # or SEQUENCE == "Diffusion":
    if len(ds.ImageType) >= 3 and ds.ImageType[2] == "MULTIECHO":
        DimensionOrganizationUID = [ProcparToDicomMap.CreateUID(
            A2D.ID_Type_DimensionIndex1, [], [],
            args.verbose), ProcparToDicomMap.CreateUID(
                A2D.UID_Type_DimensionIndex2, [], [],
                args.verbose)]
        DimOrgSeq.add_new((0x0020, 0x9164), 'UI', DimensionOrganizationUID)
        ds.DimensionOrganizationType = '3D_TEMPORAL'  # or 3D_TEMPORAL
    else:
        DimensionOrganizationUID = ProcparToDicomMap.CreateUID(
            A2D.UID_Type_DimensionIndex1, [], [], args.verbose)
        # if args.verbose:
        #    print "DimUID", DimensionOrganizationUID
        DimOrgSeq.add_new((0x0020, 0x9164), 'UI', [DimensionOrganizationUID])
        ds.DimensionOrganizationType = '3D'  # or 3D_TEMPORAL

    ds.DimensionOrganizationSequence = Sequence([DimOrgSeq])

    if len(ds.ImageType) >= 3 and ds.ImageType[2] == 'MULTIECHO':
        DimIndexSeq1 = Dataset()
        # Image position patient 20,32 or 20,12
        DimIndexSeq1.DimensionIndexPointer = (0x0020, 0x0032)

        # #DimIndexSeq1.DimensionIndexPrivateCreator=
        # #DimIndexSeq1.FunctionalGroupPointer=
        # #DimIndexSeq1.FunctionalGroupPrivateCreator=
        DimIndexSeq1.add_new(
            (0x0020, 0x9164), 'UI', DimOrgSeq.DimensionOrganizationUID[0])
        DimIndexSeq1.DimensionDescriptionLabel = 'Third Spatial dimension'

        DimIndexSeq2 = Dataset()
        DimIndexSeq2.DimensionIndexPointer = (0x0018, 0x0081)  # Echo Time
        # DimIndexSeq2.DimensionIndexPrivateCreator=
        # DimIndexSeq2.FunctionalGroupPointer=
        # DimIndexSeq2.FunctionalGroupPrivateCreator=
        DimIndexSeq2.add_new(
            (0x0020, 0x9164), 'UI', DimOrgSeq.DimensionOrganizationUID[1])
        DimIndexSeq2.DimensionDescriptionLabel = 'Fourth dimension (multiecho)'
        ds.DimensionIndexSequence = Sequence([DimIndexSeq2, DimIndexSeq1])
    else:
        DimIndexSeq1 = Dataset()
        # Image position patient 20,32 or 20,12
        DimIndexSeq1.DimensionIndexPointer = (0x0020, 0x0032)
        # #DimIndexSeq1.DimensionIndexPrivateCreator=
        # #DimIndexSeq1.FunctionalGroupPointer=
        # #DimIndexSeq1.FunctionalGroupPrivateCreator=
        DimIndexSeq1.add_new(
            (0x0020, 0x9164), 'UI', [DimensionOrganizationUID])
        DimIndexSeq1.DimensionDescriptionLabel = 'Third Spatial dimension'
        ds.DimensionIndexSequence = Sequence([DimIndexSeq1])

        # Module: Image Pixel (mandatory)
        # Reference: DICOM Part 3: Information Object Definitions C.7.6.3
        # ds.Rows                     # 0028,0010 Rows (mandatory)
        # ds.Columns                  # 0028,0011 Columns (mandatory)
        # ds.BitsStored               # 0028,0101 (mandatory)
        # ds.HighBit                  # 0028,0102 (mandatory)
        # ds.PixelRepresentation# 0028,0103 Pixel Representation (mandatory)
        # ds.PixelData               #
        # 7fe0,0010 Pixel Data (mandatory)

    FrameContentSequence = Dataset()
    # FrameContentSequence.FrameAcquisitionNumber = '1'
    # fdf_properties['slice_no']
    # FrameContentSequence.FrameReferenceDateTime
    # FrameContentSequence.FrameAcquisitionDateTime
    # FrameContentSequence.FrameAcquisitionDuration
    # FrameContentSequence.CardiacCyclePosition
    # FrameContentSequence.RespiratoryCyclePosition
    # FrameContentSequence.DimensionIndexValues = 1 #islice
    # --- fdf_properties['array_no']
    # FrameContentSequence.TemporalPositionIndex = 1
    FrameContentSequence.StackID = [str(1)]  # fourthdimid
    FrameContentSequence.InStackPositionNumber = [int(1)]  # fourthdimindex
    FrameContentSequence.FrameComments = fdf_properties['filetext']
    FrameContentSequence.FrameLabel = 'DimX'
    ds.FrameContentSequence = Sequence([FrameContentSequence])

    return ds, fdfrank, fdf_size_matrix, ImageTransformationMatrix


def Save3dFDFtoDicom(ds, procpar, image_data, fdf_properties,
                     transformation_matrix, args, outdir, filename):
    """
    Multi-dimension (3D) export of FDF format image data and metadata to DICOM

    :param ds:  Baseline pydicom dataset
    :param image_data: Pixel data of image
    :param fdf_matsize: Matrix size of image data e.g. {256,256,256}
    :param M: Image transformationm matrix
    :param args: argparse struct
    :param outdir: output directory
    :return ds: Return the pydicom dataset
    """

    print "3D export"
    voldata = numpy.reshape(image_data, fdf_properties['matrix'])

    # if procpar['recon'] == 'external':
    #
    #        pdb.set_trace()
    if procpar['recon'] == 'external' and fdf_properties['rank'] == 3:
        if procpar['seqfil'] == "epip":
            print "Transposing external recon 3D"
            voldata = numpy.transpose(voldata, (1, 2, 0))  # 1,2,0
        if procpar['seqfil'] == "fse3d":
            print "Transposing external recon 3D"
            voldata = numpy.transpose(voldata, (2, 0, 1))  # 0,2,1 works
    # readpp.m procpar('nD') == 3
    #        acq.FOVcm = [pps.lro pps.lpe pps.lpe2];
    #        acq.dims = [pps.nf pps.np/2 pps.nv2];
    #        acq.voxelmm = acq.FOVcm./acq.dims*10;

    print "Image data shape: ", str(image_data.shape)
    print "Vol data shape: ", voldata.shape
    print "fdf properties matrix: ", fdf_properties['matrix']
    print "Slice points: ", fdf_properties['matrix'][0] * fdf_properties['matrix'][1]
    #  slice_data = numpy.zeros_like(numpy.squeeze(image_data[:,:,1]))
    #   if 'ne' in procpar.keys():

    range_max = fdf_properties['matrix'][2]
    num_slicepts = fdf_properties['matrix'][0] * fdf_properties['matrix'][1]
    if procpar['recon'] == 'external' and procpar['seqfil'] == 'fse3d' and \
       fdf_properties['rank'] == 3:
        range_max = fdf_properties['matrix'][1]
        num_slicepts = fdf_properties['matrix'][
            0] * fdf_properties['matrix'][2]
        ds.Columns = fdf_properties['matrix'][2]
        ds.Rows = fdf_properties['matrix'][0]
        # FIXME FSE3d still producing bad dicoms

    if args.verbose:
        print "Columns ", ds.Columns, " Rows ", ds.Rows
        print "Range max and no slice points: ", range_max, num_slicepts
        print "Voldata[1] shape: ", voldata[:, :, 0].shape

    # Indexing in numpy matrix begins at 0, fdf/dicom filenames begin at 1
    for islice in xrange(0, range_max):
        # Reshape volume slice to 1D array
        slice_data = numpy.reshape(voldata[:, :, islice], (num_slicepts, 1))
        # Convert Pixel data to string
        ds.PixelData = slice_data.tostring()  # (7fe0,0010) Pixel Data

        # if acqndims == 3:
        if 'slice_no' in fdf_properties.keys():
            image_number = fdf_properties['slice_no']
        else:
            image_number = int(re.sub(r'^.*image(\d + ).*', r'\1', filename))

        new_filename = "slice%03dimage%03decho%03d.dcm" % (
            islice + 1, image_number, fdf_properties['echo_no'])

        if procpar['recon'] == 'external' and fdf_properties['rank'] == 3 and \
           procpar['seqfil'] == 'fse3d':
            pos = numpy.matrix([[0], [0], [islice], [1]])
        else:
            pos = numpy.matrix([[0], [0], [islice], [1]])

        Pxyz = transformation_matrix * pos
        ds.ImagePositionPatient = [
            str(Pxyz[0, 0]), str(Pxyz[1, 0]), str(Pxyz[2, 0])]

        # ds.FrameContentSequence[0].StackID = [str(volume)] # fourthdimid
        ds.FrameContentSequence[0].InStackPositionNumber = [
            int(islice)]  # fourthdimindex
        ds.FrameContentSequence[0].TemporalPositionIndex = ds.EchoNumber
        # ds.InStackPosition = islice #str(islice)

        # Save DICOM
        ds.save_as(os.path.join(outdir, new_filename))

        return ds


def Save2dFDFtoDicom(image_data, ds, fdf_properties, outdir, filename):
    """
    Export 2D image and metadata to DICOM

    """

    # Common export format
    ds.PixelData = image_data.tostring()         # (7fe0,0010) Pixel Data

    if len(ds.ImageType) >= 3 and ds.ImageType[2] == "ASL":
        image_number = fdf_properties['array_index']
        if fdf_properties["asltag"] == 1:               # Labelled
            new_filename = "slice%03dimage%03decho%03d.dcm" % (
                fdf_properties['slice_no'], image_number, 1)
        elif fdf_properties["asltag"] == -1:  # Control
            new_filename = "slice%03dimage%03decho%03d.dcm" % (
                fdf_properties['slice_no'], image_number, 2)
        else:                                            # Unknown
            new_filename = "slice%03dimage%03decho%03d.dcm" % (
                fdf_properties['slice_no'], image_number, 3)
        dcmfilename = os.path.join(outdir, new_filename)
    else:
        dcmfilename = os.path.join(
            outdir, os.path.splitext(filename)[0] + '.dcm')
    # Save DICOM
    ds.save_as(dcmfilename)

    return 1


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        usage=' ParseFDF.py -i "Input FDF directory"',
        description='''agilent2dicom is an FDF to Enhanced MR DICOM converter
        from MBI. ParseFDF takes header info from fdf files and adjusts the
        dicom dataset *ds* then rescales the image data.''')
    parser.add_argument(
        '-i', '--inputdir', help='''Input directory name. Must be an Agilent
        FDF image directory containing procpar and *.fdf files''',
        required=True)
    parser.add_argument(
        '-o', '--outputdir', help='Output directory name for DICOM files.')
    parser.add_argument(
        '-m', '--magnitude', help='Magnitude component flag.',
        action="store_true")
    parser.add_argument(
        '-p', '--phase', help='Phase component flag.', action="store_true")
    parser.add_argument(
        '-v', '--verbose', help='Verbose.', action="store_true")
    args = parser.parse_args()

    import ReadProcpar
    import RescaleFDF
    import ReadFDF as rf

    procpar, procpartext = ReadProcpar.ReadProcpar(
        os.path.join(args.inputdir, 'procpar'))
    ds, MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
    print "Rows: ", ds.Rows, " Columns: ", ds.Columns

    files = os.listdir(args.inputdir)
    fdffiles = [f for f in files if f.endswith('.fdf')]
    print "Number of FDF files ", len(fdffiles)
    rescaleintercept, rescaleslope = RescaleFDF.FindScale(
        fdffiles, ds, procpar, args)
    print "Rescale Intercept ", rescaleintercept, " Slope ", rescaleslope

    # for filename in fdffiles:
    filename = fdffiles[len(fdffiles) - 1]
    fdf_header, image_data = rf.ReadFDF(os.path.join(args.inputdir, filename))
    ds, fdfrank, matsize, tmatrix = ParseFDF(ds, fdf_header, procpar, args)
    ds, image_data = RescaleFDF.RescaleImage(
        ds, image_data, rescaleintercept, rescaleslope, args)

    print "FDF # of dims: ", fdfrank
    print "FDF matrix:  ", matsize
    print "Patient Name: ", ds.PatientName
    print "Patient ID: ", ds.PatientID
    print "Rows: ", ds.Rows
    print "Columns: ", ds.Columns
    print "Transformation matrix: ", tmatrix
    print "Stack ID: ", ds.FrameContentSequence[0].StackID
    print "InStack position: ", ds.FrameContentSequence[0].InStackPositionNumber
    print " Rescale slope: ", ds.RescaleSlope
    print " Resacle intcpt: ", ds.RescaleIntercept
    print " Image max : ", numpy.max([float("-inf"), image_data.max()])

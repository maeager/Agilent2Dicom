#!/usr/bin/env python

"""fid2dicom is used to convert Agilent FID files to DICOM format.

Version 0.1: Original code based on agilent2dicom FDF converter (Michael Eager)


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
import dicom
import numpy
import argparse

from fdf2dcm_global import *

import ReadProcpar
import RescaleFDF
import ReadFID
import ProcparToDicomMap
import ParseFDF




if __name__ == "__main__":

    # Parse command line arguments and validate img directory

    parser = argparse.ArgumentParser(usage=' fid2dicom -i "Input FID directory" [-o "Output directory"] [-m] [-p] [-v]',description='fid2dicom is an FID to Enhanced MR DICOM converter from MBI. Version '+VersionNumber)
    parser.add_argument('-i','--inputdir', help='Input directory name. Must be an Agilent FID image directory containing procpar and fid files',required=True);
    parser.add_argument('-o','--outputdir', help='Output directory name for DICOM files.');
    parser.add_argument('-m','--magnitude', help='Magnitude component flag.',action="store_true");
    parser.add_argument('-p','--phase', help='Phase component flag.',action="store_true");
    parser.add_argument('-s','--sequence', help='Sequence type (one of Multiecho, Diffusion, ASL.');
#    parser.add_argument('-d','--disable-dcmodify', help='Dcmodify flag.',action="store_true");
    parser.add_argument('-v','--verbose', help='Verbose.',action="store_true");
    
    # parser.add_argument("imgdir", help="Agilent .img directory containing procpar and fdf files")

    args = parser.parse_args()
    if args.verbose:
        print "Agilent2Dicom python converter FID to basic MR DICOM images\n Args: ", args #.accumulate(args.integers)


    # Check input folder exists
    if not os.path.exists(args.inputdir):
        print 'Error: Folder \'' + args.inputdir + '\' does not exist.'
        sys.exit(1)
        
    # Check folder contains procpar and *.fdf files
    files = os.listdir(args.inputdir)
    if args.verbose:
        print files
    if not args.sequence:
        args.sequence = ''

    if 'procpar' not in files:
        print 'Error: FID folder does not contain a procpar file'
        sys.exit(1)
        
    fidfiles = [ f for f in files if f.endswith('fid') ]
    if len(fdffiles) == 0:
        print 'Error: FID folder does not contain any fdf files'
        sys.exit(1)
    print "Number of FID files: ", len(fdffiles)
    # Check output directory
    if not args.outputdir:
        outdir = os.path.splitext(args.inputdir)[0]
        if not outdir.find('.img'):
            outdir = outdir + '.dcm'
        else:
            (dirName, imgdir) = os.path.split(outdir)
            while imdir == '':
                (dirName, imgdir) = os.path.split(outdir)

            (ImgBaseName, ImgExtension)=os.path.splitext(imgdir)
            outdir = os.path.join(dirName,ImgBaseName + '.dcm')
    else:
        outdir = args.outputdir
    if args.verbose:
        print 'Output directory: ', outdir

    if os.path.exists(outdir):
        if args.verbose:
            print 'Output folder ' + outdir + ' already exists'
        #sys.exit(1)
    else:
        if args.verbose:
            print 'Making output folder: ' + outdir    
        os.makedirs(outdir)
        
    

    ## Read in data procpar
    procpar, procpartext = ReadProcpar.ReadProcpar(os.path.join(args.inputdir,'procpar'))
    # if args.verbose:
    #     print procpar


    ## Map procpar to DICOM and create pydicom struct
    ds,MRAcquisitionType = ProcparToDicomMap.ProcparToDicomMap(procpar,args)


    # Calculate the max and min throughout all fdf iles in dataset;
    # calculate the intercept and slope for casting to UInt16
    RescaleIntercept,RescaleSlope = RescaleFDF.FindScale(fdffiles,ds,procpar,args)
 

    ## Per frame implementation
    # Read in data from fid file, if 3D split frames    
    volume=1

    filename = fidfiles[len(fidfiles)-1]
    procpar,hdr,dims,image_data_real,image_data_imag=ReadFID.readfid(args.inputdir,procpar)
    image_data,ksp=ReadFID.recon(procpar,dims,hdr,image_data_real,image_data_imag)

    
    if args.verbose:
        print 'Image_data shape:', str(image_data.shape)

    # Change dicom for specific FID header info
    ds,matsize,ImageTransformationMatrix = ReadFID.ParseFID(ds,fdf_properties,procpar,args)

    # Rescale image data
    ds,image_data =ReadFID.RescaleImage(ds,image_data,RescaleIntercept,RescaleSlope,args)

    # Export dicom to file
    if MRAcquisitionType == '3D':
        ds=ParseFDF.Save3dFDFtoDicom(ds,procpar,image_data,fdf_properties,ImageTransformationMatrix,args,outdir,filename)
    else:
        ParseFDF.Save2dFDFtoDicom(ds,image_data, outdir, filename)
        
    if ds.ImageType[2]=="MULTIECHO" or re.search('slab|img_',filename):
        print ds.FrameContentSequence[0].StackID, ds.FrameContentSequence[0].StackID[0]
        print type(ds.FrameContentSequence[0].StackID), type(ds.FrameContentSequence[0].StackID[0])
        ds.FrameContentSequence[0].StackID = str(int(ds.FrameContentSequence[0].StackID[0])+1)
        if args.verbose:
            print "Incrementing volume StackID ", ds.FrameContentSequence[0].StackID


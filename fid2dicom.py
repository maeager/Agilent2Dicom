#!/usr/bin/env python

"""fid2dicom is used to convert Agilent FID files to DICOM format.

  Version 0.1: Original code based on agilent2dicom FDF converter (Michael Eager)
  Version 0.2: Cplx filters upgraded and arguments improved for extra variables
  Version 0.3: Exporting DICOMs with correct parsing of Procpar and fid headers and rearrangement of image data

  $Id: fid2dicom.py,v d34a9504aac0 2015/01/23 03:40:06 michael $

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

from agilent2dicom_globalvars import *

import ReadProcpar
import ProcparToDicomMap
import ReadFID as FID
import cplxfilter as CPLX
import kspace_filter as KSP
import ParseFDF
import numpy,math,scipy
import nibabel as nib

def SaveNifti(image,basefilename):
    affine = numpy.eye(4)
    if image.ndim ==5:
        for i in xrange(0,image.shape[4]):
            raw_image=nib.Nifti1Image(numpy.abs(image[:,:,:,0,i]),affine)
            raw_image.set_data_dtype(numpy.float32)
            nib.save(raw_image,basefilename+'_0'+str(i)+'.nii.gz')
    else:
        raw_image=nib.Nifti1Image(numpy.abs(image),affine)
        raw_image.set_data_dtype(numpy.float32)
        nib.save(raw_image,basefilename+'.nii.gz')


def sigmaparse(s):
    try:
        print "Sigma parse ",s
        if s.find(','):
            print "Mapping "
            x, y, z = map(float, s.split(','))
            print x,y,z
            return (x, y, z)
        else:
            return (float(s),)
    except:
        raise argparse.ArgumentTypeError("Sigma must be single value or comma-seperated (eg. 1.0,1.0,3.0)")
        
if __name__ == "__main__":

    # Parse command line arguments and validate img directory

    parser = argparse.ArgumentParser(usage=' fid2dicom.py -i "Input FID directory" [-o "Output directory"] [-m] [-p] [-r] [-k] [-N] [-v] [[-g -s 1.0 [-go 0 -gm wrap]] [-l -s 1.0] [-d -n 5] [-w -n 5] [-y -s 1.0 -n 3]] [[-D] [-G -s 1.0] [-L -s 1.0][-Y -s 1.0 -n 3]]',description='fid2dicom is an FID to Enhanced MR DICOM converter from MBI. Complex filtering enabled with -g, -l, -n or -w arguments.  FID2DICOM Version '+VersionNumber)
    parser.add_argument('-i','--inputdir', help='Input directory name. Must be an Agilent FID image directory containing procpar and fid files',required=True);
    parser.add_argument('-o','--outputdir', help='Output directory name for DICOM files.');
    parser.add_argument('-m','--magnitude', help='Save Magnitude component. Default output of filtered image outputput',action="store_true");
    parser.add_argument('-p','--phase', help='Save Phase component.',action="store_true");
    parser.add_argument('-k','--kspace', help='Save Kspace data in outputdir-ksp.mat file',action="store_true");
    parser.add_argument('-r','--realimag', help='Save real and imaginary data in outputdir-real and outputdir-imag.',action="store_true");
    parser.add_argument('-N','--nifti',help='Save filtered outputs to NIFTI.',action="store_true");
    parser.add_argument('-D','--double-resolution',help='Zero pad k-space data before recnstruction to double resolution in image space.',action="store_true");
#    parser.add_argument('-s','--sequence', help='Sequence type (one of Multiecho, Diffusion, ASL).',choices={"MULTIECHO", "DIFFUSION", "ASL"});
#    parser.add_argument('-d','--disable-dcmodify', help='Dcmodify flag.',action="store_true");
    parser.add_argument('-g','--gaussian_filter', help='Gaussian filter smoothing of reconstructed RE and IM components.',action="store_true");
    parser.add_argument('-G','--FT_gaussian_filter', help='Gaussian filter smoothing in Fourier domain. Fourier transform of Gaussian filter applied to Kspace before reconstruction.',action="store_true");
    parser.add_argument('-l','--gaussian_laplace', help='Gaussian Laplace filter smoothing of reconstructed RE and IM components. -s/--sigma variable must be declared.',action="store_true");
    parser.add_argument('-s','--sigma',help='Gaussian and Laplace-Gaussian sigma variable. Default 1/srqt(2).',type=sigmaparse) # ,default='0.707',type=float
    parser.add_argument('-go','--gaussian_order',help='Gaussian and Laplace-Gaussian order variable. Default 0.',type=int,choices=[0,1,2,3],default=0)
    parser.add_argument('-gm','--gaussian_mode',help='Gaussian and Laplace-Gaussian mode variable. Default nearest.',choices=['reflect','constant','nearest','mirror','wrap'],default='nearest');
    parser.add_argument('-d','--median_filter', help='Median filter smoothing of reconstructed RE and IM components. ',action="store_true");
    parser.add_argument('-w','--wiener_filter', help='Wiener filter smoothing of reconstructed RE and IM components.',action="store_true");
    parser.add_argument('-y','--epanechnikov_filter', help='Epanechnikov filter smoothing of reconstructed RE and IM components.',action="store_true");
    parser.add_argument('-n','--window_size',type=int,help='Window size of Wiener and Median filters. Default 5.',default=5);
    parser.add_argument('-wn','--wiener_noise',help='Wiener filter noise. Estimated variance of image. If none or zero, local variance is calculated. Default 0=None.',default=0);
    parser.add_argument('-v','--verbose', help='Verbose.',action="store_true");
    
    # parser.add_argument("imgdir", help="Agilent .img directory containing procpar and fdf files");

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
#    if not args.sequence:
#        args.sequence = ''

    if 'procpar' not in files:
        print 'Error: FID folder does not contain a procpar file'
        sys.exit(1)
        
    fidfiles = [ f for f in files if f.endswith('fid') ]
    if len(fidfiles) == 0:
        print 'Error: FID folder does not contain any fid files'
        sys.exit(1)
    print "Number of FID files: ", len(fidfiles)
    ## Check output directory
    if not args.outputdir:
        outdir = os.path.splitext(args.inputdir)[0]
        if not outdir.find('.img'):
            outdir = outdir + '.dcm'
        else:
            (dirName, imgdir) = os.path.split(outdir)
            while imdir == '':
                (dirName, imgdir) = os.path.split(outdir)

            (ImgBaseName, ImgExtension)=os.path.splitext(imgdir)
            outdir = os.path.join(dirName,ImgBaseName+'.dcm')
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

 
    # Read in data from fid file
    volume=1

    filename = fidfiles[len(fidfiles)-1]
    procpar,hdr,dims,image_data_real,image_data_imag=FID.readfid(args.inputdir,procpar,args)

    # Change dicom for specific FID header info
    ds,matsize,ImageTransformationMatrix = FID.ParseFID(ds,hdr,procpar,args)

    image_data,ksp=FID.recon(procpar,dims,hdr,image_data_real,image_data_imag,args)

    
    if args.verbose:
        print 'Image_data shape:', str(image_data.shape)
        
    # Rescale image data
    #ds,image_scaled =FID.RescaleImage(ds,image_data,RescaleIntercept,RescaleSlope,args)

    FID.Save3dFIDtoDicom(ds,procpar,numpy.absolute(image_data),hdr,ImageTransformationMatrix,args,outdir)
    

    if args.gaussian_filter:
        if not args.sigma:
            args.sigma=0.707
        if not args.gaussian_order:
            args.gaussian_order=0
        if not args.gaussian_mode:
            args.gaussian_mode='nearest'
        if args.verbose:
            print "Computing Gaussian filtered image from Original image"
        image_filtered = CPLX.cplxgaussian_filter(image_data.real,image_data.imag,args.sigma,args.gaussian_order,args.gaussian_mode)
        ds.DerivationDescription='%s\nAgilent2Dicom Version: %s - %s\nScipy version: %s\nComplex Gaussian filter: sigma=%f order=%d mode=%s.' % (Derivation_Description,VersionNumber,DVCSstamp,scipy.__version__,args.sigma, args.gaussian_order,args.gaussian_mode)
        FID.SaveFIDtoDicom(ds,procpar,image_filtered,hdr,ImageTransformationMatrix,args,re.sub('.dcm','-gaussian.dcm',outdir))
        if args.nifti:
            SaveNifti(image_filtered,re.sub('.dcm','-gaussian',outdir))
            

    if args.gaussian_laplace:
        if not args.sigma:
            args.sigma=0.707
        if args.verbose:
            print "Computing Gaussian Laplace filtered image from Gaussian filtered image"
        image_filtered = CPLX.cplxgaussian_filter(image_data.real,image_data.imag,args.sigma)
        image_filtered = CPLX.cplxgaussian_laplace(image_filtered.real,image_filtered.imag,args.sigma)
        ds.DerivationDescription='%s\nAgilent2Dicom Version: %s - %s\nScipy version: %s\nComplex Gaussian filter: sigma=%f order=0 mode=nearest. Complex Laplace filter: sigma.' % (Derivation_Description,VersionNumber,DVCSstamp,scipy.__version__,args.sigma,args.sigma)
        FID.SaveFIDtoDicom(ds,procpar,image_filtered,hdr,ImageTransformationMatrix,args,re.sub('.dcm','-laplacegaussian.dcm',outdir))
        if args.nifti:
            SaveNifti(image_filtered,re.sub('.dcm','-laplacegauss',outdir))

    if args.median_filter:
        if not args.window_size:
            args.window_size=3
        if args.verbose:
            print "Computing Median filtered image from Original image"
        image_filtered = CPLX.cplxmedian_filter(image_data.real,image_data.imag,args.window_size)
        ds.DerivationDescription='%s\nAgilent2Dicom Version: %s - %s\nScipy version: %s\nComplex Median filter: windown size=%d.' % (Derivation_Description,VersionNumber,DVCSstamp,scipy.__version__,args.window_size)
        FID.SaveFIDtoDicom(ds,procpar,image_filtered,hdr,ImageTransformationMatrix,args,re.sub('.dcm','-median.dcm',outdir))
        if args.nifti:
            SaveNifti(image_filtered,re.sub('.dcm','-median',outdir))

         
    if args.wiener_filter:
        if not args.window_size:
            args.window_size=3
        if not args.wiener_noise:
            args.wiener_noise=None
        if args.verbose:
            print "Computing Wiener filtered image from Original image"
        image_filtered = CPLX.cplxwiener_filter(image_data.real,image_data.imag,args.window_size,args.wiener_noise)
        ds.DerivationDescription='%s\nAgilent2Dicom Version: %s - %s\nScipy version: %s\nComplex Wiener filter: window size=%d, noise=%f.' % (Derivation_Description,VersionNumber,DVCSstamp,scipy.__version__,args.window_size, args.wiener_noise)
        FID.SaveFIDtoDicom(ds,procpar,image_filtered,hdr,ImageTransformationMatrix,args,re.sub('.dcm','-wiener.dcm',outdir))
        if args.nifti:
            SaveNifti(image_filtered,re.sub('.dcm','-wiener',outdir))



    if args.epanechnikov_filter:
        if not args.sigma:
            args.sigma=np.sqrt(7.0/2.0)
        if not args.window_size:
            args.window_size=3
        if args.verbose:
            print "Computing Epanechnikov filtered image from Original image"
        image_filtered = CPLX.cplxepanechnikov_filter(image_data.real,image_data.imag,args.sigma,(3,3,3))
        ds.DerivationDescription='%s\nAgilent2Dicom Version: %s - %s\nScipy version: %s\nComplex Epanechnikov filter: sigma=%f window=%d order=0 mode=reflect.' % (Derivation_Description,VersionNumber,DVCSstamp,scipy.__version__,args.sigma,args.window_size)
        FID.SaveFIDtoDicom(ds,procpar,image_filtered,hdr,ImageTransformationMatrix,args,re.sub('.dcm','-epanechnikov.dcm',outdir))
        if args.nifti:
            SaveNifti(image_filtered,re.sub('.dcm','-epanechnikov',outdir))
            

    if args.FT_gaussian_filter:
        if not args.sigma:
            args.sigma=1.0/np.sqrt(2.0)
        print "Computing FT Gaussian filtered image"
        
        kspgauss1 = KSP.kspacegaussian_filter2(ksp1,sigma_=float(512.0)/float(1.0/np.sqrt(2.0)))
        from scipy.fftpack import ifftn,fftshift,ifftshift
        image_filtered = fftshift(ifftn(ifftshift(kspgauss1)));

        # image_filtered = CPLX.cplxgaussian_filter(image_data.real,image_data.imag,args.sigma,args.gaussian_order,args.gaussian_mode)
        ds.DerivationDescription='%s\nAgilent2Dicom Version: %s - %s\nScipy version: %s\nComplex Fourier Gaussian filter: sigma=%f.' % (Derivation_Description,VersionNumber,DVCSstamp,scipy.__version__,args.sigma)
        FID.SaveFIDtoDicom(ds,procpar,image_filtered,hdr,ImageTransformationMatrix,args,re.sub('.dcm','-kspgaussian.dcm',outdir))
        
        if args.nifti:
            SaveNifti(image_filtered,re.sub('.dcm','-kspgaussian',outdir))
              
    if args.FT_epanechnikov_filter:
        if not args.sigma:
            args.sigma=1.0/np.sqrt(2.0)
        if not args.window_size:
            args.window_size = 3
        print "Computing FT Epanechnikov filtered image"
        
        kspepa = KSP.kspaceepanechnikov_filter(ksp,sigma_=float(ksp.shape[0]])/float(1.0/np.sqrt(2.0)))
        from scipy.fftpack import ifftn,fftshift,ifftshift
        image_filtered = fftshift(ifftn(ifftshift(kspepa)));

        # image_filtered = CPLX.cplxgaussian_filter(image_data.real,image_data.imag,args.sigma,args.gaussian_order,args.gaussian_mode)
        ds.DerivationDescription='%s\nAgilent2Dicom Version: %s - %s\nScipy version: %s\nComplex Fourier Epanechnikov filter: sigma=%f.' % (Derivation_Description,VersionNumber,DVCSstamp,scipy.__version__,args.sigma)
        FID.SaveFIDtoDicom(ds,procpar,image_filtered,hdr,ImageTransformationMatrix,args,re.sub('.dcm','-kspepa.dcm',outdir))
        
        if args.nifti:
            SaveNifti(image_filtered,re.sub('.dcm','-kspepa',outdir))
              

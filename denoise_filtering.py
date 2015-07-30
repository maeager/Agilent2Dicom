#!/usr/bin/env python
""" Denoising techniques

Noiseest:  Noise estimation from two images
Phaserecon:  homodyne filter phase reconstruction

 - Michael Eager  (michael.eager@monash.edu)
"""
"""
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

import numpy as np
from scipy import ndimage

def tworepnoiseest(img1,img2):
    if np.iscomplexobj(img1) and np.iscomplexobj(img2):
        real_STD = noiseest(img1.real,img2.real)
        imag_STD = noiseest(img1.imag,img2.imag)
        return np.sqrt( (real_STD**2 + imag_STD**2) /2.0)
    else:
        # Normalise image
        nimg1 = (img1 - ndimage.minimum(img1)) / (ndimage.maximum(img1)- ndimage.minimum(img1)) 
        nimg2 = (img2 - ndimage.minimum(img2)) / (ndimage.maximum(img2) - ndimage.minimum(img2))

        #nimg1*=256.0
        #nimg2*=256.0
        
        return np.sqrt(0.5)*(nimg1-nimg2) 

def non_local_means_denoising(img,noise_estimate):

    if noise_estimate in None:
        #calculate noise estimate
        noise_estimate=0.05
    
    if np.iscomplexobj(img1):
        denoised_real,denoised_imag = NLmeans(img,noise_estimate)
        denoised_image = np.complex(denoised_real,denoised_imag)
    else:
        denoised = NLmeans(img,noise_estimate)


def NLmeansfilter(image,t,f,h):
    """
    input: image to be filtered
    t: radio of search window
    f: radio of similarity window
    h: degree of filtering
    
    Author: Jose Vicente Manjon Herrera & Antoni Buades
    Date: 09-03-2006
    
    Implementation of the Non local filter proposed for A. Buades, B. Coll and J.M. Morel in
    "A non-local algorithm for image denoising"
    
    Implemented in python, Michael Eager, Monash Biomedical Imaging 2015
    """

    # Size of the image
    sz=image.shape
 
    # Memory for the output
    output = np.empty_like(image)
    
    # Replicate the boundaries of the input image
    
    input2 = np.zeros(sz+f)
    input2[np.floor(f/2.0)-1:sz[0]+np.floor(f/2.0),np.floor(f/2.0)-1:sz[1]+np.floor(f/2.0)] = image
 
    # Used kernel
    kernel = make_kernel(f)
    kernel = kernel / np.sum(kernel)
 
    h=h*h
    j1=i1=0
    for i in xrange(0,sz[0]):
        for j in xrange(0,sz[1]):
            i1 = i+ f
            j1 = j+ f
            W1= input2[i1-f:i1+f , j1-f:j1+f]
            wmax=0
            average=0
            sweight=0
            rmin = np.max([i1-t,f+1]);
            rmax = np.min([i1+t,m+f]);
            smin = np.max([j1-t,f+1]);
            smax = np.min([j1+t,n+f]);
            
            for r in xrange(rmin,rmax+1):
                for s in xrange(smin,smax+1):
                    if (r==i1 && s==j1):
                        break
                        
                    W2 = input2(r-f:r+f , s-f:s+f);                
                 
                    d = np.sum(kernel*(W1-W2)*(W1-W2)))
                                               
                    w = np.exp(-d/h);                 
                                 
                    if w>wmax:                
                        wmax=w                   
                    sweight = sweight + w
                    average = average + w*input2[r,s]
                    
             
            average = average + wmax*input2[i1,j1]
            sweight = sweight + wmax
                   
            if sweight > 0:
                output[i,j] = average / sweight
            else:
                output[i,j] = image[i,j]

    return output

def make_kernel(fsize):
    """
    """
    kernel = np.zeros((2.0*fsize+1.0,2.0*fsize+1.0,2.0*fsize+1.0))
    for d in xrange(0,fsize):    
        value = 1.0 / (2.0*d+2)**2
        for ii in xrange(-d,d+1):
            for jj in xrange(-d,d+1):
                kernel[fsize-ii,fsize-jj] = kernel[fsize-ii,fsize-jj] + value
    kernel = kernel / fsize
    return kernel
        
def pipeline_one(magn_image):
    """
    automatic NL denoising

    """
    denoised_magn_image = non_local_means_denoising(magn_image)
    return denoised_magn_image

def pipeline_two(magn_image,magn_image_secondary):
    """
    Noise estimated NL denoising

    Estimate noise from two scans then implement NL means denoising
    """
    #calculate noise estimate
    noise_estimate = noiseest(magn_image,magn_image_secondary)
    denoised_magn_image = non_local_means_denoising(magn_image,noise_estimate)
    return denoised_magn_image


def pipeline_three(cmplx_ksp):
    """
    Complex NL filtering
    """

    if not np.iscomplexobj(cmplx_ksp):
        raise ValueError('pipeline three requires a complex k-space image')

    #Calculate homodyne filtered image
    mag_hdyne,pha_hdyne = phaserecon(cmplx_ksp,10.0)
    hdyne_image = mag_hdyne * np.exp(-1i* pha_hdyne)



if __name__ == "__main__":

    import os
    # import sys
    # import math
    # import re
    import argparse
    import ReadProcpar as Procpar
    import ProcparToDicomMap
    # import RescaleFDF
    import nibabel as nib
    from ReadFID import *

    parser = argparse.ArgumentParser(
        usage=' kspace_filters.py -i "Input FDF directory"',
        description='''kspace_filter algorithms for improving image
        quality.''')
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
        '-s', '--sequence', help='''Sequence type (one of Multiecho, Diffusion,
        ASL).''')
    parser.add_argument('-a', '--axis_order', help='Axis order eg 1,0,2.')
    parser.add_argument(
        '-v', '--verbose', help='Verbose comments.', action="store_true")
    args = parser.parse_args()

    # import ProcparToDicomMap as ptd

    procpar, procpartext = Procpar.ReadProcpar(
        os.path.join(args.inputdir, 'procpar'))
    ds, MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
    print "Rows: ", ds.Rows, " Columns: ", ds.Columns

    files = os.listdir(args.inputdir)
    fidfiles = [f for f in files if f.endswith('fid')]
    print "Number of FID files ", len(fidfiles)

    # for filename in fidfiles:
    print "Reading FID"
    filename = fidfiles[len(fidfiles) - 1]
    pp, hdr, dims, data_real, data_imag = readfid(args.inputdir,
                                                  procpar, args)
    print "Echoes: ", hdr['nEchoes'], " Channels: ", hdr['nChannels']
    affine = np.eye(4)
    # # affine[:3, :3]= np.arange(9).reshape((3, 3))
    # raw_data = nib.Nifti1Image(normalise(image_data_real), affine)
    # nib.save(raw_data, 'raw_data.nii.gz')

    print "Computing Original image (reconstruction)"
    image, ksp = recon(pp, dims, hdr,
                       data_real,
                       data_imag, args)
    del data_real, data_imag
    print "Shift kspace centre to max point"
    ksp = kspaceshift(ksp)

    if args.axis_order:
        image = RearrangeImage(image, args.axis_order, args)
        print "Transformed image shape: ", image.shape
        # np.delete(image)
        # image = imaget
    # print "Saving raw image"
    # save_nifti(normalise(np.abs(image)), 'raw_image')

    # print "Computing Gaussian filtered image from Original image"
    # image_filtered = simpleifft(kspacegaussian_filter(ksp,
    # 0.707, 0, 'nearest'))

    # print "Saving Gaussian image"
    # save_nifti(normalise(np.abs(image_filtered)), 'gauss_fourierimage')

    print "Computing Gaussian filtered2 image from Original image"
    kspgauss = kspacegaussian_filter2(ksp, 0.707)
    image_filtered = simpleifft(procpar, dims, hdr, kspgauss, args)
    # print "Saving Gaussian image"
    save_nifti(normalise(np.abs(image_filtered)), 'gauss_kspimage')

    print "Computing Homodyne filter from Original image"
    print "sigma ", 10.0
    pha,mag = phaserecon(ksp, 10.0,intpl=1.0,thr=0.05)
    
    # print "Saving Gaussian image"
    save_nifti(mag, 'homodyne_magnitude')
    save_nifti(pha,'homodyne_phase')

    print "Computing SWI images from Homodyne filtered images"
    swi_n,swi_p = swi2(mag,pha)

    if len(ksp.shape)==5:
        pha,mag1 = phaserecon(ksp[:,:,:,0,0], 10.0,intpl=1.0,thr=0.05)
        swi1 = swi(mag1,pha)
        pha,mag2 = phaserecon(ksp[:,:,:,0,1], 10.0,intpl=1.0,thr=0.05)
        swi1 = swi(mag2,pha)
        pha,mag3 = phaserecon(ksp[:,:,:,0,2], 10.0,intpl=1.0,thr=0.05)
        swi1 = swi(mag3,pha)

        mag_mee = multiecho_enhanced(mag1,mag2,mag3)
        swi_mee = multiecho_enhanced(swi1,swi2,swi3)

        save_nifti(mag, 'homodyne_magnitude')
        save_nifti(pha,'homodyne_phase')

        kspgauss = kspacegaussian_filter2(ksp[:,:,:,0,0], 0.707)
        image_filtered = simpleifft(procpar, dims, hdr, kspgauss, args)

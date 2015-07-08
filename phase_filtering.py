#!/usr/bin/env python
"""
  Phase enhanced image processing techniques:

  phaserecon()  Homodyne filter
  swi()         Susceptibility weighted filter
  swi2()        SWI with positive and negative outputs
  basicswi()    Initial testing of SWI
  Vesselenh()   Vessel enhancement with inversion


  Copyright (C) 2015 Michael Eager  (michael.eager@monash.edu)

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

import time
import numpy as np
import scipy
if scipy.__version__[2] == 7:
    # scipy.pkgload('signal')
    scipy.pkgload('ndimage')
    scipy.pkgload('fftpack')
else:
    from scipy.fftpack import fftn, ifftn, fftshift, ifftshift
    from scipy import ndimage
    # from scipy import signal

from kspace_filter import *



def phaserecon(ksp,sigma,intpl=1.0,thr=0.05):
    """Homodyne filter reconstruction

    Calculate phase images from complex inputs
    Inputs: 
            kimg -- reconstructed kspace signals
            [obsolete] kimgsos -- sos of square on complex images
            a    -- gaussian filter size 0 to 20, 10 is good balance
            intpl -- interplation factor
            thr -- thresholding: 0 no thresholding
                                 0.05 a good choice
    Outputs:
            pha -- phase images
            mag -- magnatitue images
            swi -- phase weighted magnatitue images
    
    Original MATLAB code: Zhaolin Chen @ Howard Florey Institute
    Adapted Python code: Michael Eager @ MBI
     modified LP filter 27-04-2015
     converted to python 01-05-2015
    """    
    ksize = np.array(ksp.shape[:3]) # [Nfe,Npe,Npe2] = size(kimg);

    # creat a Gaussian LPF
    (uu, vv, ww) = fouriercoords(ksize)
    # Original Gaussian kernel
    gfilter = gaussian_fourierkernel(uu, vv, ww, sigma)
 
    #Shift kspace data to centre
    ksp = kspaceshift(ksp)
 
    if intpl is None:
        intpl=1.0
    img = fftshift(fftn(ksp,intpl*ksize))
    img_lpf = fftshift(fftn(ksp*gfilter, intpl*ksize))

    img_hpf = img / img_lpf
    mag = np.abs(img)
    pha = np.angle(img_hpf);

    if thr is None:
        thr = 0.05

    #Calculate threshold
    thd = (thr/np.sqrt(np.sqrt(intpl)))*np.max(np.abs(img))
    # Get indexes of thresholded points
    i,j,k = np.where(abs(img) > thd)
    # Apply adjustment to values above threshold
    pha[i,j,k] = pha[i,j,k]/10

    return pha, mag
    

def swi(magn, phase,order):
    """ Suseptibility weighted image
    Magnitude and phase must be unwrapped data
    """
    if order is None:
        order=4.0

    weight = phase / np.pi + 1.0
    weight = weight.clip(min=0.0, max=1.0)
    return magn * (weight**order)



def swi2(magn, phase, order=4.0):
    """ Suseptibility weighted image - Second version
         Return negative and positive components of phase
    """
    if order is None:
        order = 4.0
    wnmask = (phase / np.pi + 1.0) 
    wnmask = wnmask.clip(min=0.0, max=1.0) 
    wpmask = (1 - phase / np.pi) 
    wpmask = wpmask.clip(min=0.0, max=1.0) 

    return magn * (wnmask**order), magn ** (wpmask**order)
         



def basicswi(cmplx_input_image, mask, order=2):
    """ Suseptibility-like weighted image -
    modified to use normalised local phase as the weighting
    """
    if np.iscomplexobj(cmplx_input_image):
        magn = np.abs(cmplx_input_image)
        phase = np.angle(cmplx_input_image)
        from scipy.ndimage.filters import uniform_filter
        normphase = uniform_filter(phase, 5.0, mode='constant', origin=-2.5)
        normphase = (normphase - ndimage.minimum(normphase)) / \
            (ndimage.maximum(normphase) - ndimage.minimum(normphase))
        weight = (normphase + 1.0)
        weight = weight.clip(min=0.0, max=1.0)
        return magn * (weight ** order) * mask
    else:
        print 'Error basicswi2: input image not complex'
        return np.abs(cmplx_input_image)


def multiecho_enhanced(img1,img2,img3):
    """
    Multi-echo enhancement - Inverse sum of squares
    I = (sum(img_i.^-p,i)./N).^(-1/p), where p = 2 or 3
    For 

    SNR of combined echo = sum(img_n.^-p)/ (sigma*sum(img_n.^(-2(p+1)))
    Phase variation deltaPhi= sum(img_n^2 * phi_i * TE_i)/ sum(img_n^2 * TE_i)

    Recommended contrast enhancement suggested by Zhaolin Chen.

    This algorithm favours the later echos, which contain more information from susceptible regions and finer vessels
    For Images without susceptibility, inverse is not needed, so let p=-3.

    Implementation: Michael Eager
    Algorithm: Zhaolin Chen
    Monash Biomedical Imaging
    """
    p = 3.0
    return ((img1**(-p) + img2**(-p) + img3**(-p))/3)**(-1.0/p)
        
def getbiggestcc(mask):
    """GETBIGGESTCC get the biggest connected component in binary image
    """
    # Find connected components
    label_im, nb_labels = ndimage.label(mask)
    sizes = ndimage.sum(mask, label_im, range(nb_labels + 1))
    max_pos= ndimage.maximum_position(sizes)
    masked_sizes = sizes > (ndimage.maximum(sizes) - 1) #==max_pos#
    return masked_sizes[label_im]

    
def vessel_enhance(img):
    """ Vessel enhancement algorithm
    """
    return (1 - normalise(img)**2)**2

def vessel_mask(img,thr=0.74,reps=2):

    # Apply threshiold to vessel_enhanced img
    mask = img < thr
    # Get biggest segment
    for i in xrange(0,reps):
        mask = ndimage.binary_erosion(mask, structure = np.ones((10, 10, 10)))
        mask = getbiggestcc(mask)
        mask = ndimage.binary_dilation(pixel, structure = np.ones((10, 10, 10)))
        mask = ndimage.binary_fill_holes(mask, structure = np.ones((20, 20, 20)))
    return img * mask

    

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

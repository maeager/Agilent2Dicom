#!/usr/bin/env python
"""
cplxgaussian_filter, cplxgaussian_laplace, cplxmedian_filter, cplxwiener_filter

Methods for complex filtering of 3D k-space data


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
import os,sys,math
import re
import struct
import argparse
import ProcparToDicomMap



from scipy.fftpack import fftn,ifftn,fftshift,ifftshift
from scipy import ndimage
from scipy import signal
import numpy as np
import ReadProcpar



"""
Examples 

https://scipy-lectures.github.io/advanced/image_processing/
https://scipy-lectures.github.io/packages/scikit-image/

using PIL
http://nbviewer.ipython.org/github/mroberts3000/GpuComputing/blob/master/IPython/GaussianBlur.ipynb


"""


"""
scipy.ndimage.filters.gaussian_filter

scipy.ndimage.filters.gaussian_filter(input, sigma, order=0, output=None, mode='reflect', cval=0.0, truncate=4.0)[source]
Multidimensional Gaussian filter.

Parameters:
input : array_like
Input array to filter.
sigma : scalar or sequence of scalars
Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.
order : {0, 1, 2, 3} or sequence from same set, optional
The order of the filter along each axis is given as a sequence of integers, or as a single number. An order of 0 corresponds to convolution with a Gaussian kernel. An order of 1, 2, or 3 corresponds to convolution with the first, second or third derivatives of a Gaussian. Higher order derivatives are not implemented
output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
cval : scalar, optional
Value to fill past edges of input if mode is 'constant'. Default is 0.0
truncate : float
Truncate the filter at this many standard deviations. Default is 4.0.
Returns:
gaussian_filter : ndarray
Returned array of same shape as input.
Notes

The multidimensional filter is implemented as a sequence of one-dimensional convolution filters. The intermediate arrays are stored in the same data type as the output. Therefore, for output types with a limited precision, the results may be imprecise because intermediate results may be stored with insufficient precision.
"""


def cplxgaussian_filter(real_input, imag_input,sigma=0.707,order_=0,mode_='reflect', cval_=0.0):
    """CPLXFILTER gaussian filter of complex 3D image
     ndimage.filters.gaussian_filter is used to smooth real and imag components
    
    filtered_magnitude = cplxfilter(realimg,imagimg)

    :param RE,IM: real and imag NDimages
    :param sigma: [optional] optimal sigma = 1/sqrt(2).  Standard deviation for
    Gaussian kernel. The standard deviations of the Gaussian filter
    are given for each axis as a sequence, or as a single number, in
    which case it is equal for all axes.
    :param order: [optional] ]{0, 1, 2, 3} or sequence from same set, optional.
    The order of the filter along each axis is given as a sequence of
    integers, or as a single number. An order of 0 corresponds to
    convolution with a Gaussian kernel. An order of 1, 2, or 3
    corresponds to convolution with the first, second or third
    derivatives of a Gaussian. Higher order derivatives are not
    implemented
    :param mode: [optional] ]{'reflect', 'constant', 'nearest', 'mirror', 'wrap'},
    optional. The mode parameter determines how the array borders are
    handled, where cval is the value when mode is equal to
    'constant'. Default is 'reflect'
    :param cval: [optional] Value to fill past edges of input if mode is 'constant'. Default is 0.0

    :return filtered_image:  complex 3D array of gaussian filtered image,

    """
    # imgmag = load_untouch_nii(file1)
    # imph = load_untouch_nii(file2)
    # imcplx=(imgmag.img).*exp((1i).*(imph.img))
    # c = fftshift(ifftn(fftshift(imcplx)))
    # c = flipdim(flipdim(flipdim(c,1),2),3)
    #f = fspecial3('gaussian', [3 3 3], 1) %/sqrt(2))
    # compleximg = complex(imfilter(real(c), f), imfilter(imag(c), f))

    print "Complex Gaussian filter"
    real_img = ndimage.filters.gaussian_filter(real_input, sigma, order=order_, mode=mode_, cval=cval_) #, truncate=4.0) truncate not supported in scipy 0.10
    imag_img = ndimage.filters.gaussian_filter(imag_input, sigma, order=order_, mode=mode_, cval=cval_) #, truncate=4.0)
    filtered_image = np.empty_like(real_input, dtype=complex)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplxfilter    



"""
scipy.ndimage.filters.gaussian_laplace

scipy.ndimage.filters.gaussian_laplace(input, sigma, output=None, mode='reflect', cval=0.0, **kwargs)[source]
Multidimensional Laplace filter using gaussian second derivatives.

Parameters:
input : array_like
Input array to filter.
sigma : scalar or sequence of scalars
The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.
output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
cval : scalar, optional
Value to fill past edges of input if mode is 'constant'. Default is 0.0
Extra keyword arguments will be passed to gaussian_filter().
"""


def cplxgaussian_laplace(real_input, imag_input,sigma=0.707,mode_='reflect', cval_=0.0):
    """CPLXFILTER gaussian laplace filter of complex 3D image
     ndimage.filters.gaussian_laplace is used to process real and imag components and estimate blobs
    
    filtered_magnitude = cplxfilter(realimg,imagimg)

    :param sigma:  optimal sigma = 1/sqrt(2).  Standard deviation for kernel. The standard deviations of the Gaussian Laplace filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.
 
    :param mode: {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional. The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
    :param cval: [optional] Value to fill past edges of input if mode is 'constant'. Default is 0.0

    :return filtered_image:  complex 3D array of gaussian laplace filtered image, 
    """

    print "Complex Gaussian_laplace filter"
    real_img = ndimage.filters.gaussian_laplace(real_input, sigma,  mode=mode_, cval=cval_) #, truncate=4.0) truncate not supported in scipy 0.10
    imag_img = ndimage.filters.gaussian_laplace(imag_input, sigma,  mode=mode_, cval=cval_) #, truncate=4.0)
    filtered_image = np.empty_like(real_input, dtype=complex)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplxgaussian_laplace    


"""
scipy.ndimage.filters.laplace

scipy.ndimage.filters.laplace(input, output=None, mode='reflect', cval=0.0)[source]
N-dimensional Laplace filter based on approximate second derivatives.

Parameters:
input : array_like
Input array to filter.
output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
cval : scalar, optional
Value to fill past edges of input if mode is 'constant'. Default is 0.0

"""


"""
scipy.ndimage.fourier.fourier_gaussian

scipy.ndimage.fourier.fourier_gaussian(input, sigma, n=-1, axis=-1, output=None)[source]
Multi-dimensional Gaussian fourier filter.

The array is multiplied with the fourier transform of a Gaussian kernel.

Parameters:
input : array_like
The input array.
sigma : float or sequence
The sigma of the Gaussian kernel. If a float, sigma is the same for all axes. If a sequence, sigma has to contain one value for each axis.
n : int, optional
If n is negative (default), then the input is assumed to be the result of a complex fft. If n is larger than or equal to zero, the input is assumed to be the result of a real fft, and n gives the length of the array before transformation along the real transform direction.
axis : int, optional
The axis of the real transform.
output : ndarray, optional
If given, the result of filtering the input is placed in this array. None is returned in this case.
Returns:
fourier_gaussian : ndarray or None
The filtered input. If output is given as a parameter, None is returned.

"""

"""
scipy.ndimage.filters.median_filter(input, size=None, footprint=None, output=None, mode='reflect', cval=0.0, origin=0)[source]
Calculates a multidimensional median filter.

Parameters:
input : array_like
Input array to filter.
size : scalar or tuple, optional
See footprint, below
footprint : array, optional
Either size or footprint must be defined. size gives the shape that is taken from the input array, at every element position, to define the input to the filter function. footprint is a boolean array that specifies (implicitly) a shape, but also which of the elements within this shape will get passed to the filter function. Thus size=(n,m) is equivalent to footprint=np.ones((n,m)). We adjust size to the number of dimensions of the input array, so that, if the input array is shape (10,10,10), and size is 2, then the actual size used is (2,2,2).
output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
cval : scalar, optional
Value to fill past edges of input if mode is 'constant'. Default is 0.0
origin : scalar, optional
The origin parameter controls the placement of the filter. Default 0.0.

Returns:
median_filter : ndarray
Return of same shape as input.
"""

def cplxmedian_filter(real_input,imag_input,size_=5):
    # scipy.ndimage.filters.median_filter(input, size=None, footprint=None, output=None, mode='reflect', cval=0.0, origin=0)
    # not used footprint_=[5,5,5], output_=None, mode_='reflect', cval_=0.0, origin_=0
    print "Complex Median filter"
    real_img = ndimage.filters.median_filter(real_input, size=size_) 
    imag_img = ndimage.filters.median_filter(imag_input, size=size_)
    filtered_image = np.empty_like(real_input, dtype=complex)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplxmedian_filter        




"""
scipy.signal.wiener(im, mysize=None, noise=None)[source]
Perform a Wiener filter on an N-dimensional array.

The Wiener filter is a simple deblurring filter for denoising images. This is not the Wiener filter commonly described in image reconstruction problems but instead it is a simple, local-mean filter.

Apply a Wiener filter to the N-dimensional array im.

Parameters:
im : ndarray
An N-dimensional array.
mysize : int or arraylike, optional
A scalar or an N-length list giving the size of the Wiener filter window in each dimension. Elements of mysize should be odd. If mysize is a scalar, then this scalar is used as the size in each dimension.
noise : float, optional
The noise-power to use. If None, then noise is estimated as the average of the local variance of the input.
Returns:
out : ndarray
Wiener filtered result with the same shape as im.
"""

def cplxwiener_filter(real_input,imag_input,mysize_=5, noise_=None):
    #scipy.signal.wiener(im, mysize=None, noise=None)

    print "Complex Wiener filter"
    real_img = signal.wiener(real_input, mysize=mysize_,noise=noise_)
    imag_img = signal.wiener(imag_input, mysize=mysize_,noise=noise_)
    filtered_image = np.empty_like(real_input, dtype=complex)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplxwiener_filter        

    
    
def normalise(data):
    max = ndimage.maximum(data)
    min = ndimage.minimum(data)
    print "Normalise max %f  min %f" % (max, min)
    return data #(data - min) * (max - min) 

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage=' ParseFDF.py -i "Input FDF directory"',description='agilent2dicom is an FDF to Enhanced MR DICOM converter from MBI. ParseFDF takes header info from fdf files and adjusts the dicom dataset *ds* then rescales the image data.')
    parser.add_argument('-i','--inputdir', help='Input directory name. Must be an Agilent FDF image directory containing procpar and *.fdf files',required=True);
    parser.add_argument('-o','--outputdir', help='Output directory name for DICOM files.')
    parser.add_argument('-m','--magnitude', help='Magnitude component flag.',action="store_true");
    parser.add_argument('-p','--phase', help='Phase component flag.',action="store_true");
    parser.add_argument('-s','--sequence', help='Sequence type (one of Multiecho, Diffusion, ASL.');
    
    parser.add_argument('-v','--verbose', help='Verbose.',action="store_true");
    args = parser.parse_args()
    
    import ReadProcpar
    import RescaleFDF
    import nibabel as nib
    from ReadFID import *
    # import ProcparToDicomMap as ptd


    procpar, procpartext = ReadProcpar.ReadProcpar(os.path.join(args.inputdir,'procpar'))
    ds,MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
    print "Rows: ", ds.Rows, " Columns: ", ds.Columns
    

    files = os.listdir(args.inputdir)
    fidfiles = [ f for f in files if f.endswith('fid') ]
    print "Number of FID files ", len(fidfiles)

    # for filename in fidfiles:
    print "Reading FID"
    filename = fidfiles[len(fidfiles)-1]
    pp,hdr,dims,image_data_real,image_data_imag=readfid(args.inputdir,procpar)

    affine = np.eye(4)
    # # affine[:3,:3]= np.arange(9).reshape((3,3))
    # raw_data=nib.Nifti1Image(normalise(image_data_real),affine)
    # nib.save(raw_data,'raw_data.nii.gz')

    print "Computing Original image (reconstruction)"
    image,ksp=recon(pp,dims,hdr,image_data_real,image_data_imag)
    raw_image=nib.Nifti1Image(normalise(np.abs(image)),affine)
    nib.save(raw_image,'raw_image.nii.gz')
    raw_ksp=nib.Nifti1Image(normalise(np.abs(ksp)),affine)
    nib.save(raw_ksp,'raw_ksp.nii.gz')

    print "Computing Gaussian filtered image from Original image"
    image_filtered = cplxgaussian_filter(image.real,image.imag)
    new_image = nib.Nifti1Image(normalise(np.abs(image_filtered)),affine)
    nib.save(new_image,'new_image.nii.gz')

    print "Computing Gaussian Laplace image from Smoothed image"
    Log_filtered = cplxgaussian_laplace(image_filtered.real,image_filtered.imag)
    Log_image = nib.Nifti1Image(normalise(np.abs(Log_filtered)),affine)
    nib.save(Log_image,'Log_image.nii.gz')

    print "Computing Median filtered image"
    median_filtered = cplxmedian_filter(image.real,image.imag)
    median_image = nib.Nifti1Image(normalise(np.abs(median_filtered)),affine)
    nib.save(median_image,'median_image.nii.gz')
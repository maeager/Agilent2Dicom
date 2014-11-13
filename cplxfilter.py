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


import numpy as np
import scipy
if scipy.__version__[2] == 7:
    scipy.pkgload('signal')
    scipy.pkgload('ndimage')
    scipy.pkgload('fftpack')    
else:
    from scipy.fftpack import fftn,ifftn,fftshift,ifftshift
    from scipy import ndimage
    from scipy import signal

#from scipy.fftpack import fftn,ifftn,fftshift,ifftshift
#from scipy import ndimage
#from scipy import signal


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


def cplxgaussian_filter(real_input, imag_input,sigma=0.707,order_=0,mode_='nearest', cval_=0.0):
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

    print "Complex Gaussian filter sigma ", sigma, " order ", order_, " mode ", mode_
    if real_input.ndim == 3:
        real_img = ndimage.filters.gaussian_filter(real_input, sigma, order=order_, mode=mode_, cval=cval_) #, truncate=4.0) truncate not supported in scipy 0.10
        imag_img = ndimage.filters.gaussian_filter(imag_input, sigma, order=order_, mode=mode_, cval=cval_) #, truncate=4.0)
    else:
        real_img = np.empty_like(real_input,dtype=numpy.float32)
        imag_img = np.empty_like(real_input,dtype=numpy.float32)
        for echo in xrange(0,real_input.shape[4]):
            for n in xrange(0,real_input.shape[3]):
                real_img[:,:,:,n,echo] = ndimage.filters.gaussian_filter(real_input[:,:,:,n,echo], sigma, order=order_, mode=mode_, cval=cval_)
                imag_img[:,:,:,n,echo] = ndimage.filters.gaussian_filter(imag_input[:,:,:,n,echo], sigma, order=order_, mode=mode_, cval=cval_)              
    filtered_image = np.empty(real_input.shape, dtype=numpy.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplxfilter    



def cplx2dgaussian_filter(real_input, imag_input,sigma=0.707,order_=0,mode_='nearest', cval_=0.0):
    """CPLX2DFILTER gaussian filter of complex 3D image
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
    print "Complex 2D Gaussian filter sigma ", sigma, " order ", order_, " mode ", mode_
    if real_input.ndim == 3:
        for islice in xrange(0,real_input.shape[2]):
            real_img[:,:,islice] = ndimage.filters.gaussian_filter(real_input[:,:,islice], sigma, order=order_, mode=mode_, cval=cval_) #, truncate=4.0) truncate not supported in scipy 0.10
            imag_img[:,:,islice] = ndimage.filters.gaussian_filter(imag_input[:,:,islice], sigma, order=order_, mode=mode_, cval=cval_) #, truncate=4.0)
    else:
        real_img = np.empty_like(real_input,dtype=numpy.float32)
        imag_img = np.empty_like(real_input,dtype=numpy.float32)
        for echo in xrange(0,real_input.shape[4]):
            for n in xrange(0,real_input.shape[3]):
                for islice in xrange(0,real_input.shape[2]):
                    real_img[:,:,islice,n,echo] = ndimage.filters.gaussian_filter(real_input[:,:,islice,n,echo], sigma, order=order_, mode=mode_, cval=cval_)
                    imag_img[:,:,islice,n,echo] = ndimage.filters.gaussian_filter(imag_input[:,:,islice,n,echo], sigma, order=order_, mode=mode_, cval=cval_)              
    filtered_image = np.empty(real_input.shape, dtype=complex)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplx2dfilter_gaussian    



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

    :param sigma:  optimal sigma = 1/sqrt(2).  Standard
    deviation for kernel. The standard deviations of the Gaussian Laplace filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.
 
    :param mode: {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional. The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
    :param cval: [optional] Value to fill past edges of input if mode is 'constant'. Default is 0.0

    :return filtered_image:  complex 3D array of gaussian laplace filtered image, 
    """

    print "Complex Gaussian_laplace filter sigma ",sigma, " mode ", mode_
    if real_input.ndim ==3:
        real_img = ndimage.filters.gaussian_laplace(real_input, sigma,  mode=mode_, cval=cval_)
        imag_img = ndimage.filters.gaussian_laplace(imag_input, sigma,  mode=mode_, cval=cval_)
    else:
        real_img = np.empty_like(real_input,dtype=numpy.float32)
        imag_img = np.empty_like(real_input,dtype=numpy.float32)
        for echo in xrange(0,real_input.shape[4]):
            for n in xrange(0,real_input.shape[3]):
                real_img[:,:,:,n,echo] = ndimage.filters.gaussian_laplace(real_input[:,:,:,n,echo], sigma, mode=mode_, cval=cval_)
                imag_img[:,:,:,n,echo] = ndimage.filters.gaussian_laplace(imag_input[:,:,:,n,echo], sigma, mode=mode_, cval=cval_)
    filtered_image = np.empty(real_input.shape, dtype=complex)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplxgaussian_laplace    



def cplxlaplacian_filter(real_input, imag_input,mode_='reflect', cval_=0.0):
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
    """CPLXFILTER laplace filter of complex 3D image

    
    filtered_image = cplxlaplacian_filter(realimg,imagimg)

    :param mode: {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional. The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
    :param cval: [optional] Value to fill past edges of input if mode is 'constant'. Default is 0.0

    :return filtered_image:  complex 3D array of gaussian laplace filtered image, 
    """

    print "Complex Gaussian_laplace filter  ", " mode ", mode_
    if real_input.ndim ==3:
        real_img = ndimage.filters.laplace(real_input,   mode=mode_, cval=cval_)
        imag_img = ndimage.filters.laplace(imag_input,   mode=mode_, cval=cval_)
    else:
        real_img = np.empty_like(real_input,dtype=numpy.float32)
        imag_img = np.empty_like(real_input,dtype=numpy.float32)
        for echo in xrange(0,real_input.shape[4]):
            for n in xrange(0,real_input.shape[3]):
                real_img[:,:,:,n,echo] = ndimage.filters.laplace(real_input[:,:,:,n,echo],  mode=mode_, cval=cval_)
                imag_img[:,:,:,n,echo] = ndimage.filters.laplace(imag_input[:,:,:,n,echo],  mode=mode_, cval=cval_)
    filtered_image = np.empty(real_input.shape, dtype=complex)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplxlaplace    



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

def cplxmedian_filter(real_input,imag_input,size_=5,mode_='reflect'):
    # scipy.ndimage.filters.median_filter(input, size=None, footprint=None, output=None, mode='reflect', cval=0.0, origin=0)
    # not used footprint_=[5,5,5], output_=None, mode_='reflect', cval_=0.0, origin_=0
    filtered_image = np.empty(real_input.shape,dtype=numpy.complex64)
    print "Complex Median filter window size(s)", size_
    if real_input.ndim == 3:
        if size_.ndim == 1:
            print "Complex Median filter"
            real_img = ndimage.filters.median_filter(real_input, size=size_,mode=mode_) 
            imag_img = ndimage.filters.median_filter(imag_input, size=size_,mode=mode_)
        else:
            real_img = ndimage.filters.median_filter(real_input, footprint=size_,mode=mode_) 
            imag_img = ndimage.filters.median_filter(imag_input, footprint=size_,mode=mode_)
    else:
        real_img = np.empty_like(real_input,dtype=numpy.float32)
        imag_img = np.empty_like(real_input,dtype=numpy.float32)
        for echo in xrange(0,real_input.shape[4]):
            for n in xrange(0,real_input.shape[3]):
                real_img[:,:,:,n,echo] = ndimage.filters.median_filter(real_input[:,:,:,n,echo], (size_,size_,size_),mode=mode_) 
                imag_img[:,:,:,n,echo] = ndimage.filters.median_filter(imag_input[:,:,:,n,echo], (size_,size_,size_),mode=mode_)
                                                                               
    filtered_image = np.empty(real_input.shape, dtype=numpy.complex64)

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

def cplxwiener_filter(real_input,imag_input,size_=5, noise_=None):
    #scipy.signal.wiener(im, mysize=None, noise=None)
    # ,(size_,size_,size_)
    print "Complex Wiener filter window size ",size_, " noise ", noise_
    if real_input.ndim == 3:
        real_img = signal.wiener(real_input,(size_,size_,size),noise=noise_) #, mysize=size_,noise=noise_)
        imag_img = signal.wiener(imag_input,(size_,size_,size),noise=noise_) #, mysize=size_,noise=noise_)
    else:
        real_img = np.empty_like(real_input)
        imag_img = np.empty_like(real_input)
        for echo in xrange(0,real_input.shape[4]):
            for n in xrange(0,real_input.shape[3]):
                real_img[:,:,:,n,echo] = signal.wiener(real_input[:,:,:,n,echo], (size_,size_,size_),noise=noise_)
                imag_img[:,:,:,n,echo] = signal.wiener(imag_input[:,:,:,n,echo], (size_,size_,size_),noise=noise_)

    filtered_image = np.empty(real_input.shape,dtype=numpy.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplxwiener_filter        



def normalise(data):
    max = ndimage.maximum(data)
    min = ndimage.minimum(data)
    print "Normalise max %f  min %f" % (max, min)
    #return as float32
    return data.astype(numpy.float32) #(data - min) * (max - min) 

def save_nifti(image,basename):
    import nibabel as nib
    affine = np.eye(4)
    if image.ndim ==5:
        for i in xrange(0,image.shape[4]):
            new_image = nib.Nifti1Image(normalise(np.abs(image[:,:,:,0,i])),affine)
            new_image.set_data_dtype(numpy.float32)
            nib.save(new_image,basename+'_0'+str(i)+'.nii.gz')
    else:
        new_image = nib.Nifti1Image(normalise(np.abs(image)),affine)
        new_image.set_data_dtype(numpy.float32)
        nib.save(new_image,basename+'.nii.gz')

    
if __name__ == "__main__":

    import os,sys,math
    import re
    import argparse
    import ReadProcpar as Procpar
    import ProcparToDicomMap
    import RescaleFDF
    from ReadFID import *
    import nibabel as nib

    parser = argparse.ArgumentParser(usage=' ParseFDF.py -i "Input FDF directory"',description='agilent2dicom is an FDF to Enhanced MR DICOM converter from MBI. ParseFDF takes header info from fdf files and adjusts the dicom dataset *ds* then rescales the image data.')
    parser.add_argument('-i','--inputdir', help='Input directory name. Must be an Agilent FDF image directory containing procpar and *.fdf files',required=True);
    parser.add_argument('-o','--outputdir', help='Output directory name for DICOM files.')
    parser.add_argument('-m','--magnitude', help='Magnitude component flag.',action="store_true");
    parser.add_argument('-p','--phase', help='Phase component flag.',action="store_true");
    parser.add_argument('-s','--sequence', help='Sequence type (one of Multiecho, Diffusion, ASL).');
    parser.add_argument('-a','--axis_order', help='Axis order eg 1,0,2.');
    parser.add_argument('-v','--verbose', help='Verbose comments.',action="store_true");
    args = parser.parse_args()
    
    # import ProcparToDicomMap as ptd


    procpar, procpartext = Procpar.ReadProcpar(os.path.join(args.inputdir,'procpar'))
    ds,MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
    print "Rows: ", ds.Rows, " Columns: ", ds.Columns
    

    files = os.listdir(args.inputdir)
    fidfiles = [ f for f in files if f.endswith('fid') ]
    print "Number of FID files ", len(fidfiles)

    # for filename in fidfiles:
    print "Reading FID"
    filename = fidfiles[len(fidfiles)-1]
    pp,hdr,dims,image_data_real,image_data_imag=readfid(args.inputdir,procpar,args)
    print "Echoes: ", hdr['nEchoes'], " Channels: ", hdr['nChannels']
    affine = np.eye(4)
    # # affine[:3,:3]= np.arange(9).reshape((3,3))
    # raw_data=nib.Nifti1Image(normalise(image_data_real),affine)
    # nib.save(raw_data,'raw_data.nii.gz')

    print "Computing Original image (reconstruction)"
    if os.path.exists('raw_image_00.nii.gz'):
        nii = nib.load('raw_image_00.nii.gz')
        image[:,:,:,0,0] = nii.get_data()
        nii = nib.load('raw_image_01.nii.gz')
        image[:,:,:,0,1] = nii.get_data()
        nii = nib.load('raw_image_02.nii.gz')
        image[:,:,:,0,2] = nii.get_data()
    else:
        image,ksp=recon(pp,dims,hdr,image_data_real,image_data_imag,args)

        if args.axis_order:
            image = RearrangeImage(image,args.axis_order,args)
            print "Transformed image shape: ", image.shape
            #np.delete(image)
            #image = imaget
        print "Saving raw image"
        save_nifti(np.abs(image),'raw_image')       


        
    print "Computing Gaussian filtered image from Original image"
    image_filtered = cplxgaussian_filter(image.real,image.imag,0.707,0,'nearest')
    print "Saving Gaussian image"
    save_nifti(normalise(np.abs(image_filtered)),'gauss_image')

    print "Computing Laplacian enhanced image from Original image"
    image_filtered = image-cplxlaplacian_filter(image.real,image.imag,0.707)
    print "Saving enhanced image g(x,y,z) = f(x,y,z) - Laplacian[f(x,y,z)]"
    save_nifti(normalise(np.abs(image_filtered)),'laplace_image')

    print "Computing Gaussian Laplace image from Smoothed image"
    Log_filtered = cplxgaussian_laplace(image_filtered.real,image_filtered.imag)
    save_nifti(normalise(np.abs(Log_filtered)),'Log_image')

    print "Computing Median filtered image"
    median_filtered = cplxmedian_filter(image.real,image.imag,3)
    print "Saving Median"
    save_nifti(normalise(np.abs(median_filtered)),'median_image')

    print "Computing Wiener filtered image"
    wiener_filtered = cplxwiener_filter(image.real,image.imag,3,0.0001)
    print "Saving Wiener image"
    save_nifti(normalise(np.abs(wiener_filtered)),'wiener_image')


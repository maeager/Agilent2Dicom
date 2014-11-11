#!/usr/bin/env python
"""
kspacegaussian_filter, kspacegaussian_laplace, kspacelaplacian_filter

Methods for fourier-domain filtering of 3D k-space data


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


def kspacegaussian_filter(ksp,sigma=0.707,n_=-1, axis_=-1):
    """KSPACEFILTER gaussian filter of complex 3D image
     
    filtered_magnitude = kspacefilter(realimg,imagimg)


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

    print "Complex Gaussian filter sigma ", sigma
    if ksp.ndim == 3:
        out_img = ndimage.fourier.fourier_gaussian(ksp, sigma, n=n_, axis=axis_)
    else:
        img = np.empty_like(ksp)
        
        for echo in xrange(0,ksp.shape[4]):
            for n in xrange(0,ksp.shape[3]):
                out_img[:,:,:,n,echo] = ndimage.fourier.fourier_gaussian(ksp, sigma, n=n_, axis=axis_)

    return out_img
# end kspacegaussian_filter    


def fourierlaplace(ksp_shape):
    """
    Laplacian operator in Fourier domain is very simple:
      D^2 g(x,y,z) => -(4pi^2) (u^2 + v^2 + w^2) G(u,v,w) 
    """
    siz=ksp_shape[:3]
    sz = (np.array(siz))/2
    x=np.array(range(-sz[0],sz[0]))
    y=np.array(range(-sz[1],sz[1]))
    z=np.array(range(-sz[2],sz[2]))
    mult_fact=np.ones((len(y),len(x),len(z)))
    (u,v,w)= (x[np.newaxis,:,np.newaxis]*mult_fact, \
        y[:,np.newaxis,np.newaxis]*mult_fact, \
    z[np.newaxis,np.newaxis,:]*mult_fact)
    return -(4*np.pi*np.pi)*(u*u+v*v+w*w)/(siz[0]*siz[1]*siz[2])

def fouriergauss(siz,sigma):
    """
    Gaussian operator in Fourier domain is another Gaussian :
      g(x,y,z)=(1/sqrt(2*pi).sigma).exp(-(x^2+y^2+z^2)/2sigma^2)
          => A.(exp(-2*pi*pi*(u^2 + v^2 + w^2)*(sigma^2)) 
    """
    sz = (np.array(siz))/2
    xx = np.array(range(-int(sz[0]),int(sz[0])))
    yy = np.array(range(-int(sz[1]),int(sz[1])))
    zz = np.array(range(-int(sz[2]),int(sz[2])))
    mult_fact = np.ones((len(yy),len(xx),len(zz)))
    uu = xx[np.newaxis,:,np.newaxis]*mult_fact
    vv = yy[:,np.newaxis,np.newaxis]*mult_fact
    ww = zz[np.newaxis,np.newaxis,:]*mult_fact
    if type(sigma) is float:
        return np.exp(-np.pi*np.pi*(uu*uu+vv*vv+ww*ww)*(2*sigma*sigma))  
    else:
        return np.exp(-np.pi*np.pi*((sigma[0]*sigma[0]*uu*uu)+(sigma[1]*sigma[1]*vv*vv)+(ww*ww*sigma[2]*sigma[2])))

def fwhm(sigma):
    """
    The full width at half maximum is therefore given by
    FWHM=2sqrt(2ln2)sigma approx 2.3548sigma.
    """
    return sigma*2*sqrt(2*np.log(2))


    
def inhomogeneousCorrection(ksp,siz,sigma):
    """
    Gaussian operator in Fourier domain is another Gaussian :
      g(x,y,z)=(A/sqrt(2*pi).sigma).exp(-(x^2+y^2+z^2)/2sigma^2)
          => A.(exp(-pi*(u^2 + v^2 + w^2)*(2.sigma^2)) 

    g(x,y,z)=exp(-(x^2+y^2+z^2)/2*sigma^2)
          => sqrt(pi*2*sigma^2)(exp(-pi^2*(u^2 + v^2 + w^2)*(2.sigma^2))

    Use large sigma for smoothing the MR image
    """
     
    sz = (np.array(siz))/2
    xx = np.array(range(-int(sz[0]),int(sz[0])))
    yy = np.array(range(-int(sz[1]),int(sz[1])))
    zz = np.array(range(-int(sz[2]),int(sz[2])))
    mult_fact = np.ones((len(yy),len(xx),len(zz)))
    uu = xx[np.newaxis,:,np.newaxis]*mult_fact
    vv = yy[:,np.newaxis,np.newaxis]*mult_fact
    ww = zz[np.newaxis,np.newaxis,:]*mult_fact
    del xx,yy,zz
    
#    arg   = -(xx.*xx + yy.*yy)/(2*sigma*sigma)
    # arg = -(uu*uu+vv*vv+ww*ww)/(2*sigma*sigma)
    # hG   = np.exp(arg)
    # #h(h<eps*max(h(:))) = 0;
    # del arg, xx,yy,zz,uu,vv,ww
    # sumh = sum(hG(:));
    # if sumh != 0:
    #     hG  = hG/sumh
    # HG=fftn(ifftshift(hG))
#    HGh=sqrt(pi*2*sigma*sigma)*
    HG = np.exp(-np.pi*np.pi*(uu*uu+vv*vv+ww*ww)*(2*sigma*sigma))
    del uu,vv,ww
    
    kspHG = ksp*HG
    del HG
    fsmoothed=abs(fftshift(ifftn(ifftshift(kspHG))))
    fsmoothed=fsmoothed/ndimage.mean(fsmoothed)
    return fsmoothed
# end inhomogeneousCorrection

    
def kspacegaussian_filter2(ksp,sigma_=None):
    siz=ksp.shape[0:3]
    if not sigma_:
        sigma = np.ones(3)*siz/(4*sqrt(2*np.log(2)))
    if type(sigma_) is float:
        sigma= np.ones(3)*sigma_
    else:
        sigm=sigma_.copy()
    Fgauss = fouriergauss(siz,1/sigma)
    print "Complex Laplace Gaussian filter sigma ", sigma
    if ksp.ndim == 3:
        out_img = ksp * Fgauss
    else:
        img = np.empty_like(ksp)
        
        for echo in xrange(0,ksp.shape[4]):
            for n in xrange(0,ksp.shape[3]):
                out_img[:,:,:,n,echo] = ksp[:,:,:,n,echo] * Fgauss
    return out_img
# end kspacegaussian_filter2    
                      
def kspacelaplacegaussian_filter(ksp,sigma_=None):

    siz=ksp.shape[0:3]
    if not sigma_:
        sigma = np.ones(3)*siz/(4*sqrt(2*np.log(2)))
    if type(sigma_) is float:
        sigma= np.ones(3)*sigma_
    else:
        sigma=sigma_.copy()
      
    Flaplace = fourierlaplace(siz)
    Fgauss = fouriergauss(siz,1/sigma)
    print "Complex Laplace Gaussian filter sigma ", sigma
    if ksp.ndim == 3:
        out_img = ksp * Fgauss * Flaplace
    else:
        img = np.empty_like(ksp)
        
        for echo in xrange(0,ksp.shape[4]):
            for n in xrange(0,ksp.shape[3]):
                #out_img[:,:,:,n,echo] = ndimage.filters.fourier_gaussian(ksp, sigma, n=n_, axis=axis_)
                out_img[:,:,:,n,echo] = ksp[:,:,:,n,echo] * Flaplace * Fgauss
    return out_img
# end kspacelaplacegaussian_filter    

def kspaceshift(ksp):
    kmax = np.array(ndimage.maximum_position(ksp))
    siz=np.array(ksp.shape[0:3])
    sub=int(siz.astype(float)/2.) - int(kmax)
    print "Shifting kspace ", sub
    for x in xrange(0,3):
        ksp=np.roll(ksp,sub[x],axis=x)
    return 
#end kspaceshift




def normalise(data):
    max = ndimage.maximum(data)
    min = ndimage.minimum(data)
    print "Normalise max %f  min %f" % (max, min)
    #return as float32
    return data.astype(numpy.float32) #(data - min) * (max - min) 

    
if __name__ == "__main__":

    import os,sys,math
    import re
    import argparse
    import ReadProcpar as Procpar
    import ProcparToDicomMap
    import RescaleFDF
    import nibabel as nib
    from ReadFID import *

    parser = argparse.ArgumentParser(usage=' kspace_filters.py -i "Input FDF directory"',description='kspace_filter algorithms for improving image quality.')
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
    image,ksp=recon(pp,dims,hdr,image_data_real,image_data_imag,args)

    if args.axis_order:
        image = RearrangeImage(image,args.axis_order,args)
        print "Transformed image shape: ", image.shape
        #np.delete(image)
        #image = imaget
    print "Saving raw image"
    # if image.ndim ==5:
    #     for i in xrange(0,image.shape[4]):
    #         raw_image=nib.Nifti1Image(normalise(np.abs(image[:,:,:,0,i])),affine)
    #         raw_image.set_data_dtype(numpy.float32)
    #         nib.save(raw_image,'raw_image_0'+str(i)+'.nii.gz')
    # else:
    #     raw_image=nib.Nifti1Image(normalise(np.abs(image)),affine)
    #     raw_image.set_data_dtype(numpy.float32)
    #     nib.save(raw_image,'raw_image.nii.gz')
       


        
#     print "Computing Gaussian filtered image from Original image"
#     image_filtered = fftshift(fftn(ifftshift(kspacegaussian_filter(ksp,0.707,0,'nearest'))))
#     print "Saving Gaussian image"
#     if image_filtered.ndim ==5:
#         for i in xrange(0,image_filtered.shape[4]):
#             new_image = nib.Nifti1Image(normalise(np.abs(image_filtered[:,:,:,0,i])),affine)
#             new_image.set_data_dtype(numpy.float32)
#             nib.save(new_image,'gauss_image_0'+str(i)+'.nii.gz')
#     else:
#         new_image = nib.Nifti1Image(normalise(np.abs(image_filtered)),affine)
#         new_image.set_data_dtype(numpy.float32)
#         nib.save(new_image,'gauss_image.nii.gz')


        
    print "Computing Gaussian filtered2 image from Original image"
    kspgauss =kspacegaussian_filter2(ksp,256/0.707)
    image_filtered = fftshift(fftn(ifftshift(kspgauss)))
    print "Saving Gaussian image"
    if image_filtered.ndim ==5:
        for i in xrange(0,image_filtered.shape[4]):
            new_image = nib.Nifti1Image(normalise(np.abs(image_filtered[:,:,:,0,i])),affine)
            new_image.set_data_dtype(numpy.float32)
            nib.save(new_image,'gauss_image2_0'+str(i)+'.nii.gz')
    else:
        new_image = nib.Nifti1Image(normalise(np.abs(image_filtered)),affine)
        new_image.set_data_dtype(numpy.float32)
        nib.save(new_image,'gauss_image2.nii.gz')

    # inhomogeneousCorrection
    image_corr = inhomogeneousCorrection(ksp,ksp.shape,3.0/60.0)
    # print "Saving Correction image"
    
    new_image = nib.Nifti1Image(normalise(np.abs(image_filtered/image_corr)),affine)
    new_image.set_data_dtype(numpy.float32)
    nib.save(new_image,'image_inhCorr3.nii.gz')


    print "Computing Laplacian enhanced image"
    laplacian = fftshift(fftn(ifftshift(kspgauss * fourierlaplace(ksp.shape))))
    alpha=ndimage.mean(np.abs(image_filtered))/ndimage.mean(np.abs(laplacian))
    image_lfiltered = image_filtered - alpha*laplacian
    print "Saving enhanced image g(x,y,z) = f(x,y,z) - Laplacian[f(x,y,z)]"
    if image_lfiltered.ndim ==5:
        for i in xrange(0,image_filtered.shape[4]):
            new_image = nib.Nifti1Image(normalise(np.abs(image_lfiltered[:,:,:,0,i])),affine)
            new_image.set_data_dtype(numpy.float32)
            nib.save(new_image,'laplacian_image_0'+str(i)+'.nii.gz')
    else:
        new_image = nib.Nifti1Image(normalise(np.abs(limage_filtered)),affine)
        new_image.set_data_dtype(numpy.float32)
        nib.save(new_image,'laplacian_enhanced.nii.gz')
    del image_filtered, image_lfiltered
        
    print "Computing Gaussian Laplace image from Smoothed image"
    image_filtered = fftshift(fftn(ifftshift(kspacelaplacegaussian_filter(ksp,256.0/0.707))))
#    Log_filtered = kspacegaussian_laplace(image_filtered.real,image_filtered.imag)
    Log_image = nib.Nifti1Image(normalise(np.abs(image_filtered)),affine)
    Log_image.set_data_dtype(numpy.float32)
    nib.save(Log_image,'Log_image.nii.gz')



    print "Double res and gaussian filter"
    ksplarge=np.zeros(np.array(ksp.shape)*2, dtype=complex )
    szmin = np.array(ksp.shape)/2 - 1
    szmax = np.array(ksp.shape) + szmin 
    ksplarge[szmin[0]:szmax[0],szmin[1]:szmax[1],szmin[2]:szmax[2]]=kspacegaussian_filter2(ksp,256/0.707)
    image_filtered = fftshift(fftn(ifftshift(ksplarge)))
    print "Saving Gaussian image"
    if image_filtered.ndim ==5:
        for i in xrange(0,image_filtered.shape[4]):
            new_image = nib.Nifti1Image(normalise(np.abs(image_filtered[:,:,:,0,i])),affine)
            new_image.set_data_dtype(numpy.float32)
            nib.save(new_image,'gauss_large_0'+str(i)+'.nii.gz')
    else:
        new_image = nib.Nifti1Image(normalise(np.abs(image_filtered)),affine)
        new_image.set_data_dtype(numpy.float32)
        nib.save(new_image,'gauss_large.nii.gz')




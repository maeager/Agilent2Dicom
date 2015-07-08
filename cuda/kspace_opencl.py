#!/usr/bin/env python
"""
   kspacegaussian_filter, kspacegaussian_laplace, kspacelaplacian_filter

   Methods for fourier-domain filtering of 3D k-space data

   OPENCL derivative

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

# import reikna.cluda as cl    
# from reikna.cluda import dtypes, any_api
# from reikna.fft import FFT
# from reikna.core import Annotation, Type, Transformation, Parameter
from pyfft.cl import Plan
import pyopencl as cl
import pyopencl.array as cl_array

    
from pyopencl.tools import clear_first_arg_caches,get_gl_sharing_context_properties

def clinit():
    """Initialize OpenCL with GL-CL interop.
    """
    clear_first_arg_caches()
    plats = cl.get_platforms()
    # handling OSX
    if sys.platform == "darwin":
        ctx = cl.Context(properties=get_gl_sharing_context_properties(),
                             devices=[])
    else:
        ctx = cl.Context(properties=[(cl.context_properties.PLATFORM, plats[0])])
    queue = cl.CommandQueue(ctx)
    return ctx, queue

def clfftn(data):
    """ OpenCL FFT 3D
    """
    clear_first_arg_caches()
    #ctx = cl.create_some_context(interactive=False)
    #queue = cl.CommandQueue(ctx)
    ctx,queue = clinit()
    plan = Plan(data.shape, normalize=True, queue=queue)
    # forward transform on device
    gpu_data = cl_array.to_device(queue, data)
    # forward transform
    plan.execute(gpu_data.data) 
    #result = gpu_data.get()
    result = gpu_data.get()
    return result
    
def clifftn(data):
    clear_first_arg_caches() 
    ctx = cl.create_some_context(interactive=False)
    queue = cl.CommandQueue(ctx)
    plan = Plan(data.shape, normalize=True, queue=queue)
    #Inverse transform:
    plan.execute(gpu_data.data, inverse=True) 
    result = gpu_data.get()
    return result
    
# sigma = np.ones(3)

def kspacegaussian_filter_pyfftCL(ksp,sigma):
    clear_first_arg_caches()
    sz = ksp.shape
    dtype = np.complex64
    ftype = np.float32
    #api = cluda.ocl_api()
    ctx = cl.create_some_context(interactive=False)
    queue = cl.CommandQueue(ctx)
    queue.flush()
    data_dev = cl_array.to_device(queue,ksp)
    w = h = k = 512
    plan = Plan((w,h,k), normalize=True, queue=queue)
    FACTOR=1.0
    program = cl.Program(ctx,"""
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "pyopencl-complex.h" 
__kernel void gauss_kernel(__global cfloat_t *dest) //, __global cfloat_t *src)
{
  uint x = get_global_id(0);uint y = get_global_id(1);uint z = get_global_id(2);
  uint dim1= %d;
  uint dim2= %d;
  uint dim3= %d;                    
  float sigma[3];
  sigma[0]=%f;sigma[1]=%f;sigma[2]=%f;
  float factor = %f;            
  float TWOPISQ = 19.739208802178716; //6.283185307179586;  //2*3.141592;

  ulong idx = z*dim1*dim2 + y*dim1 + x;
  float i = (float)(x);  //(x / dim3) / dim2);
      i = (i - (float)floor((float)(dim1)/2.0f))/(float)(dim1);
  float j = (float)y; //(x / dim3);
      //if((int)j > dim2) {j=(float)fmod(j, (float)dim2);};
      j = (j - (float)floor((float)(dim2)/2.0f))/(float)(dim2);
  //Account for large global index (stored as ulong) before performing modulus
  //double pre_k=fmod((double)(x) , (double) dim3);
  float k = (float) z; // pre_k;
      k = (k - (float)floor((float)(dim3)/2.0f))/(float)(dim3);

  float weight = exp(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  dest[idx].x = dest[idx].x * weight;
  dest[idx].y = dest[idx].y * weight; 
  
}
""" % (sz[0],sz[1],sz[2],sigma[0],sigma[1],sigma[2],FACTOR)).build()
    gauss_kernel = program.gauss_kernel
    #data_dev = thr.empty_like(ksp_dev)
    gauss_kernel(queue, sz, None, data_dev.data).wait() #, data_dev.data
    ksp_out = data_dev.get()
    queue.flush()
    ctx = cl.create_some_context(interactive=False)
    queue = cl.CommandQueue(ctx)
    w = h = k = 512
    plan = Plan((w,h,k), normalize=True, queue=queue)
    data2_dev = cl_array.to_device(queue, ksp_out)
    plan.execute(data2_dev.data, inverse=True)
    result = data2_dev.get()
    result = np.fft.fftshift(result)
    queue.finish()
    return result  #,ksp_out


    

    
def fouriercoords(siz):
    """fouriercoords
    Create x,y,z mesh of Fourier domain space
    """
    sz = np.ceil(np.array(siz) / 2.0)
    xx = np.array(range(-int(sz[0]), int(sz[0])))
    yy = np.array(range(-int(sz[1]), int(sz[1])))
    maxlen = ndimage.maximum(np.array(siz))
    if len(siz) == 3:
        zz = np.array(range(-int(sz[2]), int(sz[2])))
        mult_fact = np.ones((len(xx), len(yy), len(zz)))
        uu = xx[:, np.newaxis, np.newaxis] * mult_fact / maxlen  # * voxmm[0]
        vv = yy[np.newaxis, :, np.newaxis] * mult_fact / maxlen  # * voxmm[0]
        ww = zz[np.newaxis, np.newaxis, :] * mult_fact / maxlen  # * voxmm[0]
        if np.prod(siz) != np.prod(sz * 2):
            uu = uu[:siz[0], :siz[1], :siz[2]]
            vv = vv[:siz[0], :siz[1], :siz[2]]
            ww = ww[:siz[0], :siz[1], :siz[2]]
        return (uu, vv, ww)
    else:
        mult_fact = np.ones((len(xx), len(yy)))
        uu = xx[:, np.newaxis] * mult_fact / maxlen  # * voxmm[0]
        vv = yy[np.newaxis, :] * mult_fact / maxlen  # * voxmm[0]
        if np.prod(siz) != np.prod(sz * 2):
            uu = uu[:siz[0], :siz[1]]
            vv = vv[:siz[0], :siz[1]]
        return (uu, vv, [])


def kspacegaussian_filter(ksp,sigma):
    from reikna import cluda
    from reikna.cluda import functions, dtypes
    sz = ksp.shape
    dtype = np.complex64
    ftype = np.float32
    api = cluda.ocl_api()
    thr = api.Thread.create()
    FACTOR = 1.0
    program = thr.compile("""
KERNEL void gauss_kernel(
    GLOBAL_MEM ${ctype} *dest,
    GLOBAL_MEM ${ctype} *src)
{
  const ${ultype} x = (${ultype}) get_global_id(0);
  const SIZE_T dim1= %d;
  const SIZE_T dim2= %d;
  const SIZE_T dim3= %d;                    
  ${ftype} sigma[3];
  sigma[0]=%f;sigma[1]=%f;sigma[2]=%f;
  ${ftype} factor = %f;            
  const double TWOPISQ = 19.739208802178716; //6.283185307179586;  //2*3.141592;
  const ${ftype} SQRT2PI = 2.5066282746;
  const double CUBEDSQRT2PI = 15.749609945722419;
  const ${ultype} idx = x;
  ${ftype} i = (${ftype})((x / dim3) / dim2);
      i = (i - (${ftype})floor((${ftype})(dim1)/2.0))/(${ftype})(dim1);
  ${ftype} j = (${ftype})(x / dim3);
      if((SIZE_T)j > dim2) {j=(${ftype})fmod(j, (${ftype})dim2);};
      j = (j - (${ftype})floor((${ftype})(dim2)/2.0f))/(${ftype})(dim2);
  //Account for large global index (stored as ulong) before performing modulus
  double pre_k=fmod((double)(x) , (double) dim3);
  ${ftype} k = (${ftype}) pre_k;
      k = (k - (${ftype})floor((${ftype})(dim3)/2.0f))/(${ftype})(dim3);

  ${ftype} weight = exp(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  //${ftype} weight = expm1(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]))+1;
  //${ftype} weight= ${exp}(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  dest[idx].x = src[idx].x * weight;
  dest[idx].y = src[idx].y * weight; 
  
}
""" % (sz[0],sz[1],sz[2],sigma[0],sigma[1],sigma[2],FACTOR),
                          render_kwds=dict(ctype=dtypes.ctype(dtype),
                                           ftype=dtypes.ctype(ftype), ultype=dtypes.ctype(np.uint64),
                                           exp=functions.exp(ftype)),fast_math=True)

    gauss_kernel = program.gauss_kernel
    data_dev = thr.to_device(ksp)
    gauss_kernel(data_dev, data_dev, global_size=sz[0]*sz[1]*sz[2])
    ksp_out = data_dev.get()
    return ksp_out
        
def kspacegaussian_filter_CL(ksp,sigma):
    from reikna import cluda
    from reikna.cluda import functions, dtypes
    sz = np.array(ksp.shape)
    dtype = np.complex64
    ftype = np.float32
    api = cluda.ocl_api()
    thr = api.Thread.create()
    FACTOR = 1.0
    program = thr.compile("""
KERNEL void gauss_kernel(
    GLOBAL_MEM ${ctype} *dest,
    GLOBAL_MEM ${ctype} *src)
{
  const ${ultype} x = (${ultype}) get_global_id(0);
  const SIZE_T dim1= %d;
  const SIZE_T dim2= %d;
  const SIZE_T dim3= %d;                    
  ${ftype} sigma[3];
  sigma[0]=%f;sigma[1]=%f;sigma[2]=%f;
  ${ftype} factor = %f;            
  const double TWOPISQ = 19.739208802178716; //6.283185307179586;  //2*3.141592;
  const ${ftype} SQRT2PI = 2.5066282746;
  const double CUBEDSQRT2PI = 15.749609945722419;
  const ${ultype} idx = x;
  ${ftype} i = (${ftype})((x / dim3) / dim2);
      i = (i - (${ftype})floor((${ftype})(dim1)/2.0))/(${ftype})(dim1);
  ${ftype} j = (${ftype})(x / dim3);
      if((SIZE_T)j > dim2) {j=(${ftype})fmod(j, (${ftype})dim2);};
      j = (j - (${ftype})floor((${ftype})(dim2)/2.0f))/(${ftype})(dim2);
  //Account for large global index (stored as ulong) before performing modulus
  double pre_k=fmod((double)(x) , (double) dim3);
  ${ftype} k = (${ftype}) pre_k;
      k = (k - (${ftype})floor((${ftype})(dim3)/2.0f))/(${ftype})(dim3);

  ${ftype} weight = exp(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  //${ftype} weight = expm1(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]))+1;
  //${ftype} weight= ${exp}(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  dest[idx].x = src[idx].x * weight;
  dest[idx].y = src[idx].y * weight; 
  
}
""" % (sz[0],sz[1],sz[2],sigma[0],sigma[1],sigma[2],FACTOR),
                          render_kwds=dict(ctype=dtypes.ctype(dtype),
                                           ftype=dtypes.ctype(ftype), ultype=dtypes.ctype(np.uint64),
                                           exp=functions.exp(ftype)),fast_math=True)

    gauss_kernel = program.gauss_kernel
    data_dev = thr.empty_like(ksp)
    gauss_kernel(data_dev, data_dev, global_size=sz[0]*sz[1]*sz[2])
    ksp_out = data_dev.get()

    ifft = FFT(data_dev)
    cifft = ifft.compile(thr)
    cifft(data_dev, data_dev,inverse=0)
    result = np.fft.fftshift(data_dev.get() / sz[0]*sz[1]*sz[2])
    result = result[::-1,::-1,::-1]
    result = np.roll(np.roll(np.roll(result,1,axis=2),1,axis=1),1,axis=0)
    
    return ksp_out


def kspacegaussian_filter_CL2(ksp,sigma):
    """ Kspace gaussian filter and recon using GPU OpenCL

    1. GPU intialisation
    2. push KSP complex matrix to GPU
    3. declare FFT program
    4. declare Complex Gaussian GPU filter program
    5. Execute Gaussian GPU program
    6. GPU sync
    7. Execute FFT Recon
    8. Execute FFTshift
    9. Retrieve reconstruced complex image from GPU
    10. Reorganise image to standard (mimic numpy format)
    
    """
    sz=ksp.shape
    dtype = np.complex64
    ftype = np.float32
    ultype = np.uint64
    #api = cluda.ocl_api()
    api = any_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(ksp)
    ifft = FFT(data_dev)
    FACTOR=1.0
    program = thr.compile("""
KERNEL void gauss_kernel(
    GLOBAL_MEM ${ctype} *dest,
    GLOBAL_MEM ${ctype} *src)
{
  const ulong x = get_global_id(0);
  const SIZE_T dim1= %d;
  const SIZE_T dim2= %d;
  const SIZE_T dim3= %d;                    
  ${ftype} sigma[3];
  sigma[0]=%f;sigma[1]=%f;sigma[2]=%f;
  ${ftype} factor = %f;            
  const double TWOPISQ = 19.739208802178716; //6.283185307179586;  //2*3.141592;
  const ${ftype} SQRT2PI = 2.5066282746;
  const double CUBEDSQRT2PI = 15.749609945722419;
  const ulong idx = x;
  ${ftype} i = (${ftype})((x / dim3) / dim2);
      i = (i - (${ftype})floor((${ftype})(dim1)/2.0f))/(${ftype})(dim1);
  ${ftype} j = (${ftype})(x / dim3);
      if((SIZE_T)j > dim2) {j=(${ftype})fmod(j, (${ftype})dim2);};
      j = (j - (${ftype})floor((${ftype})(dim2)/2.0f))/(${ftype})(dim2);
  // Account for large global index (stored as ulong) before performing modulus
  double pre_k=fmod((double)(x), (double)dim3);
  ${ftype} k = (${ftype}) pre_k;
      k = (k - (${ftype})floor((${ftype})(dim3)/2.0f))/(${ftype})(dim3);

  ${ftype} weight = exp(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  // ${ftype} weight = expm1(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]))+1;
  // ${ftype} weight= ${exp}(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  dest[idx].x = src[idx].x * weight;
  dest[idx].y = src[idx].y * weight; 
  
}
""" % (sz[0],sz[1],sz[2],sigma[0],sigma[1],sigma[2],FACTOR),
                          render_kwds=dict(ctype=dtypes.ctype(dtype),
                                           ftype=dtypes.ctype(ftype),
                                           exp=functions.exp(ftype)),fast_math=True)
    gauss_kernel = program.gauss_kernel
    gauss_kernel(data_dev, data_dev, global_size=sz[0]*sz[1]*sz[2])
    thr.synchronize()
    ## Recon
    #data_dev = thr.to_device(ksp)
    ifftobj = FFT(data_dev)
    cifft = ifftobj.compile(thr)
    fftshiftobj = FFTShift(data_dev)
    cfftshift = fftshiftobj.compile(thr)
    cifft(data_dev, data_dev,inverse=0)
    thr.synchronize()
    cfftshift(data_dev,data_dev)
    thr.synchronize()
    result2 = data_dev.get() / np.prod(np.array(ksp.shape))
    result2 = result2[::-1,::-1,::-1]
    thr.release()
    return result2


    

def gaussian_fourierkernel(siz, sigma_):
    """
    Create Gaussian Fourier filter kernel with GPU
    """
    if not hasattr(sigma, "__len__"):  # type(sigma) is float:
        sigma = np.ones(3)*sigma_
    elif len(sigma) == 2:
        sigma[2] = 0.0
    

    sz=siz
    ctype = np.complex64
    ftype = np.float32
    #api = cluda.ocl_api()
    api = any_api()
    thr = api.Thread.create()
    base = np.ones(siz,ctype)
    data_dev = thr.to_device(base)
    FACTOR=1.0
    program = thr.compile("""
KERNEL void gauss_kernel(
    GLOBAL_MEM ${ctype} *dest,
    GLOBAL_MEM ${ctype} *src)
{
  const ulong x = get_global_id(0);
  const SIZE_T dim1= %d;
  const SIZE_T dim2= %d;
  const SIZE_T dim3= %d;                    
  ${ftype} sigma[3];
  sigma[0]=%f;sigma[1]=%f;sigma[2]=%f;
  ${ftype} factor = %f;            
  const double TWOPISQ = 19.739208802178716; //6.283185307179586;  //2*3.141592;
  const ${ftype} SQRT2PI = 2.5066282746;
  const double CUBEDSQRT2PI = 15.749609945722419;
  const ulong idx = x;
  ${ftype} i = (${ftype})((x / dim3) / dim2);
      i = (i - (${ftype})floor((${ftype})(dim1)/2.0))/(${ftype})(dim1);
  ${ftype} j = (${ftype})(x / dim3);
      if((SIZE_T)j > dim2) {j=(${ftype})fmod(j, (${ftype})dim2);};
      j = (j - (${ftype})floor((${ftype})(dim2)/2.0f))/(${ftype})(dim2);
  //Account for large global index (stored as ulong) before performing modulus
  double pre_k=fmod((double)(x) , (double) dim3);
  ${ftype} k = (${ftype}) pre_k;
      k = (k - (${ftype})floor((${ftype})(dim3)/2.0f))/(${ftype})(dim3);

  ${ftype} weight = exp(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  //${ftype} weight = expm1(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]))+1;
  //${ftype} weight= ${exp}(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  dest[idx].x = src[idx].x * weight;
  dest[idx].y = src[idx].y * weight; 
  
}
""" % (sz[0],sz[1],sz[2],sigma[0],sigma[1],sigma[2],FACTOR),
                          render_kwds=dict(ctype=dtypes.ctype(ctype),
                                           ftype=dtypes.ctype(ftype),
                                           exp=functions.exp(ftype)),fast_math=True)
    gauss_kernel = program.gauss_kernel
    #data_dev = thr.empty_like(ksp_dev)
    gauss_kernel(data_dev, data_dev, global_size=sz[0]*sz[1]*sz[2])
    gfilter = data_dev.get()
    thr.synchronize()
    thr.release()

    return gfilter
    

def fourierlaplace(siz):
    """
    Laplacian operator in Fourier domain is very simple:
      D^2 g(x,y,z) => -(4pi^2) (u^2 + v^2 + w^2) G(u,v,w)
    """
#    (uu, vv, ww) = fouriercoords(siz)
#    # / (siz[0] * siz[1] * siz[2])
#    return -(4 * np.pi * np.pi) * (uu * uu + vv * vv + ww * ww)

    
    dtype = np.complex64
    ftype = np.float32
    #api = cluda.ocl_api()
    ctx = cl.create_some_context(interactive=False)
    queue = cl.CommandQueue(ctx)
    queue.flush()
    laplace=np.ones(siz)
    laplace_dev = cl_array.to_device(queue,laplace)
    
    FACTOR =  (4 * np.pi * np.pi).astype(ftype)
    program = cl.Program(ctx,"""
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "pyopencl-complex.h" 
__kernel void gauss_kernel(__global cfloat_t *dest) //, __global cfloat_t *src)
{
  uint x = get_global_id(0);uint y = get_global_id(1);uint z = get_global_id(2);
  uint dim1= %d;
  uint dim2= %d;
  uint dim3= %d;                    
  float factor = %f;            
  float TWOPISQ = 19.739208802178716; //6.283185307179586;  //2*3.141592;

  ulong idx = z*dim1*dim2 + y*dim1 + x;
  float i = (float)(x);  //(x / dim3) / dim2);
      i = (i - (float)floor((float)(dim1)/2.0f))/(float)(dim1);
  float j = (float)y; //(x / dim3);
      //if((int)j > dim2) {j=(float)fmod(j, (float)dim2);};
      j = (j - (float)floor((float)(dim2)/2.0f))/(float)(dim2);
  //Account for large global index (stored as ulong) before performing modulus
  //double pre_k=fmod((double)(x) , (double) dim3);
  float k = (float) z; // pre_k;
      k = (k - (float)floor((float)(dim3)/2.0f))/(float)(dim3);

  float weight = -factor*((i*i) + (j*j) + (k*k));
  dest[idx].x = dest[idx].x * weight;
  dest[idx].y = dest[idx].y * weight; 
  
}
""" % (sz[0],sz[1],sz[2])).build()
    gauss_kernel = program.laplace_kernel
    #data_dev = thr.empty_like(ksp_dev)
    laplace_kernel(queue, sz, None, laplace_dev.data).wait() #, data_dev.data
    laplace = laplace_dev.get()
    return laplace

    

def fourierlaplaceinhom(siz, sigma):
    """
    Laplacian operator in Fourier domain is very simple:
      D^2 g(x,y,z) => -(4pi^2) (u^2 + v^2 + w^2) G(u,v,w)
    """
    (uu, vv, ww) = fouriercoords(siz)
    laplace = -(uu * uu / (sigma[0] * sigma[0]) + vv * vv /
                (sigma[1] * sigma[1]) + ww * ww / (sigma[2] * sigma[2]))
    laplace = (laplace - ndimage.minimum(laplace)) / \
              (ndimage.maximum(laplace) - ndimage.minimum(laplace))
    return laplace


def fourierepanechnikov(siz, sigma):
    """
    Epanechnikov kernel in Fourier domain is
     A.(1-|x|^2)  => (3/2*w^3)(sin(w) - w*cos(w)/2)
    """
    # (uu, vv, ww) = fouriercoords(siz)
    # uu = uu + np.spacing(1)
    # vv = vv + np.spacing(1)
    # ww = ww + np.spacing(1)

    # if not hasattr(sigma, "__len__"):
    # #if type(sigma) is float or type(sigma) is numpy.float64:
    #     return ((3.0*sigma/16.0)/(np.pi*(uu + vv +
    #     ww)/(sigma))**3)*(np.sin(2*np.pi*(uu + vv + ww)/(sigma)) - np.pi*(uu
    #     + vv + ww)/(sigma)*np.cos(2*np.pi*(uu + vv + ww)/(sigma))/2)
    # else:
    # return ((3.0/16.0)/(np.pi*((uu**3)/sigma[0]**4 + (vv**3)/sigma[1]**4 +
    # (ww**3)/sigma[2]**4)))*(np.sin(2*np.pi*(uu/sigma[0] + vv/sigma[1] +
    # ww/sigma[2])) - np.pi*(uu/sigma[0] + vv/sigma[1] +
    # ww/sigma[2])*np.cos(2*np.pi*(uu/sigma[0] + vv/sigma[1] + ww/sigma[2])))

    from cplxfilter import epanechnikov_kernel
    if not hasattr(sigma, "__len__"):
        Kepa = epanechnikov_kernel(
            (np.ceil(sigma) + 1, np.ceil(sigma) + 1, np.ceil(sigma) + 1),
            sigma)
    else:
        print (np.ceil(sigma[0]) + 1,
               np.ceil(sigma[1]) + 1, np.ceil(sigma[2]) + 1)
        print sigma
        Kepa = epanechnikov_kernel((np.ceil(sigma[0]) + 1, np.ceil(sigma[1]) + 1,
                                    np.ceil(sigma[2]) + 1), sigma)
    Kfilter = np.zeros(np.array(siz), dtype=np.float32)
    szmin = np.floor(
        np.array(siz) / 2.0 - np.floor(np.array(Kepa.shape) / 2.0) - 1)
    szmax = np.floor(szmin + np.array(Kepa.shape))
    print "Epa filter size ", siz, " image filter ", Kepa.shape, " szmin ", szmin, " szmax ", szmax
    Kfilter[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]] = Kepa
    #return np.abs(fftshift(clfftn(Kfilter)))
    
    api = any_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(ksp)
    ifft = FFT(data_dev)
    cifft = ifft.compile(thr)
    fftshift = FFTShift(data_dev)
    cfftshift = fftshift.compile(thr)
    cifft(data_dev, data_dev)
    thr.synchronize()
    cfftshift(data_dev,data_dev)
    thr.synchronize()
    result2 = data_dev.get() / np.prod(np.array(ksp.shape))
    #result2 = result2[::-1,::-1,::-1]
    thr.release()
    return np.abs(result2)




def fwhm(sigma):
    """
    The full width at half maximum is therefore given by
    FWHM=2sqrt(2ln2)sigma approx 2.3548sigma.
    """
    return sigma * 2 * np.sqrt(2 * np.log(2))


def inhomogeneouscorrection(ksp, siz, sigma):
    """
    Gaussian operator in Fourier domain is another Gaussian :
      g(x,y,z)=(A/sqrt(2*pi).sigma).exp(-(x^2+y^2+z^2)/2sigma^2)
          => A.(exp(-pi*(u^2 + v^2 + w^2)*(2.sigma^2))

    g(x, y, z)=exp(-(x^2 + y^2 + z^2)/2*sigma^2)
          => sqrt(pi*2*sigma^2)(exp(-pi^2*(u^2 + v^2 + w^2)*(2.sigma^2))

    Use large sigma for smoothing the MR image
    """

    sz = np.ceil((np.array(siz)) / 2.0)
    xx = np.array(range(-int(sz[0]), int(sz[0])))
    yy = np.array(range(-int(sz[1]), int(sz[1])))
    zz = np.array(range(-int(sz[2]), int(sz[2])))
    mult_fact = np.ones((len(yy), len(xx), len(zz)))
    uu = xx[np.newaxis, :, np.newaxis] * mult_fact
    vv = yy[:, np.newaxis, np.newaxis] * mult_fact
    ww = zz[np.newaxis, np.newaxis, :] * mult_fact
    del xx, yy, zz, mult_fact
    if np.prod(siz) != np.prod(sz * 2):
        uu = uu[:siz[0], :siz[1], :siz[2]]
        vv = vv[:siz[0], :siz[1], :siz[2]]
        ww = ww[:siz[0], :siz[1], :siz[2]]

#    arg   = -(xx.*xx + yy.*yy)/(2*sigma*sigma)
    # arg = -(uu*uu + vv*vv + ww*ww)/(2*sigma*sigma)
    # hG   = np.exp(arg)
    # #h(h < eps*max(h(:))) = 0;
    # del arg, xx, yy, zz, uu, vv, ww
    # sumh = sum(hG(:));
    # if sumh != 0:
    #     hG  = hG/sumh
    # HG = fftn(ifftshift(hG))
#    HGh = sqrt(pi*2*sigma*sigma)*
    HG = np.expm1(-np.pi * np.pi * (uu * uu + vv * vv + ww * ww)
                * (2 * sigma * sigma))+1
    del uu, vv, ww

    kspHG = ksp * HG
    del HG
    fsmoothed = abs(fftshift(clifftn(ifftshift(kspHG))))
    fsmoothed = fsmoothed / ndimage.mean(fsmoothed)
    return fsmoothed
# end inhomogeneouscorrection


def kspaceepanechnikov_filter(ksp, bdwidth_=None):
    """
    Apply Epanechnikov filter in Fourier domain to kspace data
    """
    siz = ksp.shape[0:3]
    bdwidth = 0
    if 'bdwidth_' not in locals():
        bdwidth = np.sqrt(7) * np.array(siz) / (4 * np.sqrt(2 * np.log(2)))
    else:
        if not hasattr(bdwidth_, "__len__"):
            bdwidth = np.ones(3) * bdwidth_
        else:
            bdwidth = bdwidth_.copy()
    Fepanechnikov = fourierepanechnikov(siz, bdwidth)
    out_ksp = np.empty_like(ksp, dtype=np.complex64)
    print "Complex Epanechnikov filter bandwidth ", bdwidth
    if ksp.ndim == 3:
        out_ksp.real = Fepanechnikov * ksp.real
        out_ksp.imag = Fepanechnikov * ksp.imag
    else:
        for echo in xrange(0, ksp.shape[4]):
            for n in xrange(0, ksp.shape[3]):
                out_ksp[:, :, :, n, echo].real = ksp[
                    :, :, :, n, echo].real * Fepanechnikov
                out_ksp[:, :, :, n, echo].imag = ksp[
                    :, :, :, n, echo].imag * Fepanechnikov

    return out_ksp
# end kspaceepanechnikov_filter


from cplxfilter import epanechnikov_kernel

def kspaceepanechnikov_filter_CL2(ksp,sigma):
    sz = ksp.shape
    dtype = np.complex64
    ftype = np.float32
    clear_first_arg_caches()
    fsiz = (5,5,5)
    print (np.ceil(sigma[0]) + 2,
           np.ceil(sigma[1]) + 2, np.ceil(sigma[2]) + 2)
    print sigma
    fsiz =  (np.ceil(sigma)+2).astype(int)
    for i in xrange(0,fsiz.size):
        if not fsiz[i] & 0x1:
            fsiz[i]+=1
    ## Create image-domain Epanechikov kernel            
    Kepa = epanechnikov_kernel(fsiz, sigma)
    ## Place kernel at centre of ksp-sized matrix
    Kfilter = np.zeros(np.array(sz), dtype=np.complex64)
    szmin = np.floor(
        np.array(sz) / 2.0 - np.floor(np.array(Kepa.shape) / 2.0) - 1)
    szmax = np.floor(szmin + np.array(Kepa.shape))
    print "Epa filter size ", sz, " image filter ", Kepa.shape, " szmin ", szmin, " szmax ", szmax
    Kfilter[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]] = Kepa
    Kfilter[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]].imag = Kepa
    ## Create fourier-domain Epanechnikov filter
    api = any_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(Kfilter)
    rfft = FFT(data_dev)
    crfft = rfft.compile(thr)
    fftshift = FFTShift(data_dev)
    cfftshift = fftshift.compile(thr)
    crfft(data_dev, data_dev)
    thr.synchronize()
    cfftshift(data_dev,data_dev)
    Fepanechnikov =  np.abs(data_dev.get()) # / np.prod(np.array(ksp.shape))
    #result2 = result2[::-1,::-1,::-1]
    thr.synchronize()
    #result = np.zeros(np.array(siz), dtype=np.complex64)
    #result.real = np.abs(result2) / np.sqrt(2)
    #result.imag = np.abs(result2) / np.sqrt(2)
    del data_dev, rfft, crfft, fftshift, cfftshift
    # Multiply Epanechnikov filter to real and imag ksp data
    program = thr.compile("""
KERNEL void multiply_them(
    GLOBAL_MEM ${ctype} *dest,
    GLOBAL_MEM ${ctype} *a,
    GLOBAL_MEM ${ftype} *f)
{
  const SIZE_T i = get_local_id(0);
  dest[i].x = a[i].x * f[i];
  dest[i].y = a[i].y * f[i];
}""", render_kwds=dict(ctype=dtypes.ctype(dtype),ftype=dtypes.ctype(ftype)))
    
    data_dev = thr.to_device(ksp)
    filter_dev = thr.to_device(Fepanechnikov)
    multiply_them = program.multiply_them
    multiply_them(data_dev,data_dev,filter_dev, global_size=512*512*512)
    thr.synchronize()
    del filter_dev, program
    FACTOR=1.0

    ## Recon
    ifft = FFT(data_dev)
    cifft = ifft.compile(thr)
    fftshiftobj = FFTShift(data_dev)
    cfftshift = fftshiftobj.compile(thr)
    cifft(data_dev, data_dev,inverse=0)
    thr.synchronize()
    cfftshift(data_dev,data_dev)
    thr.synchronize()
    result2 = data_dev.get() / np.prod(np.array(ksp.shape))
    result2 = result2[::-1,::-1,::-1]
    thr.release()
    return result2

    


def kspaceshift(ksp):
    """kspaceshift
    Shift k-space data to centre maximum
    """
    print "K-space shift ", ksp.shape
    if len(ksp.shape) == 3:
        kmax = np.array(ndimage.maximum_position(np.abs(ksp)))
        siz = np.array(ksp.shape[0:3])
        sub = (siz / 2.).astype(int) - kmax
        print "Shifting kspace ", sub
        for x in xrange(0, 3):
            if sub[x] != 0:
                ksp = np.roll(ksp, sub[x], axis=x)
            print ""
    else:
        kmax = np.array(
            ndimage.maximum_position(np.squeeze(np.abs(ksp[:, :, :, 0, 0]))))
        siz = np.array(ksp.shape[0:3])
        sub = (siz / 2.).astype(int) - kmax
        for echo in xrange(0, ksp.shape[4]):
            for nchannel in xrange(0, ksp.shape[3]):
                print "Shifting kspace ", sub
                for x in xrange(0, 3):
                    if sub[x] != 0:
                        ksp[:, :, :, nchannel, echo] = np.roll(
                            ksp[:, :, :, nchannel, echo], sub[x], axis=x)
                    # print ""
    return ksp
# end kspaceshift


def imageshift(image1, image2):
    """imageshift
    Shift image2 into same arrangement as image 1 using the maximum
    value position. Return the shifted image2
    """
    print "Image shift ", image.shape, image2.shape
    i1max = np.array(ndimage.maximum_position(np.abs(image1)))
    i2max = np.array(ndimage.maximum_position(np.abs(image2)))
    siz = np.array(image1.shape[0:3])
    sub = i1max - i2max
    print "Shifting image ", sub, "size ", siz
    for x in xrange(0, 3):
        image2 = np.roll(image2, sub[x], axis=x)
    print ""
    return image2
# end imageshift


def open_image(image_filtered):
    """open_image example ndimage.grey_opening
    """
    c = ndimage.grey_opening(np.abs(image_filtered), size=(5, 5, 5))
    new_image = nib.Nifti1Image(normalise(c), affine)
    new_image.set_data_dtype(np.float32)
    nib.save(new_image, 'image_open.nii.gz')


def close_image(image_filtered):
    """close_image example ndimage.grey_closing
    """
    closed = ndimage.grey_closing(np.abs(image_filtered), size=(5, 5, 5))
    new_image = nib.Nifti1Image(normalise(closed), affine)
    new_image.set_data_dtype(np.float32)
    nib.save(new_image, 'image_fill.nii.gz')


def sobel_image(image):
    """sobel_image example ndimage.filters.sobel
    """
    d = ndimage.filters.sobel(image, axis=0)
    e = ndimage.filters.sobel(image, axis=1)
    f = ndimage.filters.sobel(image, axis=2)
    new_image = nib.Nifti1Image(np.abs(d) + np.abs(e) + np.abs(f), affine)
    new_image.set_data_dtype(np.float32)
    nib.save(new_image, 'sobel.nii.gz')


def normalise(data):
    """Normalise ndimage
    (must not be complex)
    """
    _max = ndimage.maximum(data)
    _min = ndimage.minimum(data)
    print "Normalise max %f  min %f" % (_max, _min)
    # return as float32
    data = ((data - _min) * (_max - _min))
    return data.astype(np.float32)


def save_nifti(image, basename):
    """save_nifti
    Save image as NIFTI
    """
    import nibabel as nib
    affine = np.eye(4)
    if image.ndim == 5:
        for echo in xrange(0, image.shape[4]):
            for channel in xrange(0, image.shape[3]):
                new_image = nib.Nifti1Image(
                    np.abs(image[:, :, :, channel, echo]), affine)
                new_image.set_data_dtype(np.float32)
                nib.save(new_image, basename + '_' + str(channel) + str(echo)
                         + '.nii.gz')
    else:
        new_image = nib.Nifti1Image(np.abs(image), affine)
        new_image.set_data_dtype(np.float32)
        nib.save(new_image, basename + '.nii.gz')


def save_nifti_int(image, basename):
    """save_nifti_int
    Save image as NIFTI in int format
    """
    import nibabel as nib
    affine = np.eye(4)
    if image.ndim == 5:
        for echo in xrange(0, image.shape[4]):
            for channel in xrange(0, image.shape[3]):
                new_image = nib.Nifti1Image(
                    np.abs(image[:, :, :, channel, echo]).astype(int), affine)
                new_image.set_data_dtype(np.int32)
                nib.save(new_image, basename + '_' + str(channel) + str(echo)
                         + '.nii.gz')
    else:
        new_image = nib.Nifti1Image(np.abs(image).astype(int), affine)
        new_image.set_data_dtype(np.int32)
        nib.save(new_image, basename + '.nii.gz')


def test_double_resolution(ksp, basename):
    """test_double_resolution
    """
    print "Double res and " + basename + " filter"
    # two 32-bit float
    ksplarge = np.zeros(np.array(ksp.shape) * 2, dtype=np.complex64)
    szmin = np.array(ksp.shape) / 2 - 1
    szmax = np.array(ksp.shape) + szmin
    ksplarge[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]] = ksp
    image_filtered = fftshift(clifftn(ifftshift(ksplarge)))
    print "Saving Double res image: " + basename + " filtered"
    save_nifti(np.abs(image_filtered), basename + '_large')


def test_depth_algorithm(image_filtered, basename='depth'):
    """test_depth_algorithm
    Depth algorithm testing
    """
    print "Testing depth algorithm"
    # t1 = time.time()

    # Close gaussian filtered image
    c = ndimage.grey_closing(np.abs(image_filtered), size=(5, 5, 5))
    # Mask closed image
    # cm = c * (c>8000).astype(float)
    cm = c / (ndimage.maximum(c))
    # avoid div by zero
    cm = 0.99 * cm + 0.00001
    # Regularise gaussian filtered image
    gm = (np.abs(image_filtered) / cm)  # * (c>8000).astype(float)
    # Depth = difference between closed image and regularised gaussian
    depth = c / ndimage.maximum(c) - gm / ndimage.maximum(gm)
    # mask regularised image
    depth = depth * (c > 0.00015).astype(float)
    # Normalise
    # depth = (depth -
    # ndimage.minimum(depth))/(ndimage.maximum(depth)-ndimage.minimum(depth))
    # save to nifti
    new_image = nib.Nifti1Image(np.abs(depth), affine)
    new_image.set_data_dtype(np.float32)
    nib.save(new_image, basename + '.nii.gz')


def test_double_resolution_depth(ksp, basename):
    """
    Test depth alg on double res images
    """
    print "Double res and gaussian filter"
    ksplarge = np.zeros(np.array(ksp.shape) * 2, dtype=np.complex64)
    szmin = np.array(ksp.shape) / 2 - 1
    szmax = np.array(ksp.shape) + szmin
    ksplarge[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]] = ksp
    image_filtered = fftshift(ifftn(ifftshift(ksplarge)))
    test_depth_algorithm(image_filtered, basename)


def double_resolution(ksp, basename):
    """double_resolution creates double resolution image from k-space data
    based on super-resolution methods for multiple averages, this just expands the
    kspace data and reconstructs image.  Equivalent to interpolation in image space.

    ifft process works for isotropic images only. use ReadFID.simpleifft for other images
    """
    print "Double res and " + basename + " filter"
    # two 32-bit float
    ksplarge = np.zeros(np.array(ksp.shape) * 2, dtype=np.complex64)
    szmin = np.array(ksp.shape[:3]) / 2 - 1
    szmax = np.array(ksp.shape[:3]) + szmin
    ksplarge[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]] = ksp[:3]
    print "Double resolution k-space created. Starting reconstruction ...(may take some time)"
    image_filtered = fftshift(clifftn(ifftshift(ksplarge)))
    print 'Double res recon ' + basename
    del ksplarge

    print "Saving Double res image: " + basename + " filtered"
    save_nifti(np.abs(image_filtered), basename + '-super')


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


    # initialize the CL context
    #ctx, queue = clinit()
    # create a pure read-only OpenCL buffer
    #clbuf = cl.Buffer(ctx,
    #                cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
    #                hostbuf=ksp)

    
        
    # print "Saving Gaussian image"
    # save_nifti(normalise(np.abs(image_filtered)), 'gauss_fourierimage')

    print "Computing Gaussian filtered2 image from Original image"
    kspgauss = kspacegaussian_filter2(ksp, 0.707)
    image_filtered = simpleifft(procpar, dims, hdr, kspgauss, args)
    # print "Saving Gaussian image"
    save_nifti(normalise(np.abs(image_filtered)), 'gauss_kspimage')

    print "Computing Complex Gaussian filtered image from Original image"
    kspcgauss = kspacecplxgaussian_filter(ksp, 0.707)
    image_filtered = simpleifft(procpar, dims, hdr, kspcgauss, args)
    # print "Saving Gaussian image"
    save_nifti(normalise(np.abs(image_filtered)), 'gauss_kspimage2')

    print "Computing Gaussian 1.0 filtered2 image from Original image"
    kspgauss1 = kspacegaussian_filter2(ksp, 1)
    image_filtered = simpleifft(procpar, dims, hdr, kspgauss1, args)
    # print "Saving Gaussian image"
    save_nifti(normalise(np.abs(image_filtered)), 'gauss_kspimage1')

    print "Computing Gaussian 2.0 filtered2 image from Original image"
    kspgauss2 = kspacegaussian_filter2(ksp, 2)
    image_filtered = simpleifft(procpar, dims, hdr, kspgauss2, args)
    # print "Saving Gaussian image"
    save_nifti(normalise(np.abs(image_filtered)), 'gauss_kspimage2')

    print "Computing Gaussian sub-band1.5 image from Original image"
    Fsubband = fouriergausssubband15(ksp.shape, 0.707)
    image_filtered = simpleifft(procpar, dims, hdr, (ksp * Fsubband), args)
    # print "Saving Gaussian image"
    save_nifti(normalise(np.abs(image_filtered)), 'gauss_subband')

    # inhomogeneous correction
    image_corr = inhomogeneouscorrection(ksp, ksp.shape, 3.0 / 60.0)
    # print "Saving Correction image"
    save_nifti(np.abs(image_filtered / image_corr), 'image_inhCorr3')

    # print "Computing Laplacian enhanced image"
    laplacian = simpleifft(
        procpar, dims, hdr, (kspgauss * fourierlaplace(ksp.shape)), args)
    alpha = ndimage.mean(np.abs(image_filtered)) / \
        ndimage.mean(np.abs(laplacian))
    kspgauss = kspacegaussian_filter2(ksp, 1.707)
    image_filtered = simpleifft(procpar, dims, hdr, kspgauss, args)
    image_filtered = (np.abs(image_filtered))
    image_filtered = normalise(image_filtered)
    image_lfiltered = image_filtered - 0.5 * alpha * laplacian
    print '''Saving enhanced image g(x, y, z) = f(x, y, z) -
    Laplacian[f(x, y, z)]'''
    save_nifti(np.abs(image_lfiltered), 'laplacian_enhanced')

    # print "Computing Laplace of Gaussian smoothed image"
    Flaplace = fourierlaplace(ksp.shape)
    # Flaplace = (Flaplace - ndimage.minimum(Flaplace)) /
    #  (ndimage.maximum(Flaplace)
    #  - ndimage.minimum(Flaplace))
    Fsmooth = fouriergauss(
        ksp.shape, (4.0 * np.sqrt(2.0 * np.log(2.0))))
    # (Fgauss/ndimage.maximum(Fgauss))
    Fgauss = fouriergauss(ksp.shape, 0.707)
    laplacian = simpleifft(
        procpar, dims, hdr, (kspgauss * Flaplace *
                             (Fsmooth / ndimage.maximum(Fsmooth))), args)
    laplacian = normalise(laplacian)
    print "Saving Smoothed Gauss Laplacian"
    save_nifti(np.abs(laplacian), 'kspLog_smoothed')

    # # del image_filtered, image_lfiltered

    # #print "Computing Gaussian Laplace image from Smoothed image"
    ksplog = kspacelaplacegaussian_filter(ksp, 0.9)
    image_Log = simpleifft(procpar, dims, hdr, (ksplog), args)
    image_Log = (np.abs(image_Log))
    image_Log = normalise(image_Log)
    save_nifti(np.abs(image_Log), 'kspLog_image')
    # out_img = ndimage.fourier.fourier_gaussian(image_filtered, 2)
    # save_nifti(np.abs(out_img), 'kspLog_smooth')

    print "Computing Epanechnikov filtered image from Original image"
    kspepan = kspaceepanechnikov_filter(ksp, np.sqrt(7.0 / 2.0))
    image_filtered = simpleifft(procpar, dims, hdr, kspepan, args)
    # print "Saving Gaussian image"
    save_nifti(normalise(np.abs(image_filtered)), 'epan_kspimage')


#    test_double_resolution(ksp, 'Raw')
#    test_double_resolution(kspgauss, 'Gauss')
#    test_double_resolution(kspepan, 'Epan')
#    test_double_resolution(ksp*Fsubband, 'GaussSub')
#    test_double_resolution(ksplog, 'LoG')

#    test_depth_algorithm(simpleifft(kspgauss),
#    'gauss_kspdepth')
#    test_double_resolution_depth(kspgauss, 'gauss_depth_large')



    tic()
    imggauss = kspacegaussian_filter_pyfftCL(ksp,np.ones(3))
    print 'PyFFT +OpenCL Gaussian filter:'
    toc()


    
    tic()
    imggauss = kspacegaussian_filter_CL(ksp,np.ones(3))
    print 'Reikna OpenCL Gaussian+recon+ numpy fftshift: first run'
    toc()
    tic()
    imggauss2 = kspacegaussian_filter_CL2(ksp,np.ones(3))
    print 'Reikna OpenCL Gaussian+recon+ Reikna FFTShift: first run'
    toc()


    
#    tic()
#    imggauss = kspacegaussian_filter_reiknaCL(ksp,np.ones(3))
#    print 'Reikna FFT +OpenCL Gaussian filter:'
#    toc()

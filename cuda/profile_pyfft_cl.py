#!/usr/bin/env python
"""
   Profiling OpenCL FFT methods for fourier-domain filtering
   of 3D k-space data


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

import pstats, cProfile


import os
import argparse
import ReadProcpar as Procpar
import ProcparToDicomMap
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

args = parser.parse_args(['--inputdir','../ExampleAgilentData/kidney512iso_01.fid/'])

procpar, procpartext = Procpar.ReadProcpar(
    os.path.join(args.inputdir, 'procpar'))
ds, MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
files = os.listdir(args.inputdir)
fidfiles = [f for f in files if f.endswith('fid')]
print "Number of FID files ", len(fidfiles)
print "Reading FID"
filename = fidfiles[len(fidfiles) - 1]
args.verbose=1
pp, hdr, dims, data_real, data_imag = readfid(args.inputdir,
                                              procpar, args)

print "Echoes: ", hdr['nEchoes'], " Channels: ", hdr['nChannels']
affine = np.eye(4)
import kspace_filter as KSP

image, ksp = recon(pp, dims, hdr,
                   data_real,
                   data_imag, args)
#del data_real, data_imag
print "Shift kspace centre to max point"
ksp = KSP.kspaceshift(ksp)
# (uu,vv,ww) = KSP.fouriercoords(ksp.shape)

def tic():
    #Homemade version of matlab tic and toc functions
    #https://stackoverflow.com/questions/5849800/tic-toc-functions-analog-in-python
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print "Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds."
        return str(time.time() - startTime_for_tictoc)
    else:
        print "Toc: start time not set"
        return ""

import numpy as np

import numpy as np

from pyfft.cl import Plan
import pyopencl as cl
import pyopencl.array as cl_array

sigma= np.ones(3)

def kspacegaussian_filter_CL(ksp,sigma,ctx):
    sz=ksp.shape
    dtype = np.complex64
    ftype = np.float32
    #api = cluda.ocl_api()
    #ctx = cl.create_some_context(0)
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
    #ctx = cl.create_some_context(interactive=False)
    #queue = cl.CommandQueue(ctx)
    w = h = k = 512
    plan = Plan((w,h,k), normalize=True, queue=queue)
    data2_dev = cl_array.to_device(queue, ksp_out)
    plan.execute(data2_dev.data, inverse=True)
    result = data2_dev.get()
    result = np.fft.fftshift(result)
    queue.finish()
    return result  #,ksp_out

plats = cl.get_platforms()
ctx = cl.Context(properties=[(cl.context_properties.PLATFORM, plats[0])])

tic()
imggauss = kspacegaussian_filter_CL(ksp,np.ones(3),ctx)
print 'PyFFT +OpenCL Gaussian filter:'
toc()
tic()
print 'Complex K-space filter + Numpy IFFT'
kspgauss2 = KSP.kspacegaussian_filter2(ksp, 1)
image_filtered = simpleifft(procpar, dims, hdr, kspgauss2, args)
toc()

##  PYFFT



tic()
#ctx = cl.create_some_context(interactive=False)
#queue = cl.CommandQueue(ctx)
w = h = k = 512
plan = Plan((w,h,k), normalize=True, queue=queue)
gpu_data = cl_array.to_device(queue, ksp)
plan.execute(gpu_data.data, inverse=True) 
result = gpu_data.get()
toc()
result = np.fft.fftshift(result)
print "PyFFT OpenCL IFFT time and first three results:"
print "%s sec, %s" % (toc(), str(np.abs(result[:3,0,0])))

tic()
reference = np.fft.fftshift(np.fft.ifftn(ksp))
print "Numpy IFFTN time and first three results:"
print "%s sec, %s" % (toc(), str(np.abs(reference[:3,0,0])))


print "Calulating L1 norm "
print np.linalg.norm(result - reference) / np.linalg.norm(reference)


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# imgplot = plt.imshow(np.abs(result[:,:, 250]), aspect='auto');plt.savefig('epank.jpg')
# clear the plot


f, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(2,3, sharex='col', sharey='row')
ax1.imshow(np.abs(ksp[:,:, 250]), aspect='auto')
ax1.set_title('PyFFT + PyOpenCL Complex K-space filter')
ax2.imshow(np.log10(np.abs(result[:,:, 250]/512**3)), aspect='auto')
ax3.imshow(np.log10(np.abs(result[:,:, 250]/512**3)-np.abs(reference[:,:, 250])), aspect='auto')
ax4.imshow(np.log10(np.abs(reference[:,:, 250])), aspect='auto')
ax5.imshow(np.squeeze(np.abs(imggauss[:,:,250])), aspect='auto')
ax6.imshow(np.squeeze(np.abs(image_filtered[:,:,250])), aspect='auto')

plt.savefig('results_pyfftCL.jpg')

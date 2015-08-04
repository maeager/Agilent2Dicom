
import pstats
import cProfile


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

args = parser.parse_args(
    ['--inputdir', '../ExampleAgilentData/kidney512iso_01.fid/'])

procpar, procpartext = Procpar.ReadProcpar(
    os.path.join(args.inputdir, 'procpar'))
ds, MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
files = os.listdir(args.inputdir)
fidfiles = [f for f in files if f.endswith('fid')]
print "Number of FID files ", len(fidfiles)
print "Reading FID"
filename = fidfiles[len(fidfiles) - 1]
args.verbose = 1
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
# Fgauss = KSP.fouriergauss(ksp.shape, np.array((1,1,1)))


def tic():
    # Homemade version of matlab tic and toc functions
    # https://stackoverflow.com/questions/5849800/tic-toc-functions-analog-in-python
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
# REINKA

from reikna import cluda
from reikna.cluda import functions, dtypes
from reikna.cluda import dtypes, any_api
from reikna.fft import FFT
from reikna.fft import FFTShift
from reikna.core import Annotation, Type, Transformation, Parameter
from pyopencl.tools import clear_first_arg_caches
sigma = np.ones(3)


def kspacegaussian_filter_CL(ksp, sigma):
    sz = ksp.shape
    dtype = np.complex64
    ftype = np.float32
    #api = cluda.ocl_api()
    api = any_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(ksp)
    ifft = FFT(data_dev)
    FACTOR = 1.0
    program = thr.compile("""
KERNEL void gauss_kernel(
    GLOBAL_MEM ${ctype} *dest,
    GLOBAL_MEM ${ctype} *src)
{
  const ${ultype} x = (${ultype})get_global_id(0);
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
""" % (sz[0], sz[1], sz[2], sigma[0], sigma[1], sigma[2], FACTOR),
        render_kwds=dict(ctype=dtypes.ctype(dtype),
                         ftype=dtypes.ctype(ftype),
                         ultype=dtypes.ctype(np.uint64),
                         exp=functions.exp(ftype)), fast_math=True)
    gauss_kernel = program.gauss_kernel
    #data_dev = thr.empty_like(ksp_dev)
    gauss_kernel(data_dev, data_dev, global_size=sz[0] * sz[1] * sz[2])

    thr.synchronize()
    ##
    #api = any_api()
    #thr = api.Thread.create()
    #data_dev = thr.to_device(ksp_out)
    ifft = FFT(data_dev)
    cifft = ifft.compile(thr)
    cifft(data_dev, data_dev, inverse=0)
    result = np.fft.fftshift(data_dev.get() / sz[0] * sz[1] * sz[2])
    result = result[::-1, ::-1, ::-1]
    result = np.roll(np.roll(np.roll(result, 1, axis=2), 1, axis=1), 1, axis=0)
    return result  # ,ksp_out


def kspacegaussian_filter_CL2(ksp, sigma):
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
    sz = ksp.shape
    dtype = np.complex64
    ftype = np.float32
    ultype = np.uint64
    #api = cluda.ocl_api()
    api = any_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(ksp)
    ifft = FFT(data_dev)
    FACTOR = 1.0
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
""" % (sz[0], sz[1], sz[2], sigma[0], sigma[1], sigma[2], FACTOR),
        render_kwds=dict(ctype=dtypes.ctype(dtype),
                         ftype=dtypes.ctype(ftype),
                         exp=functions.exp(ftype)), fast_math=True)
    gauss_kernel = program.gauss_kernel
    #data_dev = thr.empty_like(ksp_dev)
    gauss_kernel(data_dev, data_dev, global_size=sz[0] * sz[1] * sz[2])

    thr.synchronize()
    # Recon
    #data_dev = thr.to_device(ksp)
    ifftobj = FFT(data_dev)
    cifft = ifftobj.compile(thr)
    fftshiftobj = FFTShift(data_dev)
    cfftshift = fftshiftobj.compile(thr)
    cifft(data_dev, data_dev, inverse=0)
    thr.synchronize()
    cfftshift(data_dev, data_dev)
    thr.synchronize()
    result2 = data_dev.get() / np.prod(np.array(ksp.shape))
    result2 = result2[::-1, ::-1, ::-1]
    thr.release()
    return result2


tic()
imggauss = kspacegaussian_filter_CL(ksp, np.ones(3))
print 'Reikna OpenCL Gaussian+recon+ numpy fftshift: first run'
toc()
tic()
imggauss2 = kspacegaussian_filter_CL2(ksp, np.ones(3))
print 'Reikna OpenCL Gaussian+recon+ Reikna FFTShift: first run'
toc()


tic()
kspgauss2 = KSP.kspacegaussian_filter2(ksp, 1)
image_filtered = simpleifft(procpar, dims, hdr, kspgauss2, args)
toc()


from reikna.cluda import dtypes, any_api
from reikna.fft import FFT
from reikna.core import Annotation, Type, Transformation, Parameter
# create two timers so we can speed-test each approach

api = any_api()
thr = api.Thread.create()
N = 512

tic()
data_dev = thr.to_device(ksp)
ifft = FFT(data_dev)
cifft = ifft.compile(thr)
cifft(data_dev, data_dev, inverse=0)
thr.synchronize()
toc()
result = np.fft.fftshift(data_dev.get() / N**3)
result = result[::-1, ::-1, ::-1]
result = np.roll(np.roll(np.roll(result, 1, axis=2), 1, axis=1), 1, axis=0)
print "Reikna IFFT time and first three results:"
print "%s sec, %s" % (toc(), str(np.abs(result[:3, 0, 0])))
thr.release()
del ifft, cifft, data_dev, thr

thr = api.Thread.create()
tic()
data_dev = thr.to_device(ksp)
ifft = FFT(data_dev)
cifft = ifft.compile(thr)
fftshiftobj = FFTShift(data_dev)
cfftshift = fftshiftobj.compile(thr)
cifft(data_dev, data_dev, inverse=0)
thr.synchronize()
toc()
cfftshift(data_dev, data_dev)
thr.synchronize()
result2 = data_dev.get() / N**3
result2 = result2[::-1, ::-1, ::-1]
#result = np.roll(np.roll(np.roll(result,1,axis=2),1,axis=1),1,axis=0)
print "Reikna IFFT time and first three results:"
print "%s sec, %s" % (toc(), str(np.abs(result2[:3, 0, 0])))
thr.release()
del ifft, cifft, data_dev, fftshiftobj, cfftshift

del thr, api

tic()
reference = np.fft.fftshift(np.fft.ifftn(ksp))
print "Numpy IFFTN time and first three results:"
print "%s sec, %s" % (toc(), str(np.abs(reference[:3, 0, 0])))

# print np.linalg.norm(imggauss-image_filtered) /
# np.linalg.norm(image_filtered)
print np.linalg.norm((np.abs(imggauss2)) - np.abs(image_filtered)) / np.linalg.norm(np.abs(image_filtered))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# imgplot = plt.imshow(np.abs(result[:,:, 250]), aspect='auto');plt.savefig('epank.jpg')
# clear the plot


f, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(
    2, 3, sharex='col', sharey='row')
ax1.imshow(np.abs(ksp[:, :, 250]), aspect='auto')
ax1.set_title('Sharing x per column, y per row')
ax2.imshow((np.abs(result[:, :, 250] / 512**3)), aspect='auto')
ax3.imshow(np.log10(np.abs(imggauss2[
           :, :, 250])) - np.log10(np.abs(image_filtered[:, :, 250])), aspect='auto')
ax4.imshow((np.abs(reference[:, :, 250])), aspect='auto')
ax5.imshow(np.squeeze(np.abs(imggauss2[:, :, 250])), aspect='auto')
ax6.imshow(np.squeeze(np.abs(image_filtered[:, :, 250])), aspect='auto')
plt.savefig('results_reiknaCL.jpg')


# Epanechnikov

from cplxfilter import epanechnikov_kernel


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

    def is_odd(num):
        return num & 0x1

    from cplxfilter import epanechnikov_kernel
    if not hasattr(sigma, "__len__"):
        Kepa = epanechnikov_kernel(
            (np.ceil(sigma) + 1, np.ceil(sigma) + 1, np.ceil(sigma) + 1),
            sigma)
    else:
        print (np.ceil(sigma[0]) + 2,
               np.ceil(sigma[1]) + 2, np.ceil(sigma[2]) + 2)
        print sigma
        fsiz = (np.ceil(sigma) + 2).astype(int)
        for i in xrange(0, fsiz.size):
            if is_odd(fsiz[i]):
                fsiz[i] += 1

        Kepa = epanechnikov_kernel((np.ceil(sigma[0]) + 2, np.ceil(sigma[1]) + 2,
                                    np.ceil(sigma[2]) + 2), sigma)

    Kfilter = np.zeros(np.array(siz), dtype=np.complex64)
    szmin = np.floor(
        np.array(siz) / 2.0 - np.floor(np.array(Kepa.shape) / 2.0) - 1)
    szmax = np.floor(szmin + np.array(Kepa.shape))
    print "Epa filter size ", siz, " image filter ", Kepa.shape, " szmin ", szmin, " szmax ", szmax
    Kfilter[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]] = Kepa
    Kfilter[szmin[0]:szmax[0], szmin[1]:szmax[
        1], szmin[2]:szmax[2]].imag = Kepa

    # return np.abs(fftshift(clfftn(Kfilter)))

    api = any_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(Kfilter)
    fft = FFT(data_dev)
    cfft = fft.compile(thr)
    fftshift = FFTShift(data_dev)
    cfftshift = fftshift.compile(thr)
    cfft(data_dev, data_dev)
    thr.synchronize()
    cfftshift(data_dev, data_dev)
    thr.synchronize()

    result2 = data_dev.get()  # / np.prod(np.array(ksp.shape))
    #result2 = result2[::-1,::-1,::-1]
    thr.release()
    result = np.zeros(np.array(siz), dtype=np.complex64)
    result.real = np.abs(result2) / np.sqrt(2)
    result.imag = np.abs(result2) / np.sqrt(2)
    return result


sigma = np.ones(3) * np.sqrt(7)
# Fepanechnikov= fourierepanechnikov(ksp.shape,sigma)

from cplxfilter import epanechnikov_kernel

# Wolfram alpha:
#     Abs[FourierTransform[(1+i)*UnitBox[x/2]*(1-x^2)*0.75]]
#    => 0.423142 abs((4 sin(omega)-4 omega cos(omega))/omega^3)


def kspaceepanechnikov_filter(ksp, sigma):
    """ Kspace gaussian filter and recon using GPU OpenCL

    1. GPU intialisation
    2. push KSP complex matrix to GPU
    3. declare FFT program
    4. declare Complex Epan GPU filter program
    5. Execute Epan GPU program
    6. GPU sync
    7. Execute FFT Recon
    8. Execute FFTshift
    9. Retrieve reconstruced complex image from GPU
    10. Reorganise image to standard (mimic numpy format)

    """
    sz = ksp.shape
    dtype = np.complex64
    ftype = np.float32
    ultype = np.uint64
    #api = cluda.ocl_api()
    api = any_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(ksp)
    ifft = FFT(data_dev)
    FACTOR = 1.0
    program = thr.compile("""
KERNEL void epan_kernel(
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
  ${ftype} omega = (i*sigma[0]+j*sigma[1]+k*sigma[2]);
  ${ftype} omega3 = ((i*sigma[0])*(i*sigma[0])*(i*sigma[0])+(j*sigma[1])*(j*sigma[1])*(j*sigma[1])+(k*sigma[2])*(k*sigma[2])*(k*sigma[2]));        
  ${ftype} weight = 0.423142 * fabs((4 * sin(omega) - 4 * omega * cos(omega)) / omega3);
  dest[idx].x = src[idx].x * weight; 
  dest[idx].y = src[idx].y * weight; 
  
}
""" % (sz[0], sz[1], sz[2], sigma[0], sigma[1], sigma[2], FACTOR),
        render_kwds=dict(ctype=dtypes.ctype(dtype),
                         ftype=dtypes.ctype(ftype)), fast_math=True)
    epan_kernel = program.epan_kernel
    #data_dev = thr.empty_like(ksp_dev)
    epan_kernel(data_dev, data_dev, global_size=sz[0] * sz[1] * sz[2])
    return data_dev()
    # thr.synchronize()
    # ## Recon
    # #data_dev = thr.to_device(ksp)
    # ifftobj = FFT(data_dev)
    # cifft = ifftobj.compile(thr)
    # fftshiftobj = FFTShift(data_dev)
    # cfftshift = fftshiftobj.compile(thr)
    # cifft(data_dev, data_dev,inverse=0)
    # thr.synchronize()
    # cfftshift(data_dev,data_dev)
    # thr.synchronize()
    # result2 = data_dev.get() / np.prod(np.array(ksp.shape))
    # result2 = result2[::-1,::-1,::-1]
    # thr.release()
    # return result2


def kspaceepanechnikov_filter_CL2(ksp, sigma):
    sz = ksp.shape
    dtype = np.complex64
    ftype = np.float32
    clear_first_arg_caches()
    fsiz = (5, 5, 5)
    print (np.ceil(sigma[0]) + 2,
           np.ceil(sigma[1]) + 2, np.ceil(sigma[2]) + 2)
    print sigma
    fsiz = (np.ceil(sigma) + 2).astype(int)
    for i in xrange(0, fsiz.size):
        if not fsiz[i] & 0x1:
            fsiz[i] += 1
    # Create image-domain Epanechikov kernel
    Kepa = epanechnikov_kernel(fsiz, sigma)
    # Place kernel at centre of ksp-sized matrix
    Kfilter = np.zeros(np.array(sz), dtype=np.complex64)
    szmin = np.floor(
        np.array(sz) / 2.0 - np.floor(np.array(Kepa.shape) / 2.0) - 1)
    szmax = np.floor(szmin + np.array(Kepa.shape))
    print "Epa filter size ", sz, " image filter ", Kepa.shape, " szmin ", szmin, " szmax ", szmax
    Kfilter[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]] = Kepa
    Kfilter[szmin[0]:szmax[0], szmin[1]:szmax[
        1], szmin[2]:szmax[2]].imag = Kepa
    # Create fourier-domain Epanechnikov filter
    api = any_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(Kfilter)
    rfft = FFT(data_dev)
    crfft = rfft.compile(thr)
    fftshift = FFTShift(data_dev)
    cfftshift = fftshift.compile(thr)
    crfft(data_dev, data_dev)
    thr.synchronize()
    cfftshift(data_dev, data_dev)
    Fepanechnikov = np.abs(data_dev.get())  # / np.prod(np.array(ksp.shape))
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
}""", render_kwds=dict(ctype=dtypes.ctype(dtype), ftype=dtypes.ctype(ftype)))

    data_dev = thr.to_device(ksp)
    filter_dev = thr.to_device(Fepanechnikov)
    multiply_them = program.multiply_them
    multiply_them(data_dev, data_dev, filter_dev, global_size=512 * 512 * 512)
    thr.synchronize()
    del filter_dev, program
    #api = cluda.ocl_api()
    #api = any_api()
    #thr = api.Thread.create()
    # Filter
    # data_dev = thr.to_device(ksp)
    # ifft = FFT(data_dev)
    FACTOR = 1.0

    # Recon
    # thr.synchronize()
    #data_dev = thr.to_device(ksp)
    ifft = FFT(data_dev)
    cifft = ifft.compile(thr)
    fftshiftobj = FFTShift(data_dev)
    cfftshift = fftshiftobj.compile(thr)
    cifft(data_dev, data_dev, inverse=0)
    thr.synchronize()
    cfftshift(data_dev, data_dev)
    thr.synchronize()
    result2 = data_dev.get() / np.prod(np.array(ksp.shape))
    result2 = result2[::-1, ::-1, ::-1]
    thr.release()
    return result2

sigma = np.ones(3) * np.sqrt(7)
tic()
kimage_epan = kspaceepanechnikov_filter_CL2(ksp, sigma)
print 'GPU Epanechnikov filter'
toc()
tic()
kspepan2 = KSP.kspacegaussian_filter2(ksp, 1)
image_filtered = simpleifft(procpar, dims, hdr, kspepan2, args)
print 'CPU Epanechnikov filter'
toc()


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# imgplot = plt.imshow(np.abs(result[:,:, 250]), aspect='auto');plt.savefig('epank.jpg')
# clear the plot


f, ((ax1, ax3, ax5), (ax2, ax4, ax6)) = plt.subplots(
    2, 3, sharex='col', sharey='row')
ax1.imshow(np.abs(ksp[:, :, 250]), aspect='auto')
ax1.set_title('KSP data')
ax3.set_title('GPU recon')
ax3.imshow(np.squeeze(np.abs(kimage_epan[250, :, :])), aspect='auto')
ax2.set_title('GPU - CPU ')
ax2.imshow((np.abs(kimage_epan[:, :, 250])) -
           (np.abs(image_filtered[:, :, 250])), aspect='auto')
ax5.set_title('CPU recon')
ax5.imshow(np.squeeze(np.abs(image_filtered[250, :, :])), aspect='auto')

ax4.imshow(np.squeeze(np.abs(kimage_epan[:, :, 250])), aspect='auto')
ax6.imshow(np.squeeze(np.abs(image_filtered[:, :, 250])), aspect='auto')
plt.savefig('results_reiknaCL_epan.jpg')

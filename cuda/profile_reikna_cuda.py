
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
## REINKA
import pycuda.driver as drv
import pycuda.tools
import pycuda.autoinit
from reikna.cluda import dtypes, cuda_api, functions
from reikna.fft import FFT
from reikna.fft import FFTShift
from reikna.core import Annotation, Type, Transformation, Parameter

#from pyopencl.tools import clear_first_arg_caches 
sigma = np.ones(3)

def kspacegaussian_filter_CL(ksp,sigma):
    sz=ksp.shape
    dtype = np.complex64
    ftype = np.float32
    #api = cluda.ocl_api()
    api = cuda_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(ksp)
    ifft = FFT(data_dev)
    FACTOR=1.0
    program = thr.compile("""
KERNEL void gauss_kernel(
    GLOBAL_MEM ${ctype} *dest,
    GLOBAL_MEM ${ctype} *src)
{
  const ulong x = get_global_id(0);const ulong y = get_global_id(1);const ulong z = get_global_id(2);
  const SIZE_T dim1= %d;
  const SIZE_T dim2= %d;
  const SIZE_T dim3= %d;                    
  ${ftype} sigma[3];
  sigma[0]=%f;sigma[1]=%f;sigma[2]=%f;
  ${ftype} factor = %f;            
  const double TWOPISQ = 19.739208802178716; //6.283185307179586;  //2*3.141592;
  const ${ftype} SQRT2PI = 2.5066282746;
  const double CUBEDSQRT2PI = 15.749609945722419;
  const ulong idx = dim1*dim2*z + dim1*y + x;
  ${ftype} i = (${ftype})(x); // )((x / dim3) / dim2);
      i = (i - (${ftype})floor((${ftype})(dim1)/2.0))/(${ftype})(dim1);
  ${ftype} j = (${ftype})(y); //(x / dim3);
      if((SIZE_T)j > dim2) {j=(${ftype})fmod(j, (${ftype})dim2);};
      j = (j - (${ftype})floor((${ftype})(dim2)/2.0f))/(${ftype})(dim2);
  //Account for large global index (stored as ulong) before performing modulus
  //double pre_k=fmod((double)(x) , (double) dim3);
  ${ftype} k = (${ftype}) (z); //pre_k;
      k = (k - (${ftype})floor((${ftype})(dim3)/2.0f))/(${ftype})(dim3);

  ${ftype} weight = exp(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  //${ftype} weight = expm1(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]))+1;
  //${ftype} weight= ${exp}(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  dest[idx].x = src[idx].x * weight;
  dest[idx].y = src[idx].y * weight; 
  
}
""" % (sz[0],sz[1],sz[2],sigma[0],sigma[1],sigma[2],FACTOR),
                          render_kwds=dict(ctype=dtypes.ctype(dtype),
                                           ftype=dtypes.ctype(ftype),
                                           exp=functions.exp(ftype)),fast_math=True)
    gauss_kernel = program.gauss_kernel
    #data_dev = thr.empty_like(ksp_dev)
    gauss_kernel(data_dev, data_dev, global_size=(sz[0],sz[1],sz[2]))
    # ksp_out = data_dev.get()
    thr.synchronize()
    ##
    #api = cuda_api()
    #thr = api.Thread.create()
    #data_dev = thr.to_device(ksp_out)
    ifft = FFT(data_dev)
    cifft = ifft.compile(thr)
    cifft(data_dev, data_dev,inverse=0)
    result = np.fft.fftshift(data_dev.get() / sz[0]*sz[1]*sz[2])
    result = result[::-1,::-1,::-1]
    result = np.roll(np.roll(np.roll(result,1,axis=2),1,axis=1),1,axis=0)
    return result  #,ksp_out

def kspacegaussian_filter_CL2(ksp,sigma):
    sz=ksp.shape
    dtype = np.complex64
    ftype = np.float32
    #api = cluda.ocl_api()
    api = cuda_api()
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
                                           ftype=dtypes.ctype(ftype),
                                           exp=functions.exp(ftype)),fast_math=True)
    gauss_kernel = program.gauss_kernel
    #data_dev = thr.empty_like(ksp_dev)
    gauss_kernel(data_dev, data_dev, global_size=(sz[0],sz[1],sz[2]))
    
    thr.synchronize()
    #data_dev = thr.to_device(ksp)
    ifft = FFT(data_dev)
    cifft = ifft.compile(thr)
    fftshift = FFTShift(data_dev)
    cfftshift = fftshift.compile(thr)
    cifft(data_dev, data_dev,inverse=0)
    thr.synchronize()
    cfftshift(data_dev,data_dev)
    thr.synchronize()
    result2 = data_dev.get() / np.prod(np.array(ksp.shape))
    result2 = result2[::-1,::-1,::-1]
    thr.release()
    return result2

pycuda.tools.clear_context_caches()
    
tic()
imggauss = kspacegaussian_filter_CL(ksp,np.ones(3))
print 'Reikna Cuda Gaussian+recon+ numpy fftshift: first run'
toc()
pycuda.tools.clear_context_caches()
tic()
imggauss2 = kspacegaussian_filter_CL2(ksp,np.ones(3))
print 'Reikna Cuda Gaussian+recon+ Reikna FFTShift: first run'
toc()
pycuda.tools.clear_context_caches()

tic()
kspgauss2 = KSP.kspacegaussian_filter2(ksp, 1)
image_filtered = simpleifft(procpar, dims, hdr, kspgauss2, args)
toc()




# create two timers so we can speed-test each approach
#start = drv.Event()
#end = drv.Event()
api = cuda_api()
thr = api.Thread.create()
N=512

tic()
data_dev = thr.to_device(ksp)
ifft = FFT(data_dev)
cifft = ifft.compile(thr)
thr.synchronize()
cifft(data_dev, data_dev,inverse=0)
result = np.fft.fftshift(data_dev.get() / N**3)
result = result[::-1,::-1,::-1]
result = np.roll(np.roll(np.roll(result,1,axis=2),1,axis=1),1,axis=0)

print "Reikna IFFT time and first three results:"
print "%s sec, %s" % (toc(), str(np.abs(result[:3,0,0])))

tic()
reference = np.fft.fftshift(np.fft.ifftn(ksp))
print "Numpy IFFTN time and first three results:"
print "%s sec, %s" % (toc(), str(np.abs(reference[:3,0,0])))

print np.linalg.norm(result - reference) / np.linalg.norm(reference)
print np.linalg.norm((np.abs(imggauss2))-np.abs(image_filtered)) / np.linalg.norm(np.abs(image_filtered))


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# imgplot = plt.imshow(np.abs(result[:,:, 250]), aspect='auto');plt.savefig('epank.jpg')
# clear the plot


f, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(2,3, sharex='col', sharey='row')
ax1.imshow(np.abs(ksp[:,:, 250]), aspect='auto')
ax1.set_title('Sharing x per column, y per row')
ax2.imshow(np.log10(np.abs(result[:,:, 250]/512**3)), aspect='auto')
ax3.imshow(np.log10(np.abs(result[:,:, 250]/512**3))-np.log10(np.abs(reference[:,:, 250])), aspect='auto')
ax4.imshow(np.log10(np.abs(reference[:,:, 250])), aspect='auto')
ax5.imshow(np.squeeze(np.abs(result[250,:,:])), aspect='auto')
ax6.imshow(np.squeeze(np.abs(reference[250,:,:])), aspect='auto')
plt.savefig('results_reiknaCUDA.jpg')


f, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(2,3, sharex='col', sharey='row')
ax1.set_title('K-space')
ax1.imshow(np.abs(ksp[:,:, 250]), aspect='auto')
ax2.set_title('Reikna CUDA Complex K-space filter\n GPU Recon')
ax2.imshow((np.abs(result[:,:, 250]/512**3)), aspect='auto')
ax3.set_title('Diff (log10)')
ax3.imshow(np.log10(np.abs(imggauss[:,:, 250])-np.abs(image_filtered[:,:, 250])), aspect='auto')
ax4.set_title('Numpy Recon')
ax4.imshow((np.abs(reference[:,:, 250])), aspect='auto')
ax5.set_title('GPU Filtered Recon')
ax5.imshow(np.squeeze(np.abs(imggauss[:,:,250])), aspect='auto')
ax6.set_title('Standard Filter')
ax6.imshow(np.squeeze(np.abs(image_filtered[:,:,250])), aspect='auto')

plt.savefig('results_reiknaCUDA.jpg')

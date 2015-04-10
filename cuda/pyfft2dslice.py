# ----------
import numpy
from pyfft.cuda import Plan
import pycuda.autoinit
import pycuda.gpuarray as gpuarray
# w,h,k are the array dimensions in a power of 2
# im1, im2 are the input 3d arrays of dtype complex64
w = h = k = 1024
im1 = numpy.random.rand(w,h,k).astype(numpy.complex64)
im2 = numpy.random.rand(w,h,k).astype(numpy.complex64)
%time plan = Plan((w,h,k/4), normalize=True)
# forward transform on device
mid = k/4
im1_ft=im1.copy()
for slab in xrange(0,4):
	%time im1_gpu = gpuarray.to_gpu(im1[:,:,slab*mid:mid*(slab+1)].astype(numpy.complex64));
	%time plan.execute(im1_gpu);
	%time im1_ft_tmp = im1_gpu.get();im1_ft[:,:,slab*mid:mid*(slab+1)]=im1_ft_tmp;del im1_gpu
im2_ft=im2.copy()
for slab in xrange(0,4):
	%time im2_gpu = gpuarray.to_gpu(im2[:,:,slab*mid:mid*(slab+1)].astype(numpy.complex64));
	%time plan.execute(im2_gpu);
	%time im2_ft_tmp = im2_gpu.get();im2_ft[:,:,slab*mid:mid*(slab+1)]=im2_ft_tmp;del im2_gpu

corr_gpu=im1.copy()
# do multiplication on host - can be done on device.
%time conv = im1_ft * im2_ft
#inverse transform on device
for slab in xrange(0,4):
	%time conv_gpu = gpuarray.to_gpu(conv[:,:,slab*mid:mid*(slab+1)].astype(numpy.complex64));
	%time plan.execute(conv_gpu);
	%time corr_ft_tmp = conv_gpu.get();corr_gpu[:,:,slab*mid:mid*(slab+1)]=corr_ft_tmp;del conv_gpu

del conv
# Reference calculation on CPU:
%time im1_ft = numpy.fft.fftn(im1)
%time im2_ft = numpy.fft.fftn(im2)
%time conv = im1_ft * im2_ft
del im1
del im2
del im1_ft
del im2_ft
%time corr_cpu = numpy.fft.ifftn(conv)
print numpy.linalg.norm(corr_cpu - corr_gpu) / numpy.linalg.norm(corr_gpu)
# ----------

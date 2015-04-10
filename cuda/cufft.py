#!/usr/bin/env python

"""
Demonstrates how to use the PyCUDA interface to CUFFT to compute a
batch of 2D FFTs.
"""

import pycuda.autoinit
import pycuda.gpuarray as gpuarray
import numpy as np

import scikits.cuda.fft as cu_fft


print 'Testing fft/ifft..'
N = 512
batch_size = 512

x = np.empty((batch_size, N, N), np.float32)
xf = np.empty((batch_size, N, N), np.complex64)
y = np.empty((batch_size, N, N), np.float32)
tic()
for i in xrange(batch_size):
    x[i, :, :] = np.asarray(np.random.rand(N, N), np.float32)
    xf[i, :, :] = np.fft.fft2(x[i, :, :])
    y[i, :, :] = np.real(np.fft.ifft2(xf[i, :, :]))
toc()
tic()
x_gpu = gpuarray.to_gpu(x)
xf_gpu = gpuarray.empty((batch_size, N, N/2+1), np.complex64)
plan_forward = cu_fft.Plan((N, N), np.float32, np.complex64, batch_size)
cu_fft.fft(x_gpu, xf_gpu, plan_forward)

y_gpu = gpuarray.empty_like(x_gpu)
plan_inverse = cu_fft.Plan((N, N), np.complex64, np.float32, batch_size)
cu_fft.ifft(xf_gpu, y_gpu, plan_inverse, True)
toc()
print 'Success status: ', np.allclose(y, y_gpu.get(), atol=1e-6)

print 'Testing in-place fft..'
tic()
x = np.empty((batch_size, N, N), np.complex64)
x = np.asarray(np.random.rand(batch_size, N, N), np.complex64)
x.imag = np.random.rand(batch_size, N, N)
x_gpu = gpuarray.to_gpu(x)
toc()
tic()
plan = cu_fft.Plan((N, N), np.complex64, np.complex64, batch_size)
cu_fft.fft(x_gpu, x_gpu, plan)
cu_fft.ifft(x_gpu, x_gpu, plan, True)
toc()
print 'Success status: ', np.allclose(x, x_gpu.get(), atol=1e-6)




# Reconstruction from complex kspace
print 'Testing recon ifft..'
tic()
x = np.empty((batch_size, N, N), np.complex64)
x = np.asarray(np.random.rand(batch_size, N, N), np.complex64)
x.imag = np.random.rand(batch_size, N, N)
x_gpu = gpuarray.to_gpu(x)
# xf_gpu = gpuarray.empty((batch_size, N, N), np.complex64)
plan = cu_fft.Plan((batch_size, N, N), np.complex64, np.complex64, 1)
toc()
tic()
cu_fft.ifft(x_gpu, x_gpu, plan, True)
toc()
tic()
y = numpy.fft.ifftn(x)
toc()
print 'Success status: ', np.allclose(y, x_gpu.get(), atol=1e-6)

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
import pycuda.autoinit
import pycuda.gpuarray as gpuarray
import scikits.cuda.fft as cu_fft


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


def gaussian_fourierkernel_elemwise(uu, vv, ww, sigma):
    """
    Create Gaussian Fourier filter kernel
    Element wise cuda implementation

    """
    import pycuda.gpuarray as gpuarray
    import pycuda.driver as cuda
    import pycuda.autoinit

    siz = np.floor(np.array(uu.shape))
    zz = uu.copy()
    u_gpu = gpuarray.to_gpu(uu)
    v_gpu = gpuarray.to_gpu(vv)
    w_gpu = gpuarray.to_gpu(ww)
    from pycuda.elementwise import ElementwiseKernel
    norm_comb = ElementwiseKernel(
        "float s, float pi, float *u, float *v, float *w, float *z",
        "z[i] = exp(-2 * (pi ^ 2) * (u[i] ^ 2 + v[i] ^ 2 + w[i] ^ 2) * (s^2))",
        "normal_combination")

    z_gpu = gpuarray.empty_like(a_gpu)
    norm_comb(sigma, np.pi, u_gpu, v_gpu, w_gpu, z_gpu)

    gfilter = (np.exp(-2 * (np.pi ** 2)) * z_gpu).get()
    return gfilter


def gaussian_fourierkernel_elemwise_flatten(uu, vv, ww, sigma):
    """
    Create Gaussian Fourier filter kernel
    Element wise cuda implementation

    """
    import pycuda.gpuarray as gpuarray
    import pycuda.driver as cuda
    import pycuda.autoinit
    import numpy as np
    siz = np.floor(np.array(uu.shape) / 2 + 1)
    zz = uu[:siz[1], :siz[2], :siz[3]].flatten()
    u_gpu = gpuarray.to_gpu(uu[:siz[1], :siz[2], :siz[3]].flatten())
    v_gpu = gpuarray.to_gpu(vv[:siz[1], :siz[2], :siz[3]].flatten())
    w_gpu = gpuarray.to_gpu(ww[:siz[1], :siz[2], :siz[3]].flatten())
    from pycuda.elementwise import ElementwiseKernel
    norm_comb = ElementwiseKernel(
        "float s, float pi, float *u, float *v, float *w, float *z",
        "z[i] = exp(-2 * (pi ^ 2) * (u[i] ^ 2 + v[i] ^ 2 + w[i] ^ 2) * (s^2))",
        "normal_combination")

    z_gpu = gpuarray.empty_like(a_gpu)
    norm_comb(sigma, np.pi, u_gpu, v_gpu, w_gpu, z_gpu)
    zz = z_gpu.get()
    gfilter = uu.copy()

    gfilter[:siz[1], :siz[2], :siz[3]] = np.reshape(zz, siz)
    gfilter[siz[1]:, :siz[2], :siz[3]] = np.roll(
        np.reshape(zz, siz), siz[1], 1)
    gfilter[:siz[1], siz[2]:, :siz[3]] = np.roll(
        np.reshape(zz, siz), siz[2], 2)
    gfilter[:siz[1], :siz[2], siz[3]:] = np.roll(
        np.reshape(zz, siz), siz[3], 3)
    gfilter[siz[1]:, siz[2]:, :siz[3]] = np.roll(
        np.roll(np.reshape(zz, siz - 1), siz[1], 1), siz[2], 2)
    gfilter[:, :, siz[3]:] = gfilter[:, :, :siz[3]]
    gfilter = gfilter / zz[-1]
    if not hasattr(sigma, "__len__"):  # type(sigma) is float:
        gfilter = np.exp(-2 * (np.pi ** 2) *
                         (uu ** 2 + vv ** 2 + ww ** 2) * (sigma ** 2))

        midpoint = np.ceil(np.array(uu.shape) / 2.0)
        maxval = ndimage.maximum(gfilter[midpoint[0] - 10:midpoint[0] + 10,
                                         midpoint[1] - 10:midpoint[1] + 10,
                                         midpoint[2] - 10:midpoint[2] + 10])
        return gfilter / maxval
    elif len(sigma) == 2:
        gfilter = np.exp(-2 * (np.pi ** 2) * ((sigma[0] ** 2) * uu ** 2 +
                                              (sigma[1] ** 2) * vv ** 2))
        midpoint = np.ceil(np.array(uu.shape) / 2.0)
        maxval = ndimage.maximum(gfilter[midpoint[0] - 10:midpoint[0] + 10,
                                         midpoint[1] - 10:midpoint[1] + 10])
        gfilter = gfilter / maxval
    else:
        gfilter = np.exp(-2 * (np.pi ** 2) * ((sigma[0] ** 2) * uu ** 2 +
                                              (sigma[1] ** 2) * vv ** 2 +
                                              (sigma[2] ** 2) * ww ** 2))
        midpoint = np.ceil(np.array(uu.shape) / 2.0)
        maxval = ndimage.maximum(gfilter[midpoint[0] - 10:midpoint[0] + 10,
                                         midpoint[1] - 10:midpoint[1] + 10,
                                         midpoint[2] - 10:midpoint[2] + 10])
        gfilter = gfilter / maxval
    return gfilter


def gaussian_fourierkernel_old(uu, vv, ww, sigma):
    """
    Create Gaussian Fourier filter kernel
    Relegated: numpy.exp too slow for values close to zero.
    """
    if not hasattr(sigma, "__len__"):  # type(sigma) is float:
        gfilter = np.exp(-2 * (np.pi ** 2) *
                         (uu ** 2 + vv ** 2 + ww ** 2) * (sigma ** 2))

        midpoint = np.ceil(np.array(uu.shape) / 2.0)
        maxval = ndimage.maximum(gfilter[midpoint[0] - 10:midpoint[0] + 10,
                                         midpoint[1] - 10:midpoint[1] + 10,
                                         midpoint[2] - 10:midpoint[2] + 10])
        return gfilter / maxval
    elif len(sigma) == 2:
        gfilter = np.exp(-2 * (np.pi ** 2) * ((sigma[0] ** 2) * uu ** 2 +
                                              (sigma[1] ** 2) * vv ** 2))
        midpoint = np.ceil(np.array(uu.shape) / 2.0)
        maxval = ndimage.maximum(gfilter[midpoint[0] - 10:midpoint[0] + 10,
                                         midpoint[1] - 10:midpoint[1] + 10])
        gfilter = gfilter / maxval
    else:
        gfilter = np.exp(-2 * (np.pi ** 2) * ((sigma[0] ** 2) * uu ** 2 +
                                              (sigma[1] ** 2) * vv ** 2 +
                                              (sigma[2] ** 2) * ww ** 2))
        midpoint = np.ceil(np.array(uu.shape) / 2.0)
        maxval = ndimage.maximum(gfilter[midpoint[0] - 10:midpoint[0] + 10,
                                         midpoint[1] - 10:midpoint[1] + 10,
                                         midpoint[2] - 10:midpoint[2] + 10])
        gfilter = gfilter / maxval
    return gfilter


def gaussian_fourierkernel(uu, vv, ww, sigma):
    """
    Create Gaussian Fourier filter kernel
    """
    if not hasattr(sigma, "__len__"):  # type(sigma) is float:
        gfilter = np.expm1(-2 * (np.pi ** 2) *
                           (uu ** 2 + vv ** 2 + ww ** 2) * (sigma ** 2)) + 1

        midpoint = np.ceil(np.array(uu.shape) / 2.0)
        maxval = ndimage.maximum(gfilter[midpoint[0] - 10:midpoint[0] + 10,
                                         midpoint[1] - 10:midpoint[1] + 10,
                                         midpoint[2] - 10:midpoint[2] + 10])
        return gfilter / maxval
    elif len(sigma) == 2:
        gfilter = np.expm1(-2 * (np.pi ** 2) * ((sigma[0] ** 2) * uu ** 2 +
                                                (sigma[1] ** 2) * vv ** 2)) + 1
        midpoint = np.ceil(np.array(uu.shape) / 2.0)
        maxval = ndimage.maximum(gfilter[midpoint[0] - 10:midpoint[0] + 10,
                                         midpoint[1] - 10:midpoint[1] + 10])
        gfilter = gfilter / maxval
    else:
        gfilter = np.expm1(-2 * (np.pi ** 2) * ((sigma[0] ** 2) * uu ** 2 +
                                                (sigma[1] ** 2) * vv ** 2 +
                                                (sigma[2] ** 2) * ww ** 2)) + 1
        midpoint = np.ceil(np.array(uu.shape) / 2.0)
        maxval = ndimage.maximum(gfilter[midpoint[0] - 10:midpoint[0] + 10,
                                         midpoint[1] - 10:midpoint[1] + 10,
                                         midpoint[2] - 10:midpoint[2] + 10])
        gfilter = gfilter / maxval
    return gfilter


def gaussian_fourierkernel_quarter(uu, vv, ww, sigma):
    """ Calculate one eighth of gaussian FT matrix
    """
    siz = np.floor(np.array(uu.shape) / 2 + 1).astype(int)
    tic()
    # ravel faster than flatten or reshape(,(N*M,))
    u_ = uu[:siz[0], :siz[1], :siz[2]].ravel()
    v_ = vv[:siz[0], :siz[1], :siz[2]].ravel()
    w_ = ww[:siz[0], :siz[1], :siz[2]].ravel()
    z_ = np.exp(-2 * (np.pi ** 2) * (u_ ** 2 + v_ ** 2 + w_ ** 2) * (sigma**2))
    gfilter = np.empty_like(uu)
    gfilter[:siz[0], :siz[1], :siz[2]] = np.reshape(z_, siz)
    gfilter[siz[0]:, :siz[1], :siz[2]] = np.roll(
        np.reshape(z_, siz), siz[0], axis=0)[1:-1, :, :]  # second quad
    gfilter[:siz[0], siz[1]:, :siz[2]] = np.roll(
        np.reshape(z_, siz), siz[1], axis=1)[:, 1:-1, :]  # third quad
    #gfilter[:siz[0],:siz[1],siz[2]:] = np.roll(np.reshape(z_,siz),siz[2],2)[:,:,1:-2]
    gfilter[siz[0]:, siz[1]:, :siz[2]] = np.roll(
        np.roll(np.reshape(z_, siz), siz[0], 0), siz[1], 1)[1:-1, 1:-1, :]

    gfilter[siz[0]:, :siz[1], :siz[2]] = gfilter[1:-1, :, :]

    gfilter[:, :, siz[2]:] = np.roll(gfilter[:, :, 1:siz[2] - 1])
    np.roll(gfilter[:, :, 1:siz[2] - 1], siz[2] - 2, axis=2)
    gfilter /= z_[-1]  # inplace division
    toc()


def gaussian_fourierkernel_quarter_v2(uu, vv, ww, sigma)
    siz = np.floor(np.array(uu.shape) / 2 + 1).astype(int)
    gfilter = np.empty_like(uu)
    alpha = np.exp(-2 * (np.pi ** 2))
    for ii in xrange(0, siz[0]):
        for jj in xrange(0, siz[1]):
            for kk in xrange(0, siz[2]):
                gfilter[ii, jj, kk] = np.exp(((sigma[0] ** 2) * uu[ii, jj, kk] ** 2) +
                                             ((sigma[1] ** 2) * vv[ii, jj, kk] ** 2) +
                                             ((sigma[2] ** 2) * ww[ii, jj, kk] ** 2))

    # Copy second quadrant
    gfilter[siz[0]:, :siz[1], :siz[2]] = gfilter[
        siz[0] - 2:0:-1, :siz[1], :siz[2]]
    # Copy third and fourth quadrant
    gfilter[:siz[0], siz[1]:, :siz[2]] = gfilter[
        :siz[0], siz[1] - 2:0:-1, :siz[2]]
    # Copy fifth to eighth quadrant
    gfilter[:, :, siz[2]:] = gfilter[:, :, siz[2] - 2:0:-1]
    maxval = gfilter[siz[0], siz[1], siz[2]]
    gfilter *= alpha / maxval


def gaussian_fourierkernel_quarter_v3(uu, vv, ww, sigma)
    siz = np.floor(np.array(uu.shape) / 2 + 1).astype(int)
    gfilter = np.empty_like(uu)
    alpha = np.exp(-2 * (np.pi ** 2))
    for ii in xrange(0, siz[0]):
        for jj in xrange(0, siz[1]):
            for kk in xrange(0, siz[2]):
                gfilter[ii, jj, kk] = np.exp(((sigma[0] ** 2) * uu[ii, jj, kk] ** 2) +
                                             ((sigma[1] ** 2) * vv[ii, jj, kk] ** 2) +
                                             ((sigma[2] ** 2) * ww[ii, jj, kk] ** 2))

    # Copy second quadrant
    m = -(siz[0] - 3) / (siz[0] - 2)
    c = 1 + uu.shape[0] * (siz[0] - 3) / (siz[0] - 2)
    for ii in xrange(siz[0], uu.shape[0]):
        for jj in xrange(0, siz[1]):
            for kk in xrange(0, siz[2]):
                gfilter[ii, jj, kk] = gfilter[m * ii + c, jj, kk]

    # Copy third and fourth quadrant
    m = -(siz[1] - 3) / (siz[1] - 2)
    c = 1 + uu.shape[1] * (siz[1] - 3) / (siz[1] - 2)
    for ii in xrange(0, uu.shape[0]):
        for jj in xrange(siz[1], uu.shape[1]):
            for kk in xrange(0, siz[2]):
                gfilter[ii, jj, kk] = gfilter[ii, m * jj + c, kk]
    # Copy fifth to eighth quadrant
    m = -(siz[2] - 3) / (siz[2] - 2)
    c = 1 + uu.shape[2] * (siz[2] - 3) / (siz[2] - 2)
    for ii in xrange(0, uu.shape[0]):
        for jj in xrange(0, uu.shape[1]):
            for kk in xrange(siz[2], uu.shape[2]):
                gfilter[ii, jj, kk] = gfilter[ii, jj, m * kk + c]
    return gfilter


def gaussian_sourcemodule(uu, vv, ww, sigma):

    import pycuda.driver as drv
    import pycuda.tools
    import pycuda.autoinit
    import numpy
    from pycuda.compiler import SourceModule
    import pycuda.gpuarray as gpuarray
    ######################
    # SourceModele SECTION
    # We write the C code and the indexing and we have lots of control
    # create two timers so we can speed-test each approach
    start = drv.Event()
    end = drv.Event()
    mod = SourceModule("""
            __global__ void gpugauss(float *dest, float *u, float *v, float *w, float sigma)
                {
                   const int i = blockDim.x*blockIdx.x + threadIdx.x;
                   dest[i] = exp( -2*pi*(sigma^2) *(u[i]^2 + v[i]^2 +w[i]^2));
                }
                       """)

    gpugauss = mod.get_function("gpugauss")
    siz = np.floor(np.array(uu.shape) / 2 + 1)
    zz = uu[:siz[1], :siz[2], :siz[3]].flatten()
    u_gpu = gpuarray.to_gpu(uu[:siz[1], :siz[2], :siz[3]].ravel())
    v_gpu = gpuarray.to_gpu(vv[:siz[1], :siz[2], :siz[3]].ravel())
    w_gpu = gpuarray.to_gpu(ww[:siz[1], :siz[2], :siz[3]].ravel())

    # create a destination array that will receive the result
    dest = numpy.zeros_like(zz)
    blocks = 4
    block_size = (zz.size) / blocks
    start.record()  # start timing
    gpugauss(drv.Out(dest), drv.In(u_gpu), drv.In(v_gpu), drv.In(
        w_gpu), grid=(blocks, 1), block=(block_size, 1, 1))
    end.record()  # end timing
    # calculate the run length
    end.synchronize()
    secs = start.time_till(end) * 1e-3
    print "SourceModule time and first three results:"
    print "%fs, %s" % (secs, str(dest[:3]))
    gfilter[:siz[0], :siz[1], :siz[2]] = np.reshape(dest.get(), siz)
    # Copy second quadrant
    gfilter[siz[0]:, :siz[1], :siz[2]] = gfilter[siz[0] - 2:-1:1, :, :]
    # Copy third and fourth quadrant
    gfilter[:, siz[1]:, :siz[2]] = gfilter[:, siz[1] - 2:-1:1, :]
    # Copy fifth to eighth quadrant
    gfilter[:, :, siz[2]:] = gfilter[:, :, siz[2] - 2:-1:1]
    maxval = gfilter[siz[0], siz[1], siz[2]]
    gfilter /= maxval
    return gfilter


def fouriergauss(siz, sigma):
    """
    Gaussian operator in Fourier domain is another Gaussian :
      g(x,y,z)=(1/sqrt(2*pi).sigma).exp(-(x^2+y^2+z^2)/2sigma^2)
          => A.(exp(-2*pi*pi*(u^2 + v^2 + w^2)*(sigma^2))

    The sigma should be relative to the voxel spacing
    """
    (uu, vv, ww) = fouriercoords(siz)
    return gaussian_fourierkernel(uu, vv, ww, sigma)


def cplxfouriergauss(siz, sigma):
    """
    Complex Gaussian in Fourier domain :
      g(x,y,z)=(1+i)(1/sqrt(2*pi).sigma).exp(-(x^2)/2sigma^2)
    .. math::
    $\mathcal{F_x}[f(x)](\omega)  = A*(1+i)*exp(-((w^2)*(sigma^2))/2 + i*mu*w)$
    $A=1/[sqrt(2*pi/sigma^2)*sigma]$
    mu is zero, so the real and imag components are:
      [exp(-((w^2)*(sigma^2))/2)] / [sqrt(2*pi/sigma^2)*sigma]


    The sigma should be relative to the voxel spacing.
    """
    (uu, vv, ww) = fouriercoords(siz)
    gfilter = gaussian_fourierkernel(uu, vv, ww, sigma)
    return gfilter + 1j * gfilter


def fourierlaplace(siz):
    """
    Laplacian operator in Fourier domain is very simple:
      D^2 g(x,y,z) => -(4pi^2) (u^2 + v^2 + w^2) G(u,v,w)
    """
    (uu, vv, ww) = fouriercoords(siz)
    # / (siz[0] * siz[1] * siz[2])
    return -(4 * np.pi * np.pi) * (uu * uu + vv * vv + ww * ww)


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
    return np.abs(fftshift(fftn(Kfilter)))


def fouriergausssubband15(siz, sigma):
    """ Subband15 Gaussian filter

    Based on 3D filtering paper:
    Max W. K. Law and Albert C. S. Chung, "Efficient Implementation for
    Spherical Flux Computation and Its Application to Vascular Segmentation,
       IEEE Transactions on Image Processing, 2009, Volume 18(3), 596V612

    http://www.cse.ust.hk/~maxlawwk/

    """
    (uu, vv, ww) = fouriercoords(siz)
    # Original Gaussian kernel
    gfilter = gaussian_fourierkernel(uu, vv, ww, sigma)
    # Subband_1.5 frequency oversampling component.

    gfilter = gfilter + gaussian_fourierkernel(uu + 1, vv, ww, sigma)
    gfilter = gfilter + gaussian_fourierkernel(uu - 1, vv, ww, sigma)
    gfilter = gfilter + gaussian_fourierkernel(uu, vv + 1, ww, sigma)
    gfilter = gfilter + gaussian_fourierkernel(uu, vv - 1, ww, sigma)
    gfilter = gfilter + gaussian_fourierkernel(uu, vv, ww + 1, sigma)
    gfilter = gfilter + gaussian_fourierkernel(uu, vv, ww - 1, sigma)
    # Normalization improves accuracy when sigma is small (e.g. sigma < 0.8
    # voxel length)
    gfilter = gfilter / ndimage.maximum(gfilter)

# End of Subband_1.5 frequency oversampling component
    return gfilter


def fouriergauss2(siz, voxmm, sigma):
    """
    Complex Gaussian in Fourier domain :
      g(x,y,z)=(1+i)(1/sqrt(2*pi).sigma).exp(-(x^2)/2sigma^2)
    .. math::
    $\mathcal{F_x}[f(x)](\omega)  => A*(1+i)*exp(-((w^2)*(sigma^2))/2 + i*mu*w)$
    $A=1/[sqrt(2*pi/sigma^2)*sigma]$
    mu is zero, so the real and imag components are:
      [exp(-((w^2)*(sigma^2))/2)] / [sqrt(2*pi/sigma^2)*sigma]


    The sigma should be relative to the voxel spacing.
    """
    (uu, vv, ww) = fouriercoords(siz)
    if not hasattr(sigma, "__len__"):
        # if type(sigma) is float or type(sigma) is np.float64:
        factor = 1 / (np.sqrt(2 * np.pi / sigma ** 2) * sigma)
        component = factor * (np.expm1(-0.5 * (uu ** 2 + vv ** 2 + ww ** 2) *
                                       (sigma ** 2)) + 1)
    else:
        factor = 1 / (np.sqrt(2 * np.pi / ((sigma[0] * sigma[0]) +
                                           (sigma[1] * sigma[1]) +
                                           (sigma[2] * sigma[2]))) *
                      np.prod(sigma))
        component = factor * (np.expm1(-0.5 * ((sigma[0] * sigma[0] * uu * uu) +
                                               (sigma[1] * sigma[1] * vv * vv) +
                                               (sigma[2] * sigma[2] * ww * ww))) + 1)
    return component


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
                  * (2 * sigma * sigma)) + 1
    del uu, vv, ww

    kspHG = ksp * HG
    del HG
    fsmoothed = abs(fftshift(ifftn(ifftshift(kspHG))))
    fsmoothed = fsmoothed / ndimage.mean(fsmoothed)
    return fsmoothed
# end inhomogeneouscorrection


def kspacegaussian_filter(ksp, sigma=np.sqrt(0.5), n_=-1, axis_=-1):
    """KSPACEFILTER gaussian filter of complex 3D image

    filtered_magnitude = kspacefilter(realimg, imagimg)


    scipy.ndimage.fourier.fourier_gaussian

    scipy.ndimage.fourier.fourier_gaussian(input, sigma, n=-1, axis=-1,
    output=None)[source]

    Multi-dimensional Gaussian fourier filter.

    The array is multiplied with the fourier transform of a Gaussian kernel.

    Parameters:
    input : array_like
    The input array.
    sigma : float or sequence

    The sigma of the Gaussian kernel. If a float, sigma is the same
    for all axes. If a sequence, sigma has to contain one value for
    each axis.

    n : int, optional

    If n is negative (default), then the input is assumed to be the
    result of a complex fft. If n is larger than or equal to zero, the
    input is assumed to be the result of a real fft, and n gives the
    length of the array before transformation along the real transform
    direction.

    axis : int, optional
    The axis of the real transform.
    output : ndarray, optional

    If given, the result of filtering the input is placed in this array. None
    is returned in this case.

    Returns:
    fourier_gaussian : ndarray or None
    The filtered input. If output is given as a parameter, None is returned.

    """

    print "Complex Gaussian filter sigma ", sigma
    if ksp.ndim == 3:
        out_ksp = ndimage.fourier.fourier_gaussian(
            ksp, sigma, n=n_, axis=axis_)
    else:
        out_ksp = np.empty_like(ksp, dtype=np.complex64)

        for echo in xrange(0, ksp.shape[4]):
            for n in xrange(0, ksp.shape[3]):
                out_ksp[:, :, :, n, echo] = ndimage.fourier.fourier_gaussian(
                    ksp, sigma, n=n_, axis=axis_)
    return out_ksp
# end kspacegaussian_filter


def kspacegaussian_filter2(ksp, sigma_=None):
    """
    Apply Gaussian filter in Fourier domain to kspace data
    """
    siz = ksp.shape[0:3]
    sigma = np.ones(3)
    if 'sigma_' not in locals():
        sigma = np.array(siz) / (4 * np.sqrt(2 * np.log(2)))
    elif not hasattr(sigma_, "__len__"):
        sigma = np.ones(3) * sigma_
    else:
        sigma = sigma_.copy()
    Fgauss = fouriergauss(siz, sigma)
    out_ksp = np.empty_like(ksp, dtype=np.complex64)
    print "Complex Gaussian filter sigma ", sigma
    if ksp.ndim == 3:
        out_ksp.real = ksp.real * Fgauss
        out_ksp.imag = ksp.imag * Fgauss
    else:
        for echo in xrange(0, ksp.shape[4]):
            for n in xrange(0, ksp.shape[3]):
                out_ksp[:, :, :, n, echo].real = ksp[
                    :, :, :, n, echo].real * Fgauss
                out_ksp[:, :, :, n, echo].imag = ksp[
                    :, :, :, n, echo].imag * Fgauss
    return out_ksp
# end kspacegaussian_filter2


def kspacecplxgaussian_filter(ksp, sigma_=None):
    """
    Apply Gaussian filter in Fourier domain to kspace real and imag data
    """
    siz = ksp.shape[0:3]
    sigma = np.ones(3)
    if 'sigma_' not in locals():
        sigma = np.array(siz) / (4 * np.sqrt(2 * np.log(2)))
    elif not hasattr(sigma_, "__len__"):
        sigma = np.ones(3) * sigma_
    else:
        sigma = sigma_.copy()
    Fgauss = cplxfouriergauss(siz, sigma)
    out_ksp = np.empty_like(ksp, dtype=np.complex64)
    print "Complex Gaussian filter sigma ", sigma
    if ksp.ndim == 3:
        out_ksp = ksp * Fgauss
    else:
        for echo in xrange(0, ksp.shape[4]):
            for n in xrange(0, ksp.shape[3]):
                out_ksp[:, :, :, n, echo] = ksp[:, :, :, n, echo] * Fgauss
    return out_ksp
# end kspacecplxgaussian_filter


def kspacelaplacegaussian_filter(ksp, sigma_=None):
    """
    Apply Laplace of Gaussian filter in Fourier domain to kspace data
    """
    siz = ksp.shape[0:3]
    sigma = np.ones(3)
    if 'sigma_' not in locals():
        sigma = np.array(siz) / (4 * np.sqrt(2 * np.log(2)))
    if not hasattr(sigma_, "__len__"):
        sigma = np.ones(3) * sigma_
    else:
        sigma = sigma_.copy()
    Flaplace = fourierlaplace(siz)
    Fgauss = fouriergauss(siz, sigma)
    out_ksp = np.empty_like(ksp, dtype=np.complex64)
    print "Complex Laplace Gaussian filter sigma ", sigma
    if ksp.ndim == 3:
        out_ksp.real = (ksp.real * Fgauss) * Flaplace
        out_ksp.imag = (ksp.imag * Fgauss) * Flaplace
    else:
        for echo in xrange(0, ksp.shape[4]):
            for n in xrange(0, ksp.shape[3]):
                out_ksp[:, :, :, n, echo].real = ksp[
                    :, :, :, n, echo].real * Flaplace * Fgauss
                out_ksp[:, :, :, n, echo].imag = ksp[
                    :, :, :, n, echo].imag * Flaplace * Fgauss
    return out_ksp
# end kspacelaplacegaussian_filter


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


def kspaceepanechnikov_filter2(ksp, bdwidth_=None):
    """
    Apply Epanechnikov filter in Fourier domain to kspace real and imag data
    """

    siz = ksp.shape[0:3]
    if 'bdwidth_' not in locals():
        bdwidth = np.sqrt(7) * np.ones(3) * siz / (4 * np.sqrt(2 * np.log(2)))
    else:
        if not hasattr(bdwidth_, "__len__"):
            bdwidth = np.ones(3) * bdwidth_
        else:
            bdwidth = bdwidth_.copy()
    Fepanechnikov = fourierepanechnikov(siz, bdwidth)
    out_ksp = np.empty_like(ksp, dtype=np.complex64)
    CmplxEpan = np.empty_like(ksp, dtype=np.complex64)
    CmplxEpan.real = Fepanechnikov
    CmplxEpan.imag = Fepanechnikov
    print "Complex Epanechnikov filter bandwidth ", bdwidth
    if ksp.ndim == 3:
        out_ksp = CmplxEpan * ksp
    else:
        for echo in xrange(0, ksp.shape[4]):
            for n in xrange(0, ksp.shape[3]):
                out_ksp[:, :, :, n, echo] = ksp[:, :, :, n, echo] * CmplxEpan
    return out_ksp
# end kspaceepanechnikov_filter


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
    image_filtered = fftshift(ifftn(ifftshift(ksplarge)))
    print "Saving Double res image: " + basename + " filtered"
    save_nifti(np.abs(image_filtered), basename + '_large')


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

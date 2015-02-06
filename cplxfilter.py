#!/usr/bin/env python
"""cplxfilter complex filtering methods for 3D MR image data

  Methods include:
  cplxgaussian_filter, cplxgaussian_laplace, cplxepanechnikov_filter
  cplxmedian_filter, cplxwiener_filter


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

# import numpy
import numpy as np
import scipy
if scipy.__version__[2] == 7:
    scipy.pkgload('signal')
    scipy.pkgload('ndimage')
    # scipy.pkgload('fftpack')
else:
    # from scipy.fftpack import fftn, ifftn, fftshift, ifftshift
    from scipy import ndimage
    from scipy import signal


def cplxgaussian_filter(real_input, imag_input, sigma=0.707, order_=0,
                        mode_='nearest', cval_=0.0):
    """CPLXFILTER gaussian filter of complex 3D image
       ndimage.filters.gaussian_filter is used to smooth real and imag
       components filtered_magnitude = cplxfilter(realimg, imagimg)

    :param RE, IM: real and imag NDimage
    :param sigma: [optional] optimal sigma = 1/sqrt(2).
    Standard deviation for Gaussian kernel. The standard deviations of
    the Gaussian filter are given for each axis as a sequence, or as a
    single number, in which case it is equal for all axes.

    :param order: [optional] ]{0, 1, 2, 3} or sequence from same set, optional.
    The order of the filter along each axis is given as a sequence of
    integers, or as a single number. An order of 0 corresponds to
    convolution with a Gaussian kernel. An order of 1, 2, or 3
    corresponds to convolution with the first, second or third
    derivatives of a Gaussian. Higher order derivatives are not
    implemented

    :param mode: [optional] ]{'reflect', 'constant', 'nearest',
    'mirror', 'wrap'}, optional. Default is 'reflect'
    The mode parameter determines how the array borders are handled,
    where cval is the value when mode is equal to 'constant'.

    :param cval: [optional]
    Value to fill past edges of input if mode is 'constant'. Default
    is 0.0

    :return filtered_image:  complex 3D array of gaussian filtered image,


    scipy.ndimage.filters.gaussian_filter

scipy.ndimage.filters.gaussian_filter(input, sigma, order=0, output=None,
mode='reflect', cval=0.0, truncate=4.0)[source] Multidimensional Gaussian
filter.

Parameters:
input : array_like
Input array to filter.
sigma : scalar or sequence of scalar:

Standard deviation for Gaussian kernel. The standard deviations of the Gaussian
filter are given for each axis as a sequence, or as a single number, in which
case it is equal for all axes.

order : {0, 1, 2, 3} or sequence from same set, optional

The order of the filter along each axis is given as a sequence of integers, or
as a single number. An order of 0 corresponds to convolution with a Gaussian
kernel. An order of 1, 2, or 3 corresponds to convolution with the first,
second or third derivatives of a Gaussian. Higher order derivatives are not
implemented

output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional

The mode parameter determines how the array borders are handled, where cval is
the value when mode is equal to 'constant'. Default is 'reflect' cval : scalar,
optional Value to fill past edges of input if mode is 'constant'. Default is
0.0 truncate : float Truncate the filter at this many standard
deviations. Default is 4.0.  Returns: gaussian_filter : ndarray Returned array
of same shape as input.  Notes

The multidimensional filter is implemented as a sequence of one-dimensional
convolution filters. The intermediate arrays are stored in the same data type
as the output. Therefore, for output types with a limited precision, the
results may be imprecise because intermediate results may be stored with
insufficient precision.

    """
    # imgmag = load_untouch_nii(file1)
    # imph = load_untouch_nii(file2)
    # imcplx=(imgmag.img).*exp((1i).*(imph.img))
    # c = fftshift(ifftn(fftshift(imcplx)))
    # c = flipdim(flipdim(flipdim(c, 1), 2), 3)
    # f = fspecial3('gaussian', [3 3 3], 1) %/sqrt(2))
    # compleximg = complex(imfilter(real(c), f), imfilter(imag(c), f))

    print "Complex Gaussian filter sigma ", sigma, " order ", order_, \
        " mode ", mode_
    if real_input.ndim == 3:
        real_img = ndimage.filters.gaussian_filter(real_input, sigma,
                                                   order=order_, mode=mode_,
                                                   cval=cval_)
        # truncate=4.0) truncate not supported in scipy 0.10
        imag_img = ndimage.filters.gaussian_filter(imag_input, sigma,
                                                   order=order_, mode=mode_,
                                                   cval=cval_)
        # truncate=4.0)
    else:
        real_img = np.empty_like(real_input, dtype=np.float32)
        imag_img = np.empty_like(real_input, dtype=np.float32)
        for echo in xrange(0, real_input.shape[4]):
            for acq in xrange(0, real_input.shape[3]):
                real_img[:, :, :, acq, echo] = ndimage.filters.gaussian_filter(
                    real_input[:, :, :, acq, echo],
                    sigma, order=order_, mode=mode_, cval=cval_)
                imag_img[:, :, :, acq, echo] = ndimage.filters.gaussian_filter(
                    imag_input[:, :, :, acq, echo],
                    sigma, order=order_, mode=mode_, cval=cval_)
    filtered_image = np.empty_like(real_input, dtype=np.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    return filtered_image
# end cplxfilter


def cplx2dgaussian_filter(real_input, imag_input, sigma=0.707, order_=0,
                          mode_='nearest', cval_=0.0):
    """CPLX2DFILTER gaussian filter of complex 3D image
     ndimage.filters.gaussian_filter is used to smooth real and imag components
       filtered_magnitude = cplxfilter(realimg, imagimg)

    :param RE, IM: real and imag NDimage

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

    :param mode: [optional] ]{'reflect', 'constant', 'nearest', 'mirror',
    'wrap'},
    optional. The mode parameter determines how the array borders are
    handled, where cval is the value when mode is equal to
    'constant'. Default is 'reflect'

    :param cval: [optional] Value to fill past edges of input if mode is
    'constant'. Default is 0.0

    :return filtered_image:  complex 3D array of gaussian filtered image,

    """
    print "Complex 2D Gaussian filter sigma ", sigma, " order ", order_, \
        " mode ", mode_
    real_img = np.empty_like(real_input, dtype=np.float32)
    imag_img = np.empty_like(real_input, dtype=np.float32)
    if real_input.ndim == 3:
        for islice in xrange(0, real_input.shape[2]):
            real_img[:, :, islice] = ndimage.filters.gaussian_filter(
                real_input[:, :, islice], sigma,
                order=order_, mode=mode_,
                cval=cval_)
            # truncate=4.0) truncate not supported in scipy 0.10
            imag_img[:, :, islice] = ndimage.filters.gaussian_filter(
                imag_input[:, :, islice], sigma,
                order=order_, mode=mode_,
                cval=cval_)
            # truncate=4.0)
    else:
        for echo in xrange(0, real_input.shape[4]):
            for acq in xrange(0, real_input.shape[3]):
                for islice in xrange(0, real_input.shape[2]):
                    real_img[:, :, islice, acq, echo] = ndimage.filters.gaussian_filter(
                        real_input[:, :, islice, acq, echo], sigma,
                        order=order_,
                        mode=mode_,
                        cval=cval_)
                    imag_img[:, :, islice, acq, echo] = ndimage.filters.gaussian_filter(
                        imag_input[:, :, islice, acq, echo], sigma,
                        order=order_,
                        mode=mode_,
                        cval=cval_)
    filtered_image = np.empty_like(real_input, dtype=np.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    return filtered_image
# end cplx2dfilter_gaussian


def cplxgaussian_laplace(real_input, imag_input, sigma=0.707, mode_='reflect',
                         cval_=0.0):
    """CPLXFILTER gaussian laplace filter of complex 3D image
       ndimage.filters.gaussian_laplace is used to process real and imag
     components and estimate blobs
    filtered_magnitude = cplxfilter(realimg, imagimg)

    :param sigma:  optimal sigma = 1/sqrt(2).  Standard
    deviation for kernel. The standard deviations of the Gaussian Laplace
    filter are given for each axis as a sequence, or as a single number, in
    which case it is equal for all axes.

    :param mode: {'reflect', 'constant', 'nearest', 'mirror', 'wrap'},
    optional. The mode parameter determines how the array borders are handled,
    where cval is the value when mode is equal to 'constant'. Default is
    'reflect'

    :param cval: [optional] Value to fill past edges of input if mode is
    'constant'. Default is 0.0

    :return filtered_image: complex 3D array of gaussian laplace filtered
    image,
   scipy.ndimage.filters.gaussian_laplace

scipy.ndimage.filters.gaussian_laplace(input, sigma, output=None,
mode='reflect', cval=0.0, **kwargs)[source] Multidimensional Laplace filter
using gaussian second derivatives.

Parameters:
input : array_like
Input array to filter.
sigma : scalar or sequence of scalars

    The standard deviations of the Gaussian filter are given for each axis as a
    sequence, or as a single number, in which case it is equal for all axes.

    output : array, optional

    The output parameter passes an array in which to store the filter output.

    mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional The
mode parameter determines how the array borders are handled, where cval is the
value when mode is equal to 'constant'. Default is 'reflect'

    cval : scalar, optional
Value to fill past edges of input if mode is 'constant'. Default is 0.0 Extra
keyword arguments will be passed to gaussian_filter().

    """

    print "Complex Gaussian_laplace filter sigma ", sigma, " mode ", mode_
    if real_input.ndim == 3:
        real_img = ndimage.filters.gaussian_laplace(real_input, sigma,
                                                    mode=mode_, cval=cval_)
        imag_img = ndimage.filters.gaussian_laplace(imag_input, sigma,
                                                    mode=mode_, cval=cval_)
    else:
        real_img = np.empty_like(real_input, dtype=np.float32)
        imag_img = np.empty_like(real_input, dtype=np.float32)
        for echo in xrange(0, real_input.shape[4]):
            for acq in xrange(0, real_input.shape[3]):
                real_img[:, :, :, acq, echo] = ndimage.filters.gaussian_laplace(
                    real_input[:, :, :, acq, echo],
                    sigma, mode=mode_, cval=cval_)
                imag_img[:, :, :, acq, echo] = ndimage.filters.gaussian_laplace(
                    imag_input[:, :, :, acq, echo],
                    sigma, mode=mode_, cval=cval_)
    filtered_image = np.empty_like(real_input, dtype=np.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    return filtered_image
# end cplxgaussian_laplace


def cplxlaplacian_filter(real_input, imag_input, mode_='reflect', cval_=0.0):
    """CPLXFILTER apply scipy's laplace filter to complex 3D image

    filtered_image = cplxlaplacian_filter(realimg, imagimg)

    :param mode: {'reflect', 'constant', 'nearest', 'mirror', 'wrap'},
    optional. The mode parameter determines how the array borders are handled,
    where cval is the value when mode is equal to 'constant'. Default is
    'reflect'
       :param cval: [optional] Value to fill past edges of input if mode is
    'constant'. Default is 0.0

    :return filtered_image: complex 3D array of gaussian laplace filtered
    image,
    scipy.ndimage.filters.laplace

    scipy.ndimage.filters.laplace(input, output=None, mode='reflect',
    cval=0.0)[source] N-dimensional Laplace filter based on approximate second
    derivatives.

    Parameters:
    input : array_like
    Input array to filter.
    output : array, optional
    The output parameter passes an array in which to store the filter output.
    mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional

    The mode parameter determines how the array borders are handled, where cval
    is the value when mode is equal to 'constant'. Default is 'reflect' cval :
    scalar, optional Value to fill past edges of input if mode is
    'constant'. Default is 0.0

    """
    print "Complex Gaussian_laplace filter  ", " mode ", mode_
    if real_input.ndim == 3:
        real_img = ndimage.filters.laplace(real_input, mode=mode_, cval=cval_)
        imag_img = ndimage.filters.laplace(imag_input, mode=mode_, cval=cval_)
    else:
        real_img = np.empty_like(real_input, dtype=np.float32)
        imag_img = np.empty_like(real_input, dtype=np.float32)
        for echo in xrange(0, real_input.shape[4]):
            for acq in xrange(0, real_input.shape[3]):
                real_img[:, :, :, acq, echo] = ndimage.filters.laplace(
                    real_input[:, :, :, acq, echo],
                    mode=mode_, cval=cval_)
                imag_img[:, :, :, acq, echo] = ndimage.filters.laplace(
                    imag_input[:, :, :, acq, echo],
                    mode=mode_, cval=cval_)
    filtered_image = np.empty_like(real_input, dtype=np.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    return filtered_image
# end cplxlaplace


def cplxmedian_filter(real_input, imag_input, size_=5, mode_='reflect'):
    """
    scipy.ndimage.filters.median_filter(input, size=None, footprint=None,
     output=None, mode='reflect', cval=0.0, origin=0) not used footprint_=[5,
     5, 5], output_=None, mode_='reflect', cval_=0.0, origin_=0

    scipy.ndimage.filters.median_filter(input, size=None, footprint=None,
output=None, mode='reflect', cval=0.0, origin=0)[source]

Calculates a multidimensional median filter.

Parameters:
input : array_like
Input array to filter.
size : scalar or tuple, optional
See footprint, below
footprint : array, optional

Either size or footprint must be defined. size gives the shape that is taken
from the input array, at every element position, to define the input to the
filter function. footprint is a boolean array that specifies (implicitly) a
shape, but also which of the elements within this shape will get passed to the
filter function. Thus size=(n, m) is equivalent to footprint = np.ones((n,
m)). We adjust size to the number of dimensions of the input array, so that, if
the input array is shape (10, 10, 10), and size is 2, then the actual size used
is (2, 2, 2).

output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional

The mode parameter determines how the array borders are handled, where cval is
the value when mode is equal to 'constant'. Default is 'reflect'

cval : scalar,
optional Value to fill past edges of input if mode is 'constant'. Default is
0.0

origin : scalar, optional The origin parameter controls the placement of
the filter. Default 0.0.

Returns:
median_filter : ndarray
Return of same shape as input.

    """
    filtered_image = np.empty(real_input.shape, dtype=np.complex64)
    print "Complex Median filter window size(s)", size_
    if not hasattr(size_, "__len__"):
        footprint_ = np.ones((size_, size_, size_))
    else:
        footprint_ = np.ones(size_[0:3])
    if real_input.ndim == 3:
        print "Complex Median filter"
        real_img = ndimage.filters.median_filter(real_input,
                                                 footprint=footprint_,
                                                 mode=mode_)
        imag_img = ndimage.filters.median_filter(imag_input,
                                                 footprint=footprint_,
                                                 mode=mode_)
    else:
        real_img = np.empty_like(real_input, dtype=np.float32)
        imag_img = np.empty_like(real_input, dtype=np.float32)
        for echo in xrange(0, real_input.shape[4]):
            for acq in xrange(0, real_input.shape[3]):
                real_img[:, :, :, acq, echo] = ndimage.filters.median_filter(
                    real_input[:, :, :, acq, echo],
                    footprint=footprint_, mode=mode_)
                imag_img[:, :, :, acq, echo] = ndimage.filters.median_filter(
                    imag_input[:, :, :, acq, echo],
                    footprint=footprint_, mode=mode_)
    filtered_image = np.empty_like(real_input, dtype=np.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    return filtered_image
# end cplxmedian_filter


def cplxwiener_filter(real_input, imag_input, mysize_=5, noise_=None):
    """cplxwiener_filter Implementation of Wiener filter on complex image

    scipy.signal.wiener(im, mysize=None, noise=None)[source]
Perform a Wiener filter on an N-dimensional array.

The Wiener filter is a simple deblurring filter for denoising images. This is
not the Wiener filter commonly described in image reconstruction problems but
instead it is a simple, local-mean filter.

Apply a Wiener filter to the N-dimensional array im.

Parameters: im : ndarray An N-dimensional array.  mysize : int or arraylike,
optional A scalar or an N-length list giving the size of the Wiener filter
window in each dimension. Elements of mysize should be odd. If mysize is a
scalar, then this scalar is used as the size in each dimension.  noise : float,
optional The noise-power to use. If None, then noise is estimated as the
average of the local variance of the input.  Returns: out : ndarray Wiener
filtered result with the same shape as im.

    """
    # scipy.signal.wiener(im, mysize=None, noise=None)
    # ,(size_, size_, size_)
    print "Complex Wiener filter window size ", mysize_, " noise ", noise_
    if not noise_:
        noise_ = ndimage.standard_deviation(real_input)
    if not hasattr(mysize_, "__len__"):
        filter_size = np.array(mysize_)
    else:
        filter_size = mysize_
    if real_input.ndim == 3:
        real_img = signal.wiener(real_input, mysize=filter_size, noise=noise_)
        imag_img = signal.wiener(imag_input, mysize=filter_size, noise=noise_)
    else:
        real_img = np.empty_like(real_input)
        imag_img = np.empty_like(real_input)
        for echo in xrange(0, real_input.shape[4]):
            for acq in xrange(0, real_input.shape[3]):
                real_img[:, :, :, acq, echo] = signal.wiener(real_input[:, :, :, acq, echo],
                                                             mysize=filter_size,
                                                             noise=noise_)
                imag_img[:, :, :, acq, echo] = signal.wiener(imag_input[:, :, :, acq, echo],
                                                             mysize=filter_size,
                                                             noise=noise_)

    filtered_image = np.empty_like(real_input, dtype=np.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    return filtered_image
# end cplxwiener_filter


def epanechnikov(filtersiz, bandwidth, order=0):
    """
    Epanechnikov filter  3/4 * (1-|u|^2), -1 <= u <= 1
    u = x/sigma

    Second order: k1(x) = 3/160 (2-x)^3 (x^2 + 6x + 4)
       Fourth order: -15/8 x^2 + 9/8   (norm 1/sqrt(5))
       Kernel Equation R(k) 4(k) eff(k)
Epanechnikov k4;1(u) = \frac{15}{8}(1-\frac{7}{3}u^2)* k1(u))
k4;1(x) = 5/2048 (2-x)^3 (7x^6 + 42x^5 + 48x^4 -160x^3 - 144x^2 + 96x + 64)
    Sixth order:
    Epanechnikov k6;1(u) = \frac{175}{64}(1-6u^2 + \frac{33}{5}u^4)* k1(u))

    Table 4: Rule of Thumb Constant
11  Kernel \nu = 2 \nu = 4 \nu = 6
Epanechnikov 2.34 3.03 3.53
Biweight 2.78 3.39 3.84
Triweight 3.15 3.72 4.13
Gaussian 1.06 1.08 1.08

/* Second Order Epanechnikov */
    /* Note that return value is preset to 0 so no ifelse necessary */
if (z_squared < 5.0) return_value =
  (double)(0.33541019662496845446-0.067082039324993690892*z_squared);
/* Fourth Order Epanechnikov */
if (z_squared < 5.0) return_value =
   (double)(0.008385254916*(-15.0 + 7.0*z_squared)*(-5.0 + z_squared));
/* Sixth Order Epanechnikov */
if (z_squared < 5.0) return_value =
   (double)(0.33541019662496845446*(2.734375 + z_squared*(-3.28125
    + 0.721875*z_squared))*(1.0-0.2*z_squared));
/* Eighth Order Epanechnikov */
if (z_squared < 5.0) return_value =
    (double)(0.33541019662496845446*(3.5888671875 + z_squared*(-7.8955078125
    + z_squared*(4.1056640625-0.5865234375*z_squared)))*(1.0-0.2*z_squared));

Nonparametric Econometrics: Theory and Practice
 By Qi Li, Jeffrey Scott Racine
    https://books.google.com.au/books?id = Zsa7ofamTIUC&pg = PT66&lpg =
    PT66&dq = fourth%2Border%2Bepanechnikov&source = bl&ots = z7PcnajR1z&sig
    = bx0-HS0AUYSx7j0gFSkQbu-5buQ&hl = en&sa = X&ei =
    C_aQVJ3WAaGomgWuroLIBw&ved
    =0CDMQ6AEwAg#v = onepage&q = fourth%2Border%2Bepanechnikov&f = false

    second order univariate: k1(u)= 3/(4*sqrt(5)) * (1-u^2) , u^2 <=5
    fourth order univariate: k4, 1(u)= 3/(4*sqrt(5)) *
    (15/8 -7/8 u^2)(1-u^2), u^2 <=5
    sizth order univariate:  k6, 1(u) = 3/(4*sqrt(5)) *
    (175/64-105/32 u^2+ 231/320 u^4)* (1-u^2) , u^2 <=5
    """
    print filtersiz
    filtersize = (1, 1, 1)
    if not hasattr(filtersiz, "__len__"):
        filtersize = np.ones(3) * filtersiz
    else:
        if len(filtersiz) != 3:
            print "epanechnikov filter: size must be 1x3"
            raise ValueError
        filtersize = np.array(filtersiz[0:3])
    print bandwidth
    if not hasattr(bandwidth, "__len__"):
        sigma = np.ones(3) * bandwidth
    else:
        if len(bandwidth) != 3:
            print "epanechnikov filter: sigma must be 1x3"
            raise ValueError
        sigma = np.array(bandwidth[0:3])

    if np.mod(np.array(filtersize), 2).any():
        sz = (np.array(filtersize)) / 2
        xx = np.array(range(-int(sz[0]), int(sz[0]) + 1))
        yy = np.array(range(-int(sz[1]), int(sz[1]) + 1))
        zz = np.array(range(-int(sz[2]), int(sz[2]) + 1))
    else:
        sz = (np.array(filtersize) - 1) / 2.0
        xx = np.array(range(-int(sz[0]), int(sz[0])))
        yy = np.array(range(-int(sz[1]), int(sz[1])))
        zz = np.array(range(-int(sz[2]), int(sz[2])))
    mult_fact = np.ones((len(yy), len(xx), len(zz)))
    uu = xx[np.newaxis, :, np.newaxis] * mult_fact
    vv = yy[:, np.newaxis, np.newaxis] * mult_fact
    ww = zz[np.newaxis, np.newaxis, :] * mult_fact
    # if not hasattr(sigma, "__len__"):
    # if type(sigma) is float or type(sigma) is np.float64:
    #    epan = (0.75)*(1- (np.abs(uu**2 + vv**2 + ww**2))/(sigma**2))
    # else:
    zsquared = ((np.abs(uu) ** 2) / sigma[0] ** 2 +
                (np.abs(vv) ** 2) / sigma[1] ** 2 + (np.abs(ww) ** 2) / sigma[2] ** 2)
    epan = (0.75) * (1 - zsquared)
    if order == 4:
        print "Fourth order Epanechnikov"
        epan = (15.0 / 8.0 - (7.0 / 8.0) * zsquared) * epan
        # epan= (0.75)*(1.5 - 2.5((np.abs(uu)**2)/sigma[0]**2 +
        # (np.abs(vv)**2)/sigma[1]**2 + (np.abs(ww)**2)/sigma[2]**2))
    if order == 6:
        print "Sixth order Epanechnikov"
        epan = (175.0 / 64.0) * (1 - 6 * zsquared +
                                 (33.0 / 5.0) * zsquared * zsquared) * epan
    # else:
    #     print "epanechnikov filter can only have order={0, 4}"
    epan = epan * (epan > 0)
    epan = epan / ndimage.sum(epan)
    return epan.astype(np.float32)


def cplxepanechnikov_filter(real_input, imag_input, sigma_=1.87,
                            size_=3, mode_='reflect', order_=0):
    """cplxepanechnikov_filter Implementation of complex
    Epanechnikov filter using ndimage's generic_filter

         scipy.ndimage.filters.generic_filter(input, function,
     size=None, footprint=None, output=None, mode='reflect', cval=0.0,
     origin=0, extra_arguments=(), extra_keywords = None)[source]
     Calculates a multi-dimensional filter using the given function.

     At each element the provided function is called. The input values
     within the filter footprint at that element are passed to the
     function as a 1D array of double values.  ,(size_, size_, size_)
    """
    print "Complex Epanechnikov filter bandwidth ", sigma_, " size ", size_
    filtersize = (1, 1, 1)
    if not hasattr(size_, "__len__"):
        filtersize = np.ones(3) * size_
    else:
        if len(size_) != 3:
            print "cplxepanechnikov_filter: size must be 1x3"
            raise ValueError
        filtersize = np.array(size_)
#    if not hasattr(sigma_,"__len__"):
#        sigma_=np.ones(3)*sigma_
    epfilter = epanechnikov(filtersize, bandwidth=sigma_, order=order_)
#    def epfunc(x, y, z):
#         return max((0.75*(1 - x*x/sigma[0]**2
#         -y*y/sigma[1]**2 - z*z/sigma[2]**2), 0))
    real_img = np.empty_like(real_input)
    imag_img = np.empty_like(real_input)
    if real_input.ndim == 3:
        #        real_img = ndimage.filters.generic_filter(real_input, epfunc,
        #        size=(size_, size_, size_))
        #        imag_img = ndimage.filters.generic_filter(imag_input, epfunc,
        #        size=(size_, size_, size_))
        real_img = ndimage.convolve(real_input, epfilter, mode='reflect')
        imag_img = ndimage.convolve(imag_input, epfilter, mode='reflect')
    else:
        for echo in xrange(0, real_input.shape[4]):
            for acq in xrange(0, real_input.shape[3]):
                real_img[:, :, :, acq, echo] = ndimage.convolve(real_input[:, :, :, acq, echo],
                                                                epfilter,
                                                                mode=mode_)
                imag_img[:, :, :, acq, echo] = ndimage.convolve(imag_input[:, :, :, acq, echo],
                                                                epfilter,
                                                                mode=mode_)

    filtered_image = np.empty_like(real_input, dtype=np.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    return filtered_image
# end cplxepanechnikov_filter


def window_stdev(image, radius=2.5):
    """Standard deviation window filter
    https://stackoverflow.com/questions/18419871/improving-code-efficiency-
    standard-deviation-on-sliding-windows
    """
    from scipy.ndimage.filters import uniform_filter
    c1 = uniform_filter(image, radius * 2, mode='constant', origin=-radius)
    c2 = uniform_filter(
        image * image, radius * 2, mode='constant', origin=-radius)
    return ((c2 - c1 * c1) ** .5)[:-radius * 2 + 1, :-radius * 2 + 1]


def phase_std_filter(phase_image, wsize=5):
    """Implementation of local stdev filter for phase images
    """
    return np.fmin(window_stdev(phase_image, float(wsize)/2.0),
                   window_stdev(phase_image + np.pi), float(wsize)/2.0)


def cplxstdev_filter(real_input, imag_input, window_size=5):
    """cplxstdev_filter
    """
    real_img = np.empty_like(real_input)
    imag_img = np.empty_like(real_input)
    if window_size:
        radius = float(window_size/2)
    else:
        radius = 2.5

    if real_input.ndim == 3:
        real_img = window_stdev(real_input, radius)
        imag_img = window_stdev(imag_input, radius)
    else:
        for echo in xrange(0, real_input.shape[4]):
            for acq in xrange(0, real_input.shape[3]):
                real_img[:, :, :, acq, echo] = window_stdev(real_input[:, :, :, acq, echo], radius)
                imag_img[:, :, :, acq, echo] = window_stdev(imag_input[:, :, :, acq, echo], radius)

    filtered_image = np.empty_like(real_input, dtype=np.complex64)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    return filtered_image
# end cplxstdev_filter

# def unwrap():
# ##NOT WORKING YET
#     p2q2 = (p**2 + q**2)
#     phidash = ifftn(
#         fftn(np.cos(phiw) * ifft(p2q2*fftn(np.sin(phiw))))
#     ) - ifftn(fftn(np.sin(phiw) * ifft(p2q2*fftn(np.cos(phiw)))))
#     phij_update = phij + 2*np.pi*np.round((phidash-phij)/2*np.pi)


def swi(cmplx_input_image):
    """ Suseptibility weighted image
    """
    if np.iscomplexobj(cmplx_input_image):
        magn = np.abs(cmplx_input_image)
        phase = np.angle(cmplx_input_image)
        weight = phase / np.pi + 1.0
        weight = weight.clip(min=0.0, max=1.0)
        return magn * weight
    else:
        print 'Error swi2: input image not complex'
        return cmplx_input_image


def swi2(cmplx_input_image, order=2):
    """ Suseptibility weighted image - modified
    """
    if np.iscomplexobj(cmplx_input_image):
        magn = np.abs(cmplx_input_image)
        phase = np.angle(cmplx_input_image)
        weight = (phase / np.pi + 1.0) ** order
        weight = weight.clip(min=0.0, max=1.0)
        return magn * weight
    else:
        print 'Error swi2: input image not complex'
        return cmplx_input_image


def basicswi(cmplx_input_image, mask, order=2):
    """ Suseptibility-like weighted image -
    modified to use normalised local phase as the weighting
    """
    if np.iscomplexobj(cmplx_input_image):
        magn = np.abs(cmplx_input_image)
        phase = np.angle(cmplx_input_image)
        from scipy.ndimage.filters import uniform_filter
        normphase = uniform_filter(phase, 5.0, mode='constant', origin=-2.5)
        normphase = (normphase - ndimage.minimum(normphase)) / (ndimage.maximum(normphase) - ndimage.minimum(normphase))
        weight = (normphase + 1.0)
        weight = weight.clip(min=0.0, max=1.0)
        return magn * (weight**order) * mask
    else:
        print 'Error basicswi2: input image not complex'
        return np.abs(cmplx_input_image)


def normalise(data):
    """Normalise image
    """
    maxval = ndimage.maximum(data)
    minval = ndimage.minimum(data)
    print "Normalise max %f  min %f" % (maxval, minval)
    # return as float32
    return data.astype(np.float32)  # (data - minval) * (maxval - minval)


def save_nifti(image, basename):
    """Save image to NIFTI format
    """
    import nibabel as nib
    affine = np.eye(4)
    if image.ndim == 5:
        for echo in xrange(0, image.shape[4]):
            for channel in xrange(0, image.shape[3]):
                new_image = nib.Nifti1Image(np.abs(image[:, :, :, channel,
                                                         echo]), affine)
                new_image.set_data_dtype(np.float32)
                nib.save(new_image,
                         basename + '_' + str(channel) + str(echo) + '.nii.gz')
    else:
        new_image = nib.Nifti1Image(np.abs(image), affine)
        new_image.set_data_dtype(np.float32)
        nib.save(new_image, basename + '.nii.gz')


if __name__ == "__main__":

    import os
    # import sys
    # import math
    # import re
    import argparse
    import ReadProcpar as Procpar
    import ProcparToDicomMap
    # import RescaleFDF
    from ReadFID import *
    import nibabel as nib

    parser = argparse.ArgumentParser(
        usage=' ParseFDF.py -i "Input FDF directory"',
        description='''agilent2dicom is an FDF to Enhanced MR
        DICOM converter from MBI. ParseFDF takes header info from fdf files and
        adjusts the dicom dataset *ds* then rescales the image data.''')
    parser.add_argument('-i', '--inputdir', help='''Input directory name.
                        Must be an Agilent FDF image directory containing
                        procpar and *.fdf files''', required=True)
    parser.add_argument('-o', '--outputdir',
                        help='''Output directory name for DICOM files.''')
    parser.add_argument('-m', '--magnitude',
                        help='Magnitude component flag.', action="store_true")
    parser.add_argument('-p', '--phase',
                        help='Phase component flag.', action="store_true")
    parser.add_argument('-s', '--sequence',
                        help='''Sequence type (one of Multiecho, Diffusion,
                        ASL).''')
    parser.add_argument('-a', '--axis_order', help='Axis order eg 1, 0, 2.')
    parser.add_argument('-v', '--verbose',
                        help='Verbose comments.', action="store_true")
    args = parser.parse_args()
    # import ProcparToDicomMap as ptd

    procpar, procpartext = Procpar.ReadProcpar(os.path.join(args.inputdir,
                                                            'procpar'))
    ds, MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
    print "Rows: ", ds.Rows, " Columns: ", ds.Columns
    files = os.listdir(args.inputdir)
    fidfiles = [f for f in files if f.endswith('fid')]
    print "Number of FID files ", len(fidfiles)

    # for filename in fidfiles:
    print "Reading FID"
    filename = fidfiles[len(fidfiles) - 1]
    pp, hdr, dims, image_data_real, image_data_imag = readfid(args.inputdir,
                                                              procpar, args)
    print "Echoes: ", hdr['nEchoes'], " Channels: ", hdr['nChannels']
    affine = np.eye(4)
    # # affine[:3, :3]= np.arange(9).reshape((3, 3))
    # raw_data = nib.Nifti1Image(normalise(image_data_real), affine)
    # nib.save(raw_data, 'raw_data.nii.gz')

    print "Computing Original image (reconstruction)"
    if os.path.exists('raw_image_00.nii.gz'):
        nii = nib.load('raw_image_00.nii.gz')
        image[:, :, :, 0, 0] = nii.get_data()
        nii = nib.load('raw_image_01.nii.gz')
        image[:, :, :, 0, 1] = nii.get_data()
        nii = nib.load('raw_image_02.nii.gz')
        image[:, :, :, 0, 2] = nii.get_data()
    else:
        image, ksp = recon(pp, dims, hdr, image_data_real,
                           image_data_imag, args)

        if args.axis_order:
            image = RearrangeImage(image, args.axis_order, args)
            print "Transformed image shape: ", image.shape
            # np.delete(image)
            # image = imaget
        print "Saving raw image"
        save_nifti(np.abs(image), 'raw_image')

    print "Computing Gaussian (sigma=0.707) filtered image from Original image"
    gauss_filtered = cplxgaussian_filter(image.real, image.imag,
                                         np.sqrt(0.5), 0, 'reflect')
    print "Saving Gaussian image"
    save_nifti(normalise(np.abs(gauss_filtered)), 'gauss_image')

#    print "Computing SWI from Gaussian filtered image"
#    swi_image = swi(gauss_filtered)
#    print "Saving SWI image"
#    save_nifti(swi_image, 'swi_image')

    print "Computing Laplacian enhanced image from Original image"
    laplace_enhanced = gauss_filtered - cplxlaplacian_filter(gauss_filtered.real,
                                                             gauss_filtered.imag)
    print "Saving enhanced image g(x, y, z)=f(x, y, z) - Laplacian[f(x, y, z)]"
    save_nifti(normalise(np.abs(laplace_enhanced)), 'laplace_enhanced')

    print "Computing Gaussian filtered (sigma=1) image from Original image"
    gauss1_filtered = cplxgaussian_filter(image.real,
                                          image.imag, 1, 0, 'nearest')
    print "Saving Gaussian image"
    save_nifti(normalise(np.abs(gauss1_filtered)), 'gauss_image1')

    print "Computing Gaussian filtered image from Original image"
    gauss10_filtered = cplxgaussian_filter(image.real, image.imag,
                                           10.0, 0, 'nearest')
    print "Saving Gaussian image"
    save_nifti(normalise(np.abs(gauss10_filtered)), 'gauss_smooth')

    print "Computing Gaussian Laplace image from Smoothed image"
    Log_filtered = cplxgaussian_laplace(gauss_filtered.real,
                                        gauss_filtered.imag, 2.0)
    save_nifti(normalise(np.abs(Log_filtered)), 'Log_image')

    print "Computing LoG smooth  image from Gaussian (0.707) image"
    Log_smoothed = cplxgaussian_filter(Log_filtered.real,
                                       Log_filtered.imag, 10.0, 0,
                                       'reflect')
    print "Saving Gaussian image"
    save_nifti(normalise(np.abs(Log_smoothed)), 'Log_smooth')

    print "Computing Median filtered image"
    median_filtered = cplxmedian_filter(image.real, image.imag, 3.0)
    print "Saving Median"
    save_nifti(normalise(np.abs(median_filtered)), 'median_image')

    print "Computing Wiener filtered image"
    wiener_filtered = cplxwiener_filter(image.real, image.imag, 3, 0.0001)
    print "Saving Wiener image"
    save_nifti(normalise(np.abs(wiener_filtered)), 'wiener_image')

    print "Computing Epanechnikov filtered image"
    epanechnikov_filtered = cplxepanechnikov_filter(image.real,
                                                    image.imag,
                                                    np.sqrt(7.0 / 2.0),
                                                    (3, 3, 3))
    print "Saving Epanechnikov image"
    save_nifti(normalise(np.abs(epanechnikov_filtered)), 'epanechnikov_image')

    print "Computing Epanechnikov (4th order) filtered image"
    epanechnikov_filtered = cplxepanechnikov_filter(image.real,
                                                    image.imag,
                                                    np.sqrt(7.0 / 2.0),
                                                    (3, 3, 3), 1)
    print "Saving Epanechnikov image"
    save_nifti(normalise(np.abs(epanechnikov_filtered)), 'epa4th_image')

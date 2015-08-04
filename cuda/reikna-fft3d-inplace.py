"""
This example illustrates how to:
- attach a transformation to an FFT computation object that will make it
  operate on real-valued inputs.
"""


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

import numpy
from reikna.cluda import dtypes, any_api
from reikna.fft import FFT
from reikna.core import Annotation, Type, Transformation, Parameter


# Pick the first available GPGPU API and make a Thread on it.
api = any_api()
thr = api.Thread.create()
N = 512

arr = numpy.random.rand(N, N, N).astype(numpy.complex64)
arr.imag = numpy.random.rand(N, N, N)
arr[::2, :, :] = 0
arr[250:260, 250:260, 250:260] = 1
data_dev = thr.to_device(arr)
#res_dev = thr.empty_like(data_dev)
# Create the FFT computation and attach the transformation above to its input.i
# (A shortcut: using the array type saved in the transformation)
ifft = FFT(data_dev)
# fft.parameter.input.connect(trf, trf.output, new_input=trf.input)
cifft = ifft.compile(thr)


# Run the computation
#arr_dev = thr.to_device(arr)
#res_dev = thr.array(arr.shape, numpy.complex64)
# cifft(data_dev, data_dev,inverse=1)
# forward = data_dev.get()


# fft_arr = numpy.fft.fftn(arr)
# print numpy.linalg.norm(forward - fft_arr) / numpy.linalg.norm(fft_arr)


cifft(data_dev, data_dev, inverse=0)


result = data_dev.get() / 512**3
np.roll(np.roll(np.roll(result, 1, axis=2), 1, axis=1), 1, axis=0)
reference = numpy.fft.ifftn(arr)


print numpy.linalg.norm(result - reference) / numpy.linalg.norm(reference)
# < 1e-6


import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# imgplot = plt.imshow(np.abs(result[:,:, 250]), aspect='auto');plt.savefig('epank.jpg')
# clear the plot


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
ax1.imshow(np.abs(arr[:, :, 250]), aspect='auto')
ax1.set_title('Sharing x per column, y per row')
ax2.imshow(np.log10(np.abs(result[:, :, 250] / 512**3)), aspect='auto')
ax3.imshow(np.log10(np.abs(result[
           :, :, 250] / 512**3)) - np.log10(np.abs(reference[:, :, 250])), aspect='auto')
ax4.imshow(np.log10(np.abs(reference[:, :, 250])), aspect='auto')
plt.savefig('results.jpg')

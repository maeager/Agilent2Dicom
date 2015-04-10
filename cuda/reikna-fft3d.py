"""
This example illustrates how to:
- attach a transformation to an FFT computation object that will make it
  operate on real-valued inputs.
"""

import numpy
from reikna.cluda import dtypes, any_api
from reikna.fft import FFT
from reikna.core import Annotation, Type, Transformation, Parameter


# Pick the first available GPGPU API and make a Thread on it.
api = any_api()
thr = api.Thread.create()
N=256

arr = numpy.random.rand(N,N,N).astype(numpy.complex64)
arr.imag = numpy.random.rand(N,N,N)
# trf = get_complex_trf(arr)

data_dev = thr.to_device(arr)
res_dev = thr.empty_like(data_dev)
# Create the FFT computation and attach the transformation above to its input.
fft = FFT(data_dev) # (A shortcut: using the array type saved in the transformation)
# fft.parameter.input.connect(trf, trf.output, new_input=trf.input)
cfft = fft.compile(thr)


# Run the computation
#arr_dev = thr.to_device(arr)
#res_dev = thr.array(arr.shape, numpy.complex64)
cfft(res_dev, data_dev)
result = res_dev.get()

reference = numpy.fft.fftn(arr)

assert numpy.linalg.norm(result - reference) / numpy.linalg.norm(reference) < 1e-6




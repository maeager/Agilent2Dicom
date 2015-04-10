import numpy
from pyfft.cl import Plan
import pyopencl as cl
import pyopencl.array as cl_array

ctx = cl.create_some_context(interactive=False)
queue = cl.CommandQueue(ctx)
w = h = k = 512
plan = Plan((w,h,k), normalize=True, queue=queue)


# im1, im2 are the input 3d arrays of dtype complex64
im1 = numpy.random.rand(w,h,k).astype(numpy.complex64)
#im2 = numpy.random.rand(w,h,k).astype(numpy.complex64)
# forward transform on device
gpu_data = cl_array.to_device(queue, im1)

# forward transform
plan.execute(gpu_data.data) 
#result = gpu_data.get()

#Inverse transform:
plan.execute(gpu_data.data, inverse=True) 
result = gpu_data.get()

error = numpy.abs(numpy.sum(numpy.abs(im1) - numpy.abs(result)) / im1.size)

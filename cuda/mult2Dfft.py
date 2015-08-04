# http://wiki.tiker.net/PyCuda/Examples/2DFFT

import numpy
import scipy.misc
import numpy.fft as nfft
import multiprocessing

from pyfft.cuda import Plan
from pycuda.tools import make_default_context
import pycuda.tools as pytools
import pycuda.gpuarray as garray
import pycuda.driver as drv


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


class GPUMulti(multiprocessing.Process):

    def __init__(self, number, input_cpu, output_cpu):
        multiprocessing.Process.__init__(self)
        self.number = number
        self.input_cpu = input_cpu
        self.output_cpu = output_cpu

    def run(self):
        drv.init()
        a0 = numpy.zeros((p,), dtype=numpy.complex64)
        self.dev = drv.Device(self.number)
        self.ctx = self.dev.make_context()
# TO VERIFY WHETHER ALL THE MEMORY IS FREED BEFORE NEXT ALLOCATION (THIS
# DOES NOT HAPPEN IN MULTITHREADING)
        print drv.mem_get_info()
        self.gpu_a = garray.empty(
            (self.input_cpu.size,), dtype=numpy.complex64)
        self.gpu_b = garray.zeros_like(self.gpu_a)
        self.gpu_a = garray.to_gpu(self.input_cpu)
        plan = Plan(a0.shape, context=self.ctx)
        plan.execute(self.gpu_a, self.gpu_b, batch=p / m)
        self.temp = self.gpu_b.get()
        self.output_cpu.put(self.temp)
        self.output_cpu.close()
        self.ctx.pop()
        del self.gpu_a
        del self.gpu_b
        del self.ctx

        print "till the end %d" % self.number


p = 8192  # INPUT IMAGE SIZE (8192 * 8192)
# TO DIVIDE THE INPUT IMAGE INTO 4* (2048 * 8192) SIZED IMAGES (Depends on
# the total memory of your GPU)
m = 4
trans = 2  # FOR TRANSPOSE-SPLIT (TS) ALGORITHM WHICH loops 2 times


# INPUT IMAGE (GENERATE A 2d SINE WAVE PATTERN)
p_n = 8000  # No. OF PERIODS OF SINE WAVES
x = numpy.arange(0, p_n, float(p_n) / float(p))
a_i = 128 + 128 * numpy.sin(2 * numpy.pi * x)
a2 = numpy.zeros([p, p], dtype=numpy.complex64)
a2[::] = a_i
# TEST THE GENERATION OF INPUT IMAGE
scipy.misc.imsave("sine.bmp", numpy.absolute(a2))

# INITIALISE THE VARIABLES
# INPUT TO THE GPU (1d ARRAY)
a2_1 = numpy.zeros([m, p * p / m], dtype=numpy.complex64)
# VERY IMPORTANT
# STORE RESULT IN GPU (MULTIPROCESSING DOES NOT ALLOW SHARING AND HENCE
# THIS IS NEEDED FOR COMMUNICATION OF DATA)
output_cpu = multiprocessing.Queue()

b2pa = numpy.zeros([p / m, p, m], dtype=numpy.complex64)  # OUTPUT FROM GPU
# RESHAPED (8192*8192) OUTPUT
b2_a = numpy.zeros([p, p], dtype=numpy.complex64)

# NOW WE ARE READY TO KICK START THE GPU

# THE NO OF GPU'S PRESENT (CHANGE ACCORDING TO THE No.OF GPUS YOU HAVE)
# I KNOW THIS IS A BAD PRACTISE, BUT I COUNDN'T FIND ANY OTHER WAY(INIT
# CANNOT BE USED HERE)
num = 2

# THE TRANSPOSE-SPLIT ALGORITHM FOR FFT
tic()
for t in range(0, trans):
    for i in range(m):
        # DIVIDE AND RESHAPE THE INPUT IMAGE INTO 1D ARRAY
        a2_1[i, :] = a2[i * p / m:(i + 1) * p / m, :].flatten()

    for j in range(m / num):

        gpu_multi_list = []

# CREATE AND START THE MULTIPROCESS
        for i in range(num):
            # FEED THE DATA INTO THE GPU
            gpu_multi = GPUMulti(i, a2_1[i + j * num, :], output_cpu)
            gpu_multi_list.append(gpu_multi)
            gpu_multi.start()  # THERE YOU GO

# COLLECT THE OUTPUT FROM THE RUNNING MULTIPROCESS AND RESHAPE

        toc()
        for gpu_pro in gpu_multi_list:
            temp_b2_1 = output_cpu.get(gpu_pro)
            b2pa[:, :, gpu_pro.number + j *
                 num] = numpy.reshape(temp_b2_1, (p / m, p))
        gpu_multi.terminate()
        toc()
# RESHAPE AGAIN TO (8192 * 8192) IMAGE
    for i in range(m):
        b2_a[i * p / m:(i + 1) * p / m, :] = b2pa[:, :, i]
    toc()

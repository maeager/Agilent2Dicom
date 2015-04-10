
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



import numpy as np
## REINKA
import pycuda.driver as drv
import pycuda.tools
import pycuda.autoinit
from reikna.cluda import dtypes, any_api
from reikna.fft import FFT
from reikna.core import Annotation, Type, Transformation, Parameter
# create two timers so we can speed-test each approach
start = drv.Event()
end = drv.Event()
api = any_api()
thr = api.Thread.create()
N=512

start.record()
data_dev = thr.to_device(ksp)
ifft = FFT(data_dev)
cifft = ifft.compile(thr)
cifft(data_dev, data_dev,inverse=0)
result = np.fft.fftshift(data_dev.get() / N**3)
result = result[::-1,::-1,::-1]
result = np.roll(np.roll(np.roll(result,1,axis=2),1,axis=1),1,axis=0)
end.record()
end.synchronize()
secs = start.time_till(end)*1e-3
print "Reikna IFFT time and first three results:"
print "%s sec, %s" % (secs, str(np.abs(result[:3,0,0])))

start.record()
reference = np.fft.fftshift(np.fft.ifftn(ksp))
end.record()
end.synchronize()
secs = start.time_till(end)*1e-3
print "Numpy IFFTN time and first three results:"
print "%s sec, %s" % (secs, str(np.abs(reference[:3,0,0])))

print np.linalg.norm(result - reference) / np.linalg.norm(reference)


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

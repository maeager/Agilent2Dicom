#!/usr/bin/env python
"""
cplxfilter, recon, get_bit, readfid

Methods for complex filtering of 3D k-space data

- (C) Michael Eager 2014  (michael.eager@monash.edu)

"""


"""
scipy.ndimage.filters.gaussian_filter

scipy.ndimage.filters.gaussian_filter(input, sigma, order=0, output=None, mode='reflect', cval=0.0, truncate=4.0)[source]
Multidimensional Gaussian filter.

Parameters:
input : array_like
Input array to filter.
sigma : scalar or sequence of scalars
Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.
order : {0, 1, 2, 3} or sequence from same set, optional
The order of the filter along each axis is given as a sequence of integers, or as a single number. An order of 0 corresponds to convolution with a Gaussian kernel. An order of 1, 2, or 3 corresponds to convolution with the first, second or third derivatives of a Gaussian. Higher order derivatives are not implemented
output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
cval : scalar, optional
Value to fill past edges of input if mode is 'constant'. Default is 0.0
truncate : float
Truncate the filter at this many standard deviations. Default is 4.0.
Returns:
gaussian_filter : ndarray
Returned array of same shape as input.
Notes

The multidimensional filter is implemented as a sequence of one-dimensional convolution filters. The intermediate arrays are stored in the same data type as the output. Therefore, for output types with a limited precision, the results may be imprecise because intermediate results may be stored with insufficient precision.
"""

"""
scipy.ndimage.filters.gaussian_laplace

scipy.ndimage.filters.gaussian_laplace(input, sigma, output=None, mode='reflect', cval=0.0, **kwargs)[source]
Multidimensional Laplace filter using gaussian second derivatives.

Parameters:
input : array_like
Input array to filter.
sigma : scalar or sequence of scalars
The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.
output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
cval : scalar, optional
Value to fill past edges of input if mode is 'constant'. Default is 0.0
Extra keyword arguments will be passed to gaussian_filter().
"""


"""
scipy.ndimage.filters.laplace

scipy.ndimage.filters.laplace(input, output=None, mode='reflect', cval=0.0)[source]
N-dimensional Laplace filter based on approximate second derivatives.

Parameters:
input : array_like
Input array to filter.
output : array, optional
The output parameter passes an array in which to store the filter output.
mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
cval : scalar, optional
Value to fill past edges of input if mode is 'constant'. Default is 0.0

"""


"""
scipy.ndimage.fourier.fourier_gaussian

scipy.ndimage.fourier.fourier_gaussian(input, sigma, n=-1, axis=-1, output=None)[source]
Multi-dimensional Gaussian fourier filter.

The array is multiplied with the fourier transform of a Gaussian kernel.

Parameters:
input : array_like
The input array.
sigma : float or sequence
The sigma of the Gaussian kernel. If a float, sigma is the same for all axes. If a sequence, sigma has to contain one value for each axis.
n : int, optional
If n is negative (default), then the input is assumed to be the result of a complex fft. If n is larger than or equal to zero, the input is assumed to be the result of a real fft, and n gives the length of the array before transformation along the real transform direction.
axis : int, optional
The axis of the real transform.
output : ndarray, optional
If given, the result of filtering the input is placed in this array. None is returned in this case.
Returns:
fourier_gaussian : ndarray or None
The filtered input. If output is given as a parameter, None is returned.

"""


"""
Examples 

https://scipy-lectures.github.io/advanced/image_processing/
https://scipy-lectures.github.io/packages/scikit-image/

using PIL
http://nbviewer.ipython.org/github/mroberts3000/GpuComputing/blob/master/IPython/GaussianBlur.ipynb


"""

import os,sys,math
import re
import struct
import argparse
import ProcparToDicomMap



from scipy.fftpack import fftn,ifftn,fftshift,ifftshift
from scipy import ndimage
import numpy as np
import ReadProcpar

def cplxfilter(real_input, imag_input,sigma=0.707,order_=0,mode_='reflect', cval_=0.0):
    """CPLXFILTER gaussian filter of complex 3D image
     ndimage.filters.gaussian_filter is used to smooth real and imag components
    
    filtered_magnitude = cplxfilter(realimg,imagimg)

    :param sigma:  optimal sigma = 1/sqrt(2).  Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.
    :param order: {0, 1, 2, 3} or sequence from same set, optional.  The order of the filter along each axis is given as a sequence of integers, or as a single number. An order of 0 corresponds to convolution with a Gaussian kernel. An order of 1, 2, or 3 corresponds to convolution with the first, second or third derivatives of a Gaussian. Higher order derivatives are not implemented
    :param mode: {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional. The mode parameter determines how the array borders are handled, where cval is the value when mode is equal to 'constant'. Default is 'reflect'
    :param cval: [optional] Value to fill past edges of input if mode is 'constant'. Default is 0.0

    :return filtered_image:  complex 3D array of gaussian filtered image, 
    """
    # imgmag = load_untouch_nii(file1)
    # imph = load_untouch_nii(file2)
    # imcplx=(imgmag.img).*exp((1i).*(imph.img))
    # c = fftshift(ifftn(fftshift(imcplx)))
    # c = flipdim(flipdim(flipdim(c,1),2),3)
    #f = fspecial3('gaussian', [3 3 3], 1) %/sqrt(2))
    # compleximg = complex(imfilter(real(c), f), imfilter(imag(c), f))
    
    real_img = ndimage.filters.gaussian_filter(real_input, sigma, order=order_, mode=mode_, cval=cval_) #, truncate=4.0) truncate not supported in scipy 0.10
    imag_img = ndimage.filters.gaussian_filter(imag_input, sigma, order=order_, mode=mode_, cval=cval_) #, truncate=4.0)
    filtered_image = np.empty_like(real_input, dtype=complex)
    filtered_image.real = real_img
    filtered_image.imag = imag_img
    
    return filtered_image
# end cplxfilter    

def get_bit(value, bit_number):
    return (value & int(1 << (bit_number-1))) != 0
# end get_bit

def readfid(folder,pp=[]):
    """
      hdr,realpart,imagpart = readfid(folder[,pp])
     read kspace data from Agilent fid and procpar files.
     output img has dimensions [freq, phase, slice, channel, echo]
     note that the image may have to be circular shifted in the phase
     direction.
    """
    # error(nargchk(1,2,nargin))

    # warning off MATLAB:divideByZero

    #  get acqcycles and TE from procpar
    if not pp:
        pp = ReadProcpar.ReadProcpar(os.path.join(folder,'procpar'))

    #Define fid headers from procpar on first occasion
    hdr=dict()
    hdr['TE'] = pp['te']
    # FOV
    hdr['volumes'] = pp['lpe3']
    hdr['nEchoes'] = pp['ne']
    # rcvrs = re.search('y',pp['rcvrs'])
    hdr['nChannels'] = 1 #len(rcvrs);
    hdr['mode'] = '%dD' % pp['nD']
    if pp['nD'] == 2:
        hdr['FOVcm'] = [pp['lro'], pp['lpe']]
        hdr['dims'] = [pp['nf']/pp['ns'], pp['np']/2, pp['ns']]
        hdr['voxelmm'] = [hdr['lro']/hdr['dims'][0], hdr['lpe']/hdr['dims'][1], pp['thk']]*10
    elif pp['nD'] == 3:
        hdr['FOVcm'] = [pp['lro'], pp['lpe'], pp['lpe2']]
        hdr['dims'] = [pp['nf'], pp['np']/2, pp['nv2']]
        hdr['voxelmm'] = np.array(hdr['FOVcm']) / np.array(hdr['dims'])*10
        
    
    
    # open fid file
    f = open(os.path.join(folder,'fid'),"rb") 
    int16size = struct.calcsize('h')
    int32size = struct.calcsize('i')
    endian='>' # > for big-endian < for little
    # Read datafileheader using: x, = struct.unpack(type,binary) method
    # unpack returns a tuple and we only want the result
    hdr['nblocks'],   = struct.unpack(endian+'i',f.read(int32size)) #fid,1,'int32')
    hdr['ntraces'],   = struct.unpack(endian+'i',f.read(int32size)) #fid,1,'int32')
    hdr['np'],        = struct.unpack(endian+'i',f.read(int32size)) #fid,1,'int32')
    hdr['ebytes'],    = struct.unpack(endian+'i',f.read(int32size)) #fid,1,'int32')
    hdr['tbytes'],    = struct.unpack(endian+'i',f.read(int32size)) #fid,1,'int32')
    hdr['bbytes'],    = struct.unpack(endian+'i',f.read(int32size)) #fid,1,'int32')
    
    hdr['vers_id'],   = struct.unpack(endian+'h',f.read(int16size)) #fid,1,'int16')
    status=f.read(int16size)
    hdr['status'],    = struct.unpack(endian+'h',status) #fid,1,'int16')
    hdr['nbheaders'], = struct.unpack(endian+'i',f.read(int32size)) #fid,1,'int32')
    print 'status : ', hdr['status'], type(hdr['status']), type(status)
    hdr['s_data']    = int(get_bit(hdr['status'],1))
    hdr['s_spec']    = int(get_bit(hdr['status'],2))
    hdr['s_32']      = int(get_bit(hdr['status'],3))
    hdr['s_float']   = int(get_bit(hdr['status'],4))
    hdr['s_complex'] = int(get_bit(hdr['status'],5))
    hdr['s_hyper']   = int(get_bit(hdr['status'],6))
    print hdr # ['s_data'],hdr['s_spec'],hdr['s_32'],hdr['s_float'],hdr['s_complex'],hdr['s_hyper']
    
    dims=[0,0,0]
    # validate dimensions
    if pp['nD'] == 2:
        dims[0] = pp['np']/2        # num phase encode lines / 2
        dims[1] = pp['nf']/pp['ns'] # num frequency lines acquired / # echoes
        dims[2] = pp['ns']          # if 2D, num slices, else ni2
        if pp['ni2'] > 1:           # fse3d sequence has nD == 2, but is a 3d acquisition???
            dims[2] = pp['ni2']
    elif pp['nD'] == 3:
        dims[0] = pp['np']/2        # NUM phase encode lines / 2
        dims[1] = pp['nf']/pp['ne'] # num frequency lines acquired / # echoes
        dims[2] = pp['ni2'] 
    else:
        raise ValueError("Can only handle 2D or 3D files (based on procpar field nD)")   

    if hdr['np'] != pp['np'] or  hdr['ntraces'] != pp['nf'] or hdr['nblocks'] != pp['arraydim']:
        print 'NP ', hdr['np'], pp['np'], ' NF ', hdr['ntraces'], pp['nf'], ' Blocks ', hdr['nblocks'], pp['arraydim']
        raise ValueError("Cannot resolve fid header with procpar. We're probably not interpreting the procpar correctly.")
    # hdr['nChannels'] = hdr['nblocks']/hdr['acqcycles']
    hdr['nPhaseEncodes'] = hdr['ntraces']/hdr['nEchoes']
                                               
    # reset output structures
    RE = np.empty([hdr['np']/2, hdr['ntraces'], hdr['nblocks']], dtype=float)
    IM = np.empty([hdr['np']/2, hdr['ntraces'], hdr['nblocks']], dtype=float)
    
    # We have to read data every time in order to increment file pointer
    nchar = 0
    for i in xrange(0,hdr['nblocks']-1):
        # fprintf(1, repmat('\b',1,nchar))
        print  'reading block ', i+1,' of ', hdr['nblocks']

        # Read a block header
        header=dict()
        header['scale'],     = struct.unpack(endian+'h',f.read(int16size)) #fid,1,'int16')
        header['bstatus'],   = struct.unpack(endian+'h',f.read(int16size)) #fid,1,'int16')
        header['index'],     = struct.unpack(endian+'h',f.read(int16size)) #fid,1,'int16')
        header['mode'],      = struct.unpack(endian+'h',f.read(int16size)) #fid,1,'int16')
        header['ctcount'],   = struct.unpack(endian+'i',f.read(int32size)) #fid,1,'int32')
        header['lpval'],     = struct.unpack(endian+'f',f.read(int32size)) #fid,1,'float32')
        header['rpval'],     = struct.unpack(endian+'f',f.read(int32size)) #fid,1,'float32')
        header['lvl'],       = struct.unpack(endian+'f',f.read(int32size)) #fid,1,'float32')
        header['tlt'],       = struct.unpack(endian+'f',f.read(int32size)) #fid,1,'float32')
        print header
        if hdr['s_float'] == 1:
            # data = fread(fid,hdr.np*hdr.ntraces,'*float32')
            dtype_str=np.dtype(endian+'f4') # 'float32'
            print 'reading 32bit floats '
        elif hdr['s_32'] == 1:
            # data = fread(fid,hdr.np*hdr.ntraces,'*int32')
            dtype_str='int32'
            print 'reading 32bit int'
        else:
            # data = fread(fid,hdr.np*hdr.ntraces,'*int16')
            dtype_str ='int16' 
            str='reading 16bit int'
        data = np.fromfile(f,count=hdr['np']*hdr['ntraces'],dtype=dtype_str)

        data = np.reshape(data, [hdr['ntraces'],hdr['np']])
        RE[:,:,i] = data[:,:hdr['np']:2]   # hdr['np']
        IM[:,:,i] = data[:,1:hdr['np']:2]  # hdr['np']
        #break
    f.close()
    print i
    #hdr.pp = pp
    print "Data Row:   %.15g %.15g %.15g %.15g %.15g %.15g %.15g" % (data[0,0],data[0,1],  data[0,2],data[0,3],  data[0,4],data[0,1022],data[511,1022])
    print "Data Col:   %.15g %.15g %.15g %.15g %.15g %.15g %.15g" % (data[0,0],data[1,0],  data[2,0],data[3,0],  data[4,0],data[0,1023],data[511,1023])
    print "RE : %.15g  %.15g %.15g %.15g %.15g %.15g %.15g" % (RE[0,0,0],RE[0,1,0],  RE[0,2,0],RE[1,0,0],  RE[2,0,0],RE[0,511,0],RE[511,511,0])
    print "IM : %.15g  %.15g %.15g %.15g %.15g %.15g %.15g" % (IM[0,0,0],IM[0,1,0],  IM[0,2,0],IM[1,0,0],  IM[2,0,0],IM[0,511,0],IM[511,511,0])
    return pp,hdr,dims,RE,IM
# end readfid


def recon(pp,dims,hdr,RE,IM):
    print 'Performing fourier transform...'

    ksp = np.empty([dims[0], dims[1], dims[2], hdr['nChannels'], hdr['nEchoes']], dtype=complex) #float32
    img = np.empty([dims[0], dims[1], dims[2], hdr['nChannels'], hdr['nEchoes']], dtype=complex) #float32

    if pp['nD'] == 2 and pp['ni2'] == 1:
        for echo in xrange(0,int(hdr['nEchoes']-1)):
            for channel in xrange(0,int(hdr['nChannels']-1)):
                for islice in xrange(0,dims(3)-1):
                    # ksp(:,:,islice,channel,echo) = complex(RE(:,echo:hdr['nEchoes']:end,channel:hdr['nChannels']:end), IM(:,echo:hdr['nEchoes']:end,channel:hdr['nChannels']:end))
                    ksp[:,:,islice,channel,echo].real = RE[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    ksp[:,:,islice,channel,echo].imag = IM[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    
                    img[:,:,islice,channel,echo] = fftshift(ifftn(ifftshift(ksp[:,pp['pelist']-minimum(pp['pelist']),islice,channel,echo])))
    else: #if pp.nD == 3
        if hdr['nEchoes'] == 1 and hdr['nChannels']== 1:
            ksp = np.empty([dims[0], dims[1], dims[2]], dtype=complex) #float32
            img = np.empty([dims[0], dims[1], dims[2]], dtype=complex) #float32
            ksp[:,:,:].real = RE  #[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
            ksp[:,:,:].imag = IM  #[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
            if 'pelist' in pp.keys():
                img[:,:,:] = fftshift(ifftn(ifftshift(ksp[:,pp['pelist']-min(pp['pelist']),:])))
            else:
                img[:,:,:] = fftshift(ifftn(ifftshift(ksp[:,:,:])))
        else:
            for echo in xrange(0,int(hdr['nEchoes']-1)):
                for n in xrange(0,int(hdr['nChannels']-1)):
                    ksp[:,:,:,n,echo].real = RE[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    ksp[:,:,:,n,echo].imag = IM[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    if 'pelist' in pp.keys():
                        img[:,:,:,n,echo] = fftshift(ifftn(ifftshift(ksp[:,pp['pelist']-min(pp['pelist']),:,n,echo])))
                    else:
                        img[:,:,:,n,echo] = fftshift(ifftn(ifftshift(ksp[:,:,:,n,echo])))
        
    return img,ksp
#end recon
    
def normalise(data):
    max = ndimage.maximum(data)
    min = ndimage.minimum(data)
    print "Normalise max %f  min %f" % (max, min)
    return data #(data - min) * (max - min) 

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage=' ParseFDF.py -i "Input FDF directory"',description='agilent2dicom is an FDF to Enhanced MR DICOM converter from MBI. ParseFDF takes header info from fdf files and adjusts the dicom dataset *ds* then rescales the image data.')
    parser.add_argument('-i','--inputdir', help='Input directory name. Must be an Agilent FDF image directory containing procpar and *.fdf files',required=True);
    parser.add_argument('-o','--outputdir', help='Output directory name for DICOM files.')
    parser.add_argument('-m','--magnitude', help='Magnitude component flag.',action="store_true");
    parser.add_argument('-p','--phase', help='Phase component flag.',action="store_true");
    parser.add_argument('-s','--sequence', help='Sequence type (one of Multiecho, Diffusion, ASL.');
    
    parser.add_argument('-v','--verbose', help='Verbose.',action="store_true");
    args = parser.parse_args()
    
    import ReadProcpar
    import RescaleFDF
    import nibabel as nib
    # import ReadFDF as rf
    # import ProcparToDicomMap as ptd


    procpar, procpartext = ReadProcpar.ReadProcpar(os.path.join(args.inputdir,'procpar'))
    ds,MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
    print "Rows: ", ds.Rows, " Columns: ", ds.Columns
    

    files = os.listdir(args.inputdir)
    fidfiles = [ f for f in files if f.endswith('fid') ]
    print "Number of FID files ", len(fidfiles)

    # for filename in fidfiles:
    filename = fidfiles[len(fidfiles)-1]
    pp,hdr,dims,image_data_real,image_data_imag=readfid(args.inputdir,procpar)

    affine = np.eye(4)
    # # affine[:3,:3]= np.arange(9).reshape((3,3))
    raw_data=nib.Nifti1Image(normalise(image_data_real),affine)
    nib.save(raw_data,'raw_data.nii.gz')

    
    image,ksp=recon(pp,dims,hdr,image_data_real,image_data_imag)

    raw_image=nib.Nifti1Image(normalise(np.abs(image)),affine)
    nib.save(raw_image,'raw_image.nii.gz')

    raw_ksp=nib.Nifti1Image(normalise(np.abs(ksp)),affine)
    nib.save(raw_ksp,'raw_ksp.nii.gz')

    
    image_filtered = cplxfilter(image.real,image.imag,1.0)


    new_image = nib.Nifti1Image(normalise(np.abs(image_filtered)),affine)
    nib.save(new_image,'new_image.nii.gz')
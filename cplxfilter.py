import os,sys,math
import argparse
import ProcparToDicomMap
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

from scipy.fftpack import fftn,ifftn,fftshift,ifftshift
from scipy import ndimage
import numpy as np
import ReadProcpar

def cplxfilter(real_input, imag_input):
    """
    # file1,file2,file3):
    #
    # cplxfilter('nifti/mag1.nii.gz','nifti/pha1.nii.gz','nifti/abs1.nii.gz')
    #
    """
    # imgmag = load_untouch_nii(file1)
    # imph = load_untouch_nii(file2)
    # imcpl=(imgmag.img).*exp((1i).*(imph.img))
    # c = fftshift(ifftn(fftshift(imcpl)))
    # c = flipdim(flipdim(flipdim(c,1),2),3)
    #f = fspecial3('gaussian', [3 3 3], 1) %/sqrt(2))
    # compleximg = complex(imfilter(real(c), f), imfilter(imag(c), f))
    
    real_img = ndimage.filters.gaussian_filter(real_input, sigma, order=0, mode='reflect', cval=0.0, truncate=4.0)
    imag_img = ndimage.filters.gaussian_filter(imag_input, sigma, order=0, mode='reflect', cval=0.0, truncate=4.0)
    result = np.empty(real_input.shape[:-1], dtype=complex)
    result.real = real_img; result.imag = imag_img
    #compleximg = complex(imfilter(real(imcpl), f), imfilter(imag(imcpl), f))
    #nii=imgmag
    #nii.img=abs(compleximg)
    #save_untouch_nii(nii,file3)
    return result
    

def get_bitos(value, bit_number):
    return (value & (1 << bit_number)) != 0
   #     )))(byteval,idx):
   # return ((byteval&(1<<idx))!=0);))))

def readfid(folder,pp=[]):
    """
    function [img, hdr, ksp] = readfid(folder)
    % [img, hdr, ksp] = readfid(folder)
    % read kspace data from Agilent fid and procpar files.
    % output img has dimensions [freq, phase, slice, channel, echo]
    % note that the image may have to be circular shifted in the phase
    % direction.
    """
    # error(nargchk(1,2,nargin))

    # warning off MATLAB:divideByZero

    #  get acqcycles and TE from procpar
    if empty(pp):
        pp = ReadProcpar.ReadProcpar(folder+'/procpar')
    '''
    % fid = fopen([folder '/procpar'],'r')
    % %hdr.acqcycles = []
    % %hdr.nEchoes = 1
    % while ~feof(fid)
    %     line = fgetl(fid)
    %     [par,line] = strtok(line)
    %     if strcmp(par,'acqcycles')
    %         line = fgetl(fid)
    %         [~,hdr.acqcycles] = strtok(line)
    %         hdr.acqcycles = str2num(hdr.acqcycles)
    %     elseif strcmp(par,'TE')
    %         line = fgetl(fid)
    %         hdr.nEchoes = str2num(strtok(line))
    %     end
    %     if ~isempty(hdr.acqcycles) && ~isempty(hdr.nEchoes)
    %         break
    %     end
    % end

    % fclose(fid)
    '''
    # open fid file
    f = open(folder+'/fid','r') #'b'

    # Read datafileheader
    hdr = dict()
    hdr['nblocks']   = f.read(4) #fid,1,'int32')
    hdr['ntraces']   = f.read(4) #fid,1,'int32')
    hdr['np']        = f.read(4) #fid,1,'int32')
    hdr['ebytes']    = f.read(4) #fid,1,'int32')
    hdr['tbytes']    = f.read(4) #fid,1,'int32')
    hdr['bbytes']    = f.read(4) #fid,1,'int32')
    hdr['vers_id']   = f.read(2) #fid,1,'int16')
    hdr['status']    = f.read(2) #fid,1,'int16')
    hdr['nbheaders'] = f.read(4) #fid,1,'int32')

    hdr['s_data']    = get_bit(hdr['status'],1)
    hdr['s_spec']    = get_bit(hdr['status'],2)
    hdr['s_32']      = get_bit(hdr['status'],3)
    hdr['s_float']   = get_bit(hdr['status'],4)
    hdr['s_complex'] = get_bit(hdr['status'],5)
    hdr['s_hyper']   = get_bit(hdr['status'],6)
    dims=[0,0,0]
    # validate dimensions
    if pp['nD'] == 2:
        dims[1] = pp['np']/2 # num phase encode lines / 2
        dims[2] = pp['nf']/pp['ns'] # num frequency lines acquired / # echoes
        dims[3] = pp['ns'] # if 2D, num slices, else ni2
        if pp['ni2'] > 1: # fse3d sequence has nD == 2, but is a 3d acquisition???
            dims[3] = pp['ni2']
    elif pp['nD'] == 3:
        dims[1] = pp['np']/2  # NUM phase encode lines / 2
        dims[2] = pp['nf']/pp['ne'] # num frequency lines acquired / # echoes
        dims[3] = pp['ni2'] 
    else:
        error('Can only handle 2D or 3D files (based on procpar field nD)')
    

    if hdr['np'] != pp['np'] or  hdr['ntraces'] != pp['nf'] or hdr['nblocks'] != pp['arraydim']:
        error('Cannot resolve fid header with procpar. We''re probably not interpreting the procpar correctly.')


    # reset output structures
    RE = np.empty([hdr['np']/2, hdr['ntraces'], hdr['nblocks']], dtype=float32)
    IM = np.empty([hdr['np']/2, hdr['ntraces'], hdr['nblocks']], dtype=float32)
    # RE = zeros(hdr['np']/2, hdr['ntraces'], hdr['nblocks'], 'single')
    # IM = zeros(hdr['np']/2, hdr['ntraces'], hdr['nblocks'], 'single')

    # We have to read data every time in order to increment file pointer
    nchar = 0
    for i in xrange(0,hdr['nblocks']-1):
        # fprintf(1, repmat('\b',1,nchar))
        print  'reading block %d of %d' % i+1, hdr['nblocks']

        # Read a block header
        header=dict()
        header['scale']     = f.read(2) #fid,1,'int16')
        header['bstatus']   = f.read(2) #fid,1,'int16')
        header['index']     = f.read(2) #fid,1,'int16')
        header['mode']      = f.read(2) #fid,1,'int16')
        header['ctcount']   = f.read(4) #fid,1,'int32')
        header['lpval']     = f.read(4) #fid,1,'float32')
        header['rpval']     = f.read(4) #fid,1,'float32')
        header['lvl']       = f.read(4) #fid,1,'float32')
        header['tlt']       = f.read(4) #fid,1,'float32')
        print header
        if hdr['s_float'] == 1:
            # data = fread(fid,hdr.np*hdr.ntraces,'*float32')
            dtype_str='float32'
            print 'reading floats'
        elif hdr['s_32'] == 1:
            # data = fread(fid,hdr.np*hdr.ntraces,'*int32')
            dtype_str='int32'
            print 'reading 32bit'
        else:
            # data = fread(fid,hdr.np*hdr.ntraces,'*int16')
            dtype_str ='int16' 
            str='reading 16bit'
        data = numpy.fromfile(hdr['np']*hdr['ntraces'],dtype=dtype_str)

    data = numpy.reshape(data,[hdr['np'], hdr['ntraces']])
    RE[:,:,i] = data[:hdr['np']:2,:]
    IM[:,:,i] = data[1:hdr['np']:2,:]
    f.close(fid)
    return hdr,RE,IM



def recon(pp,dims,hdr,RE,IM):
    print '\nPerforming fourier transform...'

    #hdr['nChannels'] = hdr['nblocks']/hdr['acqcycles']
    hdr['nPhaseEncodes'] = hdr['ntraces']/hdr['nEchoes']

    #img = zeros(size(RE,1), size(RE,2)/hdr['nEchoes'], size(RE,3)/hdr['nChannels'], hdr['nChannels'], hdr['nEchoes'], 'single')
    ksp = np.empty([dims(1), dims(2), dims(3), hdr['nChannels'], hdr['nEchoes']], dtype=float32)
    ksp = np.empty([dims(1), dims(2), dims(3), hdr['nChannels'], hdr['nEchoes']], dtype=float32)
    #ksp = zeros(dims(1), dims(2), dims(3), hdr['nChannels'], hdr['nEchoes'], 'single')
    #img = zeros(dims(1), dims(2), dims(3), hdr['nChannels'], hdr['nEchoes'], 'single')

    # img = complex(RE, IM)
    # return

    if pp['nD'] == 2 and pp['ni2'] == 1:
        for echo in xrange(1,hdr['nEchoes']):
            for channel in xrange(1,hdr['nChannels']):
                for islice in xrange(1,dims(3)):
                    # ksp(:,:,islice,channel,echo) = complex(RE(:,echo:hdr['nEchoes']:end,channel:hdr['nChannels']:end), IM(:,echo:hdr['nEchoes']:end,channel:hdr['nChannels']:end))
                    ksp[:,:,islice,channel,echo].real = RE[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    ksp[:,:,islice,channel,echo].imag = IM[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    
                    img[:,:,islice,channel,echo] = fftshift(ifftn(ifftshift(ksp[:,pp['pelist']-minimum(pp['pelist'])+1,islice,channel,echo])))
    else: #if pp.nD == 3
        for echo in xrange(0,hdr['nEchoes']-1):
            for n in xrange(0,hdr['nChannels']-1):
                ksp[:,:,:,n,echo].real = RE[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                ksp[:,:,:,n,echo].imag = IM[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                if 'pelist' in pp.keys():
                    img[:,:,:,n,echo] = fftshift(ifftn(ifftshift(ksp[:,pp['pelist']-min(pp['pelist'])+1,:,n,echo])))
                else:
                    img[:,:,:,n,echo] = fftshift(ifftn(ifftshift(ksp[:,:,:,n,echo])))

            #hdr.pp = pp

    return img


if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage=' ParseFDF.py -i "Input FDF directory"',description='agilent2dicom is an FDF to Enhanced MR DICOM converter from MBI. ParseFDF takes header info from fdf files and adjusts the dicom dataset *ds* then rescales the image data.')
    parser.add_argument('-i','--inputdir', help='Input directory name. Must be an Agilent FDF image directory containing procpar and *.fdf files',required=True);
    parser.add_argument('-o','--outputdir', help='Output directory name for DICOM files.')
    
    parser.add_argument('-v','--verbose', help='Verbose.',action="store_true");
    args = parser.parse_args()
    
    import ReadProcpar
    import RescaleFDF
    # import ReadFDF as rf
    # import ProcparToDicomMap as ptd


    procpar, procpartext = ReadProcpar.ReadProcpar(os.path.join(args.inputdir,'procpar'))
    ds,MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
    print "Rows: ", ds.Rows, " Columns: ", ds.Columns
    

    files = os.listdir(args.inputdir)
    fidfiles = [ f for f in files if f.endswith('.fid') ]
    print "Number of FID files ", len(fidfiles)

    # for filename in fidfiles:
    filename = fidfiles[len(fidfiles)-1]
    fid_properties,image_real,image_imag=readfid(os.path.join(args.inputdir,filename))

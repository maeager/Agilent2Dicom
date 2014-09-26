#!/usr/bin/env python

"""ReadFID is used to read Agilent FID image files

   (c) 2014  Michael Eager (michael.eager@monash.edu)
"""

import os,sys
import struct
import numpy
import argparse
import ReadProcpar

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
        hdr['voxelmm'] = numpy.array(hdr['FOVcm']) / numpy.array(hdr['dims'])*10
        
    
    
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
    RE = numpy.empty([hdr['np']/2, hdr['ntraces'], hdr['nblocks']], dtype=float)
    IM = numpy.empty([hdr['np']/2, hdr['ntraces'], hdr['nblocks']], dtype=float)
    
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
            dtype_str=numpy.dtype(endian+'f4') # 'float32'
            print 'reading 32bit floats '
        elif hdr['s_32'] == 1:
            # data = fread(fid,hdr.np*hdr.ntraces,'*int32')
            dtype_str='int32'
            print 'reading 32bit int'
        else:
            # data = fread(fid,hdr.np*hdr.ntraces,'*int16')
            dtype_str ='int16' 
            str='reading 16bit int'
        data = numpy.fromfile(f,count=hdr['np']*hdr['ntraces'],dtype=dtype_str)

        data = numpy.reshape(data, [hdr['ntraces'],hdr['np']])
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
    """recon
    Reconstruct k-space data into N-D image
    :param pp:   procpar dictionary
    :param dims: dimension array
    :param hdr: header info in fid
    :param RE: real component of image k-space
    :param IM: imaginary component of image k-space
    """
    print 'Reconstructing image'

    ksp = numpy.empty([dims[0], dims[1], dims[2], hdr['nChannels'], hdr['nEchoes']], dtype=complex) #float32
    img = numpy.empty([dims[0], dims[1], dims[2], hdr['nChannels'], hdr['nEchoes']], dtype=complex) #float32

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
            ksp = numpy.empty([dims[0], dims[1], dims[2]], dtype=complex) #float32
            img = numpy.empty([dims[0], dims[1], dims[2]], dtype=complex) #float32
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

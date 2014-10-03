#!/usr/bin/env python

"""ReadFID is used to read Agilent FID image files

   -  Michael Eager (michael.eager@monash.edu)
"""
"""
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

import os,sys
import struct
import numpy
import argparse
import ReadProcpar
from scipy.fftpack import fftn,ifftn,fftshift,ifftshift
#from scipy import ndimage


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
    print 'Dimensions: ',dims, hdr['dims']

    if hdr['np'] != pp['np'] or  hdr['ntraces'] != pp['nf'] or hdr['nblocks'] != pp['arraydim']:
        print 'NP ', hdr['np'], pp['np'], ' NF ', hdr['ntraces'], pp['nf'], ' Blocks ', hdr['nblocks'], pp['arraydim']
        raise ValueError("Cannot resolve fid header with procpar. We're probably not interpreting the procpar correctly.")
    # hdr['nChannels'] = hdr['nblocks']/hdr['acqcycles']
    hdr['nPhaseEncodes'] = hdr['ntraces']/hdr['nEchoes']
    hdr['rank'] = pp['nD']                                           

    # Axial    Orientation of the target scans is transverse (theta=0, psi=0, phi=0).
    # Coronal  Orientation of the target scans is coronal (theta=90, psi=0, phi=0).
    # Sagittal Orientation of the target scans is sagittal (theta=90, psi=90, phi=0).
#                   Table 15. Image Planning Parameters
# Parameter    Description
#              Gap between slices.
# gap
#              FOV of first phase encoding dimension, used for lope, vox1.
# lpe
#              FOV of second phase encoding dimension, used for lpe2, vox3.
# lpe2
#              FOV of read-out dimension, used for lro and vox2.
# lro
#              Number of slices.
# ns
#              Euler angle, used for phi, sphi, vphi.
# phi
#              Euler angle, used for psi, spsi, vpsi.
# psi
#              Shift of stack center along first phase encoding dimension.
# ppe
#              Used for ppe, pos1.
#              Shift of stack center along readout dimension. Used for pro, pos2.
# pro
#              Shifts of slices along z axis.
# pss
#              Shift of stack center along z, the axis perpendicular to the plane.
# pss0
#              Used for pss0 and pos3.
#              Fan angle of radial slices.
# radialAngles
#              Euler angle, used for theta, stheta, vtheta.
# theta
#              Thickness of slices or saturation bands (satthk).
# thk

    # Create FDF-like header variables
    hdr['location'] = [-pp['pro'], pp['ppe'], pp['pss'] ]
    hdr['span'] = [pp['lro'], pp['lpe'], pp['lpe2']]
    hdr['roi'] = [pp['lro'], pp['lpe'], pp['lpe2']]
    hdr['origin']=[-pp['pro']-pp['lro']/2,  pp['ppe']-pp['lpe']/2,  pp['pss']-pp['lpe2']/2]
    if pp['orient']=="sag":
        hdr['orientation']= [0,0,1,1,0,0,0,1,0]
    else:
        hdr['orientation']= [1,0,0,0,1,0,0,0,1]

    # reset output structures
    RE = numpy.empty([ hdr['np']/2,hdr['ntraces'], hdr['nblocks']], dtype=float)
    IM = numpy.empty([ hdr['np']/2,hdr['ntraces'], hdr['nblocks']], dtype=float)
    
    # We have to read data every time in order to increment file pointer
    nchar = 0
    for iblock in xrange(0,hdr['nblocks']):
        # fprintf(1, repmat('\b',1,nchar))
        print  'reading block ', iblock+1,' of ', hdr['nblocks']
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
        data = numpy.fromfile(f,count=hdr['np']*hdr['ntraces'],dtype=dtype_str)
        print data.ndim, data.shape
        data = numpy.reshape(data, [hdr['ntraces'],hdr['np']])
        RE[:,:,iblock] = numpy.matrix(data[:,:hdr['np']:2]).T   # hdr['np'] #[::2,:] #
        IM[:,:,iblock] = numpy.matrix(data[:,1:hdr['np']:2]).T  # hdr['np'] #[1::2,:]      #
        #break
    f.close()
    print iblock
    #hdr.pp = pp
    print "Data Row 1:   %.5g %.5g %.5g %.5g %.5g  %.5g ... %.5g %.5g %.5g %.5g" % (data[0,0],data[0,1],  data[0,2],data[0,3],  data[0,4], data[0,5], data[0,-4], data[0,-3], data[0,-2], data[0,-1])
    print "Data Col 1:   %.5g %.5g %.5g %.5g %.5g %.5g  ... %.5g %.5g %.5g %.5g" % (data[0,0],data[1,0],  data[2,0],data[3,0],  data[4,0], data[5,0], data[-4,0],data[-3,0], data[-2,0], data[-1,0])
    print "Data Row -1:   %.5g %.5g %.5g %.5g %.5g %.5g ... %.5g %.5g %.5g %.5g" % (data[-1,0],data[-1,1],  data[-1,2],data[-1,3],  data[-1,4],data[-1,5], data[-1,-4], data[-1,-3], data[-1,-2],data[-1,-1])
    print "Data Col -1:   %.5g %.5g %.5g %.5g %.5g %.5g ... %.5g %.5g %.5g %.5g" % (data[0,-1],data[1,-1],  data[2,-1],data[3,-1],  data[4,-1], data[5,-1],data[-4,-1], data[-3,-1], data[-2,-1],data[-1,-1])
    print "RE : %.5g  %.5g %.5g | %.5g %.5g | %.5g %.5g |%.5g" % (RE[0,0,iblock],RE[0,1,iblock],  RE[0,2,iblock],RE[1,0,iblock],  RE[2,0,iblock], RE[0,-1,iblock], RE[-1,0,iblock], RE[-1,-1,iblock])
    print "IM : %.5g  %.5g %.5g | %.5g %.5g | %.5g %.5g |%.5g" % (IM[0,0,iblock],IM[0,1,iblock],  IM[0,2,iblock],IM[1,0,iblock],  IM[2,0,iblock], IM[0,-1,iblock], IM[-1,0,iblock], IM[-1,-1,iblock])
    return pp,hdr,dims,RE,IM
# end readfid


def recon(pp,dims,hdr,RE,IM):
    """recon Reconstruct k-space image data into N-D image
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
        for echo in xrange(0,int(hdr['nEchoes'])):
            for channel in xrange(0,int(hdr['nChannels'])):
                for islice in xrange(0,dims(3)):
                    # ksp(:,:,islice,channel,echo) = complex(RE(:,echo:hdr['nEchoes']:end,channel:hdr['nChannels']:end), IM(:,echo:hdr['nEchoes']:end,channel:hdr['nChannels']:end))
                    ksp[:,:,islice,channel,echo].real = RE[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    ksp[:,:,islice,channel,echo].imag = IM[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    
                    img[:,:,islice,channel,echo] = fftshift(ifftn(ifftshift(ksp[:,pp['pelist']-minimum(pp['pelist']),islice,channel,echo])))
    else: #if pp.nD == 3
        if hdr['nEchoes'] == 1 and hdr['nChannels'] == 1:
            ksp = numpy.empty([dims[0], dims[1], dims[2]], dtype=complex) #float32
            img = numpy.empty([dims[0], dims[1], dims[2]], dtype=complex) #float32
            ksp[:,:,:].real = RE  #[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
            ksp[:,:,:].imag = IM  #[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
            if 'pelist' in pp.keys():
                img[:,:,:] = fftshift(ifftn(ifftshift(ksp[:,pp['pelist']-min(pp['pelist']),:])))
            else:
                img[:,:,:] = fftshift(ifftn(ifftshift(ksp[:,:,:])))
        else:
            for echo in xrange(0,int(hdr['nEchoes'])):
                for n in xrange(0,int(hdr['nChannels'])):
                    print "Processing echo ",echo," channel ",n
                    ksp[:,:,:,n,echo].real = RE[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    ksp[:,:,:,n,echo].imag = IM[:,echo::hdr['nEchoes'],n::hdr['nChannels']]
                    if 'pelist' in pp.keys():
                        img[:,:,:,n,echo] = fftshift(ifftn(ifftshift(ksp[:,pp['pelist']-min(pp['pelist']),:,n,echo])))
                    else:
                        img[:,:,:,n,echo] = fftshift(ifftn(ifftshift(ksp[:,:,:,n,echo])))
        
    return img,ksp
#end recon


def RescaleImage(ds,image_data,args):

        # Read in data from all files to determine scaling
    datamin = float("inf")
    datamax = float("-inf")

    # RescaleSlope of phase imgs set to [-pi,pi]
    if ds.ImageType[2] == "PHASE MAP" or \
            (hasattr(ds,'ComplexImageComponent') and ds.ComplexImageComponent == 'PHASE'):
        # this implies either args.phase is on or procpar['imPH']=='y'
        datamin = -math.pi
        datamax = math.pi
    else:
        datamin = numpy.min([datamin,image_data.min()])
        datamax = numpy.max([datamax,image_data.max()])


    RescaleIntercept = datamin
    # Numpy recast to int16 with range (-32768 or 32767)
    if ds.PixelRepresentation == 0:
        slope_factor = UInt16MaxRange
    else:
        slope_factor = Int16MaxRange
    RescaleSlope = (datamax - datamin) / slope_factor 

    
    if args.verbose:
        print "Rescale data to uint16"
        print "Intercept: ", RescaleIntercept, "  Slope: ", RescaleSlope
        print "Current data min: ", image_data.min(), " max ", image_data.max()
    image_data = (image_data - RescaleIntercept) / RescaleSlope
    image_data = image_data.astype(numpy.int16)

    ## Adjusting Dicom parameters for rescaling
        # Rescale intercept string must not be longer than 16
    if len(str(RescaleIntercept))>16:
        print "Cropping rescale intercept from ", str(RescaleIntercept), " to ", str(RescaleIntercept)[:15]
    ds.RescaleIntercept = str(RescaleIntercept)[:15]                                #(0028,1052) Rescale Intercept
    # Rescale slope string must not be longer than 16
    if len(str(RescaleSlope))>16:
        print "Cropping rescale slope from ", str(RescaleSlope), " to ", str(RescaleSlope)[:15]
    ds.RescaleSlope = str(RescaleSlope)[:15]    #(0028,1053) Rescale Slope

    return ds,image_data
#end RescaleImage


def ParseFID(ds,fid_properties,procpar,args):
    """
    ParseFID modify the dicom dataset structure based on FID
    header information and procpar information.
    
    Comment text copied from VNMRJ Programming.pdf

    :param ds:       Dicom dataset
    :param fid_properties: Dict of fid header label/value pairs
    :param procpar:  Dict of procpar label/value pairs
    :param args:     Argparse object
    :return ds:      Return updated dicom dataset struct
    :return fidrank: Number of dimensions (rank) of fid file
    """

    # if procpar['recon'] == 'external' and fid_properties['rank'] == '3' and procpar:
    #     fid_tmp=fid_properties['roi']
    #     fid_properties['roi'][0:1] = fid_tmp[1:2]
    #     fid_properties['roi'][2] = fid_tmp[0]
    #     fid_tmp=fid_properties['dims']
    #     fid_properties['dims'][0:1] = fid_tmp[1:2]
    #     fid_properties['dims'][2] = fid_tmp[0]

    #----------------------------------------------------------
    # General implementation checks
    filename = fid_properties['filename']

    # File dimensionality or Rank fields
    # rank is a positive integer value (1, 2, 3, 4,...) giving the
    # number of dimensions in the data file (e.g., int rank=2;).
    fidrank = fid_properties['rank']
    acqndims = procpar['acqdim']
    CommentStr = 'Acquisition dimensionality (ie 2D or 3D) does not match between fid and procpar'
    AssumptionStr = 'Procpar nv2 > 0 indicates 3D acquisition and fid rank property indicates dimensionality.\n'+\
        'Using local FID value '+str(fidrank)+' instead of procpar value '+str(acqndims)+'.'
    if args.verbose:
        print 'Acqdim (type): ' + ds.MRAcquisitionType + " acqndims "  + str(acqndims)
    
    AssertImplementation(acqndims != fidrank, filename, CommentStr, AssumptionStr)
        

    # matrix is a set of rank integers giving the number of data
    # points in each dimension (e.g., for rank=2, float
    # matrix[]={256,256};)
    if fidrank==3:
        fid_size_matrix = fid_properties['dims'][0:3]
    else:
        fid_size_matrix = fid_properties['dims'][0:2]
    if args.verbose:
        print "FID size matrix ",fid_size_matrix, type(fid_size_matrix)
    #fid_size_matrix=numpy.array(fid_matrix)

    # spatial_rank is a string ("none", "voxel", "1dfov", "2dfov",
    # "3dfov") for the type of data (e.g., char
    # *spatial_rank="2dfov";).
    spatial_rank = fid_properties['spatial_rank']

    #  0018,0023 MR Acquisition Type (optional)
    # Identification of spatial data encoding scheme.
    # Defined Terms: 1D 2D 3D
    fid_MRAcquisitionType = '2D'
    if spatial_rank == "3dfov":
       fid_MRAcquisitionType = '3D'
    CommentStr = 'MR Acquisition type does not match between fid and procpar'
    AssumptionStr = 'In fid, MR Acquisition type defined by spatial_rank and matrix. '+\
        'For 2D, spatial_rank="2dfov" and matrix has two elements eg. {256,256}. '+\
        'For 3D, spatial_rank="3dfov" and matrix has three elements.\n'+\
        'In procpar, MR Acquisition type is defined by nv2 > 0 or lpe2 > 0.\n'+\
        'Using local FID value '+fid_MRAcquisitionType+' instead of procpar value '+ds.MRAcquisitionType+'.'
    AssertImplementation( ds.MRAcquisitionType != fid_MRAcquisitionType,  filename, CommentStr, AssumptionStr)
    ds.MRAcquisitionType = fid_MRAcquisitionType

    # Data Content Fields
    # The following entries define the data type and size.
    #  - storage is a string ("integer", "float") that defines the data type (e.g., char
    # *storage="float";).
    #  - bits is an integer (8, 16, 32, or 64) that defines the size of the data (e.g.,
    # float bits=32;).
    # - type is a string ("real", "imag", "absval", "complex") that defines the
    # numerical data type (e.g., char *type="absval";).

    # roi is the size of the acquired data volume (three floating
    # point values), in centimeters, in the user's coordinate frame,
    # not the magnet frame (e.g., float roi[]={10.0,15.0,0.208};). Do
    # not confuse this roi with ROIs that might be specified inside
    # the data set.
    if fidrank==3:
        roi = fid_properties['FOVcm'][0:3]
    else:
        roi = fid_properties['FOVcm'][0:2]
    if args.verbose:
        print "FID roi ",roi, type(roi)
    # roi=numpy.array(roi_text)
    
    ## PixelSpacing - 0028,0030 Pixel Spacing (mandatory)
    PixelSpacing = fid_properties['voxelmm'] 
    if PixelSpacing[0] != ds.PixelSpacing[0] or PixelSpacing[1] != ds.PixelSpacing[1]:
        print "Pixel spacing mismatch, procpar ", ds.PixelSpacing , " fid spacing ", str(PixelSpacing[0]),',', str(PixelSpacing[1])
    if args.verbose:
        print "Pixel Spacing : Procpar   ", ds.PixelSpacing
        print "Pixel Spacing : FID props ", PixelSpacing 
    ds.PixelSpacing =  [str(PixelSpacing[0]),str(PixelSpacing[1])]        #(0028,0030) Pixel Spacing


    ## FID slice thickness
    if fidrank == 3:
        fidthk = fid_properties['FOVcm'][2]/fid_properties['dims'][2]*10
    else:
        fidthk = fid_properties['FOVcm'][2]*10.0

    CommentStr = 'Slice thickness does not match between fid and procpar'
    AssumptionStr = 'In fid, slice thickness defined by roi[2] for 2D or roi[2]/matrix[2].\n'+\
        'In procpar, slice thickness defined by thk (2D) or lpe2*10/(fn2/2) or lpe2*10/nv2.\n'+\
        'Using local FID value '+str(fidthk)+' instead of procpar value '+str(ds.SliceThickness)+'.'
    if args.verbose:
        print 'fidthk : ' + str(fidthk)
        print 'SliceThinkness: ' + str(ds.SliceThickness)

    SliceThickness= float(ds.SliceThickness)

    #fix me Quick hack to avoid assert errors for diffusion and 3D magnitude images
    #if not ('diff' in procpar.keys() and procpar["diff"] == 'y'): 
    #	 if MRAcquisitionType == '3D':
    #	     print 'Not testing slicethickness in diffusion and 3D MR FIDs'
    #	else:
    AssertImplementation(SliceThickness != fidthk, filename, CommentStr, AssumptionStr)




    # Slice Thickness 0018,0050 Slice Thickness (optional)
    if fidrank == 3:
        if len(PixelSpacing) != 3:
            print "Slice thickness: 3D procpar spacing not available, fidthk ", fidthk
    else:
        if PixelSpacing[2] != ds.SliceThickness:
            print "Slice Thickness mismatch, procpar ", ds.SliceThickness , " fid spacing ", PixelSpacing[2], fidthk

    # Force slice thickness to be from fid props
    ds.SliceThickness = str(fidthk)  
    SliceThickness = fidthk


    #---------------------------------------------------------------------------------
    # GROUP 0020: Relationship

    ds.ImageComments = FID2DCM_Image_Comments+'\n'+fid_properties['filetext']
    
    # For further information regarding the location, orientation, roi, span, etc 
    # properties in the FID header, see the "Agilent VNMRJ 3.2 User Programming 
    # User Guide", pgs 434-436.  Also see VNMRJ Programming.pdf Ch5 Data Location
    # and Orientation Fields p 292.
    
    # Orientation defines the user frame of reference, and is defined according to the 
    # magnet frame of reference (X,Y,Z), where
    #	Z is along the bore, from cable end to sample end
    #	Y is bottom to top, and
    #	X is right to left, looking along positive Z
    # ref: "Agilent VnmrJ 3 Imaging User Guide" pg 679
    #
    # Location defines the position of the centre of the acquired data volume, 
    # relative to the magnet centre, in the user frame of reference.
    #
    # ROI is the size of the acquired data volume in cm in the user frame of reference.
    #
    # Origin is the coordinates of the first point in the data set, in the user frame 
    # of reference.
    #
    # 'abscissa' is a set of rank strings ("hz", "s", "cm", "cm/s",
    # "cm/s2", "deg", "ppm1", "ppm2", "ppm3") that identifies the
    # units that apply to each dimension (e.g., char
    # *abscissa[]={"cm","cm"};).
    #
    # 'span' is a set of rank floating point values for the signed
    # length of each axis, in user units. A positive value means the
    # value of the particular coordinate increases going away from the
    # first point (e.g., float span[]={10.000,-15.000};).
    # 
    # ordinate is a string ("intensity", "s", "deg") that gives the units
    # that apply to the numbers in the binary part of the file (e.g.,char
    # *ordinate[]={"intensity"};).

    
    
    orientation = numpy.matrix(fid_properties['orientation']).reshape(3,3)
    location = numpy.matrix(fid_properties['location'])*10
    span = numpy.matrix(numpy.append(fid_properties['span'], 0)*10.0)

    if args.verbose:
        print "Span: ", span, span.shape
        print "Location: ", location, location.shape
    
    # diff = numpy.setdiff1d(span, location)

    if (numpy.prod(span.shape) != numpy.prod(location.shape)):
        span=numpy.resize(span,(1,3))
    # print span
    o = location - span/2.0
    
    FirstVoxel = orientation.transpose() * o.transpose()
    
    # DICOM patient coordinate system is defined such that x increases towards the
    # patient's left, y increases towards the patient's posterior, and z increases 
    # towards the patient's head. If we imageine a (miniature) human lying supine, 
    # with their head towards the cable end of the magnet, then x in the user
    # reference frame remains the same, while y and z are inverted.
    # See DICOM Standard section C.7.6.2.1.1
    
    ImagePositionPatient = FirstVoxel.flatten().tolist()[0]
    ImagePositionPatient[1] *= -1
    ImagePositionPatient[2] *= -1
    
    ImageOrientationPatient = orientation.flatten().tolist()[0]
    ImageOrientationPatient[1] *= -1
    ImageOrientationPatient[2] *= -1
    ImageOrientationPatient[4] *= -1
    ImageOrientationPatient[5] *= -1


    # (0020,0032) Image Patient Position
    ds.ImagePositionPatient = [str(ImagePositionPatient[0]),\
				   str(ImagePositionPatient[1]),\
				   str(ImagePositionPatient[2])]	      

    #(0020,0037) Image Patient Orientation
    ds.ImageOrientationPatient = [str(ImageOrientationPatient[0]),\
				      str(ImageOrientationPatient[1]),\
				      str(ImageOrientationPatient[2]),\
				      str(ImageOrientationPatient[3]),\
				      str(ImageOrientationPatient[4]),\
				      str(ImageOrientationPatient[5])]			
    if fidrank == 3:
        # Prepare to fix 3rd dimension position using transformation matrix in Save3DFIDtoDicom
        ImageTransformationMatrix = \
            numpy.matrix([[PixelSpacing[0] * ImageOrientationPatient[0], 
                           PixelSpacing[1] * ImageOrientationPatient[1], 
                           SliceThickness * ImageOrientationPatient[2],  
                           ImagePositionPatient[0]],
                          [PixelSpacing[0] * ImageOrientationPatient[3], 
                           PixelSpacing[1] * ImageOrientationPatient[4], 
                           SliceThickness * ImageOrientationPatient[5],  
                           ImagePositionPatient[1]],
                          [PixelSpacing[0] * ImageOrientationPatient[6], 
                           PixelSpacing[1] * ImageOrientationPatient[7], 
                           SliceThickness * ImageOrientationPatient[8],  
                           ImagePositionPatient[2]], 
                          [0,0,0,1]])
    else:
        ImageTransformationMatrix=[]




    # Nuclear Data Fields
    # Data fields may contain data generated by interactions between more than one nucleus
    # (e.g., a 2D chemical shift correlation map between protons and carbon). Such data requires
    # interpreting the term ppm for the specific nucleus, if ppm to frequency conversions are
    # necessary, and properly labeling axes arising from different nuclei. To properly interpret
    # ppm and label axes, the identity of the nucleus in question and the corresponding nuclear
    # resonance frequency are needed. These fields are related to the abscissa values
    # "ppm1", "ppm2", and "ppm3" in that the 1, 2, and 3 are indices into the nucleus and
    # nucfreq fields. That is, the nucleus for the axis with abscissa string "ppm1" is the
    # first entry in the nucleus field.
    #   - nucleus is one entry ("H1", "F19", same as VNMR tn parameter) for each rf
    #       channel (e.g., char *nucleus[]={"H1","H1"};).
    #   - nucfreq is the nuclear frequency (floating point) used for each rf channel (e.g.,
    #       float nucfreq[]={200.067,200.067};).

    if fid_properties['nucleus'][0] != ds.ImagedNucleus:
        print 'Imaged nucleus mismatch: ', fid_properties['nucleus'], ds.ImagedNucleus
    if math.fabs(fid_properties['nucfreq'][0] - float(ds.ImagingFrequency)) > 0.01:
        print 'Imaging frequency mismatch: ', fid_properties['nucfreq'], ds.ImagingFrequency




    # Change patient position and orientation in 
    # if procpar['recon'] == 'external' and fid_properties['rank'] == '3':
    

    #---------------------------------------------------------------------------------
    # GROUP 0028: Image Presentation
    # A good short description of this section can be found here: 
    # http://dicomiseasy.blogspot.com.au/2012/08/chapter-12-pixel-data.html


        ## Implementation check
    CommentStr = 'Number of rows does not match between fid and procpar'
    AssumptionStr = 'In FID, number of rows is defined by matrix[1]. \n'+\
        'In procpar, for 3D datasets number of rows is either fn1/2 or nv ('+str(procpar['fn1']/2.0)+','+str(procpar['nv'])+').\n'+\
        'For 2D datasets, number of rows is fn/2.0 or np ('+str(procpar['fn']/2.0)+','+str(procpar['np'])+').\n'+\
        'Using local FID value '+str(fid_properties['dims'][1])+' instead of procpar value '+str(ds.Rows)+'.'
    AssertImplementation(int(float(ds.Rows)) != int(fid_properties['dims'][1]), filename, CommentStr, AssumptionStr)
    if args.verbose:
        print 'Rows ', procpar['fn']/2.0, procpar['fn1']/2.0, procpar['nv'], procpar['np']/2.0
        print '   Procpar: rows ', ds.Rows
        print '   FID prop rows ', fid_properties['dims'][1]
    ds.Rows = fid_properties['dims'][1]                                   #(0028,0010) Rows
            


    ## Implementation check
    CommentStr = 'Number of columns does not match between fid and procpar'
    AssumptionStr = 'In FID, number of columns is defined by matrix[0]. \n'+\
        'In procpar, for 3D datasets number of columns is either fn/2 or np ('+str(procpar['fn']/2.0)+','+str(procpar['np'])+').\n'+\
        'For 2D datasets, number of rows is fn1/2.0 or nv ('+str(procpar['fn1']/2.0)+','+str(procpar['nv'])+').\n'+\
        'Using local FID value '+str(fid_properties['dims'][0])+' instead of procpar value '+str(ds.Columns)+'.'
    AssertImplementation(int(float(ds.Columns)) != int(fid_properties['dims'][0]), filename, CommentStr, AssumptionStr)
    if args.verbose:
        print 'Columns ', procpar['fn']/2.0, procpar['fn1']/2.0, procpar['nv'], procpar['np']/2.0, fid_properties['dims'][0]
        print '   Procpar: Cols ', ds.Rows
        print '   FID prop Cols ', fid_properties['dims'][0]
    ds.Columns = fid_properties['dims'][0]                                 #(0028,0011) Columns


    #---------------------------------------------------------------------------------
    # Number of frames	      
    # DICOMHDR code:
    #	  elseif $tag='(0028,0008)' then	" no of frames "
    #	    $dim = 2  "default 2D"
    #	    exists('nv2','parameter'):$ex     
    #	    if($ex > 0) then
    #	      if(nv2 > 0) then
    #		on('fn2'):$on		"3D data"
    #		if($on) then
    #		  $pe2 = fn2/2.0
    #		else
    #		  $pe2 = nv2
    #		endif	      
    #		$dim = 3
    #	      endif
    #	    endif
    #	    
    #	    if ($dim = 3) then
    #	      $f = $pe2	   "no of frames for 3D"
    #	    else 
    #	      substr(seqcon,3,1):$spe1	  
    #	      if($spe1 = 's') then
    #		$f = (ns * (arraydim/nv) * ne)	 "sems type"
    #	      else
    #		$f = (ns * arraydim * ne)	"compressed gems type"
    #	      endif	       
    #	      if($imagesout='single') then  
    #		$f = $f	"single image output: frames=(no_of_slices * array_size * ne)"
    #	      else
    #		$f = 1				" single frame"
    #	      endif
    #	    endif     
    #	   $fs='' 
    #	    format($f,0,0):$fs
    #	    $value='['+$fs+']'	    
    #	    if $DEBUG then write('alpha','    new value = "%s"',$value) endif



    #if fidrank == 3:
    #	 ds.NumberOfFrames = fid_properties['dims'][2]
    
    # dicom3tool uses frames to create enhanced MR
    # ds.NumberOfFrames = fid_properties['slices']
    # ds.FrameAcquisitionNumber = fid_properties['slice_no']

#    if 'ne' in procpar.keys() and procpar['ne'] > 1:
#	 print 'Processing multi-echo sequence image'


    if 'echo_no' in fid_properties.keys():
        volume=fid_properties['echo_no']

    if ds.ImageType[2]=="MULTIECHO" and fid_properties['echoes'] > 1:
        print 'Multi-echo sequence'
	# TE 0018,0081 Echo Time (in ms) (optional)
        if 'TE' in fid_properties.keys():
            if ds.EchoTime != str(fid_properties['TE']):
                print "Echo Time mismatch: ",ds.EchoTime, fid_properties['TE']
            ds.EchoTime	 = str(fid_properties['TE']) 
	# 0018,0086 Echo Number (optional)
    if 'echo_no' in fid_properties.keys():
        if ds.EchoNumber != fid_properties['echo_no']:
            print "Echo Number mismatch: ",ds.EchoNumber, fid_properties['echo_no']
        ds.EchoNumber = fid_properties['echo_no']		 


    if ds.ImageType[2] == "ASL":	  
        ds=ParseASL(ds,procpar,fid_properties)

    #if 'echoes' in fid_properties.keys() and fid_properties['echoes'] > 1 and fid_properties['array_dim'] == 1:
    #    ds.AcquisitionNumber = fid_properties['echo_no']
    #    ds.ImagesInAcquisition = fid_properties['echoes']
    #else:

    ds.AcquisitionNumber = fid_properties['array_index']
    ds.ImagesInAcquisition = fid_properties['array_dim']


    if ds.ImageType[2] == 'DIFFUSION':
        ds=ParseDiffusionFID(ds,procpar,fid_properties,args)


    ## Multi dimension Organisation and Index module
    DimOrgSeq = Dataset()
    #ds.add_new((0x0020,0x9164), 'UI', DimensionOrganizationUID)

    if ds.ImageType[2] == "MULTIECHO": # or SEQUENCE == "Diffusion":
        DimensionOrganizationUID = [ProcparToDicomMap.CreateUID(UID_Type_DimensionIndex1,[],[],args.verbose), ProcparToDicomMap.CreateUID(UID_Type_DimensionIndex2,[],[],args.verbose)]
        DimOrgSeq.add_new((0x0020,0x9164), 'UI',DimensionOrganizationUID)
        ds.DimensionOrganizationType='3D_TEMPORAL'  #or 3D_TEMPORAL
    else:
        DimensionOrganizationUID = ProcparToDicomMap.CreateUID(UID_Type_DimensionIndex1,[],[],args.verbose)
        #if args.verbose:
        #    print "DimUID", DimensionOrganizationUID
        DimOrgSeq.add_new((0x0020,0x9164), 'UI',[DimensionOrganizationUID])
        ds.DimensionOrganizationType='3D'  #or 3D_TEMPORAL
        
    ds.DimensionOrganizationSequence= Sequence([DimOrgSeq])


    if ds.ImageType[2] == 'MULTIECHO':
        DimIndexSeq1 = Dataset()
        DimIndexSeq1.DimensionIndexPointer = (0x0020,0x0032)  # Image position patient 20,32 or 20,12 

        # #DimIndexSeq1.DimensionIndexPrivateCreator=
        # #DimIndexSeq1.FunctionalGroupPointer=
        # #DimIndexSeq1.FunctionalGroupPrivateCreator=
        DimIndexSeq1.add_new((0x0020,0x9164), 'UI',DimOrgSeq.DimensionOrganizationUID[0])
        DimIndexSeq1.DimensionDescriptionLabel='Third Spatial dimension'

        DimIndexSeq2 = Dataset()
        DimIndexSeq2.DimensionIndexPointer=(0x0018,0x0081)  # Echo Time
        # DimIndexSeq2.DimensionIndexPrivateCreator=
        # DimIndexSeq2.FunctionalGroupPointer=
        # DimIndexSeq2.FunctionalGroupPrivateCreator=
        DimIndexSeq2.add_new((0x0020,0x9164), 'UI',DimOrgSeq.DimensionOrganizationUID[1])
        DimIndexSeq2.DimensionDescriptionLabel='Fourth dimension (multiecho)'
        ds.DimensionIndexSequence = Sequence([DimIndexSeq2, DimIndexSeq1])
    else:
        DimIndexSeq1 = Dataset()
        DimIndexSeq1.DimensionIndexPointer = (0x0020,0x0032)  # Image position patient 20,32 or 20,12 
        # #DimIndexSeq1.DimensionIndexPrivateCreator=
        # #DimIndexSeq1.FunctionalGroupPointer=
        # #DimIndexSeq1.FunctionalGroupPrivateCreator=
        DimIndexSeq1.add_new((0x0020,0x9164), 'UI',[DimensionOrganizationUID])
        DimIndexSeq1.DimensionDescriptionLabel='Third Spatial dimension'
        ds.DimensionIndexSequence = Sequence([DimIndexSeq1])
    

        # Module: Image Pixel (mandatory)
        # Reference: DICOM Part 3: Information Object Definitions C.7.6.3 
        # ds.Rows                                                              # 0028,0010 Rows (mandatory)
        # ds.Columns                                                           # 0028,0011 Columns (mandatory)
        # ds.BitsStored                                                        # 0028,0101 (mandatory)
        # ds.HighBit                                                           # 0028,0102 (mandatory)
        # ds.PixelRepresentation                                               # 0028,0103 Pixel Representation (mandatory)
        # ds.PixelData                                                        # 7fe0,0010 Pixel Data (mandatory)



    FrameContentSequence=Dataset()    
    #FrameContentSequence.FrameAcquisitionNumber = '1' #fid_properties['slice_no']
    #FrameContentSequence.FrameReferenceDateTime
    #FrameContentSequence.FrameAcquisitionDateTime
    #FrameContentSequence.FrameAcquisitionDuration
    #FrameContentSequence.CardiacCyclePosition
    #FrameContentSequence.RespiratoryCyclePosition
    #FrameContentSequence.DimensionIndexValues = 1 #islice # fid_properties['array_no']
    #FrameContentSequence.TemporalPositionIndex = 1
    FrameContentSequence.StackID = [str(1)] # fourthdimid
    FrameContentSequence.InStackPositionNumber = [int(1)] #fourthdimindex
    FrameContentSequence.FrameComments=  fid_properties['filetext']
    FrameContentSequence.FrameLabel='DimX'
    ds.FrameContentSequence = Sequence([FrameContentSequence])

    
    return ds,fidrank,fid_size_matrix,ImageTransformationMatrix

def Save3dFIDtoDicom(ds,procpar,image_data,fid_properties,M,args,outdir,filename):
    """
    Multi-dimension (3D) export of FID format image data and metadata to DICOM

    :param ds:  Baseline pydicom dataset
    :param image_data: Pixel data of image
    :param fid_matsize: Matrix size of image data e.g. {256,256,256}
    :param M: Image transformationm matrix
    :param args: argparse struct
    :param outdir: output directory
    :return ds: Return the pydicom dataset
    """

    print "3D export"
    voldata = numpy.reshape(image_data,fid_properties['dims'])

    # if procpar['recon'] == 'external':
    # 
    #        pdb.set_trace()
    if procpar['recon'] == 'external' and fid_properties['rank'] == 3:
        if procpar['seqfil'] == "epip":
            print "Transposing external recon 3D"
            voldata = numpy.transpose(voldata,(1,2,0)) # 1,2,0
        if procpar['seqfil'] == "fse3d":
            print "Transposing external recon 3D"
            voldata = numpy.transpose(voldata,(2,0,1)) # 0,2,1 works
    # readpp.m procpar('nD') == 3            
    #        acq.FOVcm = [pps.lro pps.lpe pps.lpe2];
    #        acq.dims = [pps.nf pps.np/2 pps.nv2];
    #        acq.voxelmm = acq.FOVcm./acq.dims*10;

    print "Image data shape: ", str(image_data.shape)
    print "Vol data shape: ", voldata.shape
    print "fid properties matrix: ", fid_properties['dims']
    print "Slice points: ", fid_properties['dims'][0]*fid_properties['dims'][1]
    #  slice_data = numpy.zeros_like(numpy.squeeze(image_data[:,:,1]))
    #   if 'ne' in procpar.keys():
        
    range_max = fid_properties['dims'][2]
    num_slicepts = fid_properties['dims'][0]*fid_properties['dims'][1]
    if procpar['recon'] == 'external' and procpar['seqfil'] == 'fse3d' and fid_properties['rank'] == 3: 
        range_max = fid_properties['dims'][1]
        num_slicepts = fid_properties['dims'][0]*fid_properties['dims'][2]
        ds.Columns = fid_properties['dims'][2]
        ds.Rows = fid_properties['dims'][0]
        ##FIXME FSE3d still producing bad dicoms

    if args.verbose:
        print "Columns ", ds.Columns, " Rows ", ds.Rows
        print "Range max and no slice points: ", range_max, num_slicepts
        print "Voldata[1] shape: ", voldata[:,:,0].shape

    ## Indexing in numpy matrix begins at 0, fid/dicom filenames begin at 1
    for islice in xrange(0,range_max):    
        # Reshape volume slice to 1D array
        slice_data = numpy.reshape(voldata[:,:,islice],(num_slicepts,1)) 
            # Convert Pixel data to string
        ds.PixelData = slice_data.tostring()    #(7fe0,0010) Pixel Data

        #if acqndims == 3:
        if 'slice_no' in fid_properties.keys():
            image_number = fid_properties['slice_no']
        else:
            image_number=int(re.sub(r'^.*image(\d+).*', r'\1',filename))

        new_filename = "slice%03dimage%03decho%03d.dcm" % (islice+1,image_number,fid_properties['echo_no'])

        if procpar['recon'] == 'external' and fid_properties['rank'] == 3 and procpar['seqfil'] == 'fse3d':
            pos = numpy.matrix([[0],[0],[islice],[1]])
        else:
            pos = numpy.matrix([[0],[0],[islice],[1]])

        Pxyz = M * pos
        ds.ImagePositionPatient= [str(Pxyz[0,0]),str(Pxyz[1,0]),str(Pxyz[2,0])]

        # ds.FrameContentSequence[0].StackID = [str(volume)] # fourthdimid
        ds.FrameContentSequence[0].InStackPositionNumber = [int(islice)] #fourthdimindex
        ds.FrameContentSequence[0].TemporalPositionIndex = ds.EchoNumber
        # ds.InStackPosition = islice #str(islice)

        # Save DICOM
        ds.save_as(os.path.join(outdir, new_filename))

        return ds

def Save2dFIDtoDicom(ds,image_data,outdir, filename):
    """
    Export 2D image and metadata to DICOM

    """
    
    # Common export format
    ds.PixelData = image_data.tostring()         # (7fe0,0010) Pixel Data

    if ds.ImageType[2] == "ASL":
        image_number=fid_properties['array_index']
        if fid_properties["asltag"] == 1:               # Labelled
            new_filename = "slice%03dimage%03decho%03d.dcm" % (fid_properties['slice_no'],image_number,1)
        elif  fid_properties["asltag"] == -1:            #  Control
            new_filename = "slice%03dimage%03decho%03d.dcm" % (fid_properties['slice_no'],image_number,2)
        else:                                            # Unknown
            new_filename = "slice%03dimage%03decho%03d.dcm" % (fid_properties['slice_no'],image_number,3)
        dcmfilename=os.path.join(outdir,new_filename)
    else:
        dcmfilename=os.path.join(outdir, os.path.splitext(filename)[0] + '.dcm')
    # Save DICOM
    ds.save_as(dcmfilename)
    
    return 1

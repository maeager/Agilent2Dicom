#!/usr/bin/env python
"""ReadFID is used to read Agilent FID image files

   - Michael Eager (michael.eager@monash.edu)
"""
"""
  Copyright (C) 2014 Michael Eager

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

import os
import sys
import math
import re
import struct
import numpy as np
import argparse
import ReadProcpar
import ProcparToDicomMap
from scipy.fftpack import ifftn, fftshift, ifftshift
import scipy.io
import agilent2dicom_globalvars as A2D
import RescaleFDF
# from scipy import ndimage
from dicom.sequence import Sequence
from dicom.dataset import Dataset


def get_bit(value, bit_number):
    """ get_bit read binary value at position bit_number
    """
    return (value & int(1 << (bit_number - 1))) != 0
# end get_bit


def readfid(fidfolder, procpar, args):
    """readfid
     read kspace data from Agilent fid and procpar files.
     output img has dimensions [freq, phase, slice, channel, echo]
     note that the image may have to be circular shifted in the phase
     direction.

    >>>  fid_header, realpart, imagpart = readfid(fidfolder, procpar,args)
    """

    #  get acqcycles and TE from procpar
    if not procpar:
        procpar, procpartext = ReadProcpar.ReadProcpar(
            os.path.join(fidfolder, 'procpar'))

    # Define fid headers from procpar on first occasion
    fid_header = dict()
    fid_header['TE'] = procpar['te']
    # FOV
    fid_header['volumes'] = int(procpar['lpe3'])
    fid_header['nEchoes'] = int(procpar['ne'])
    rcvrs = re.findall('y', procpar['rcvrs'])
    if rcvrs:
        fid_header['nChannels'] = len(rcvrs)
    else:
        fid_header['nChannels'] = 1
    fid_header['mode'] = '%dD' % procpar['nD']
    if procpar['nD'] == 2:
        fid_header['FOVcm'] = [procpar['lro'], procpar['lpe']]
        if 'diff' in procpar.keys() and procpar['diff'] == 'y':
            fid_header['dims'] = [
                procpar['fn'] / 2,
                procpar['fn1'] / 2,
                procpar['ns']]
        else:
            fid_header['dims'] = [
                procpar['nf'] / procpar['ns'],
                procpar['np'] / 2,
                procpar['ns']]
        # if len(procpar['thk']) > 1:
        #    print "procpar thk size greater than 1"
        fid_header['voxelmm'] = np.array(
            [procpar['lro'] / fid_header['dims'][0],
             procpar['lpe'] / fid_header['dims'][1],
             procpar['thk']]) * 10
    elif procpar['nD'] == 3:
        fid_header['FOVcm'] = [procpar['lro'],
                               procpar['lpe'],
                               procpar['lpe2']]
        fid_header['dims'] = [procpar['nf'] / procpar['ne'],
                              procpar['np'] / 2,
                              procpar['nv2']]
        fid_header['voxelmm'] = np.array(
            fid_header['FOVcm']) / np.array(fid_header['dims']) * 10
    else:
        print 'Unknown nD', procpar['nD']
    print fid_header

    # open fid file
    f = open(os.path.join(fidfolder, 'fid'), "rb")
    int16size = struct.calcsize('h')
    int32size = struct.calcsize('i')
    endian = '>'  # > for big-endian < for little
    # Read datafileheader using: x, = struct.unpack(type,binary) method
    # unpack returns a tuple and we only want the result  (p.285-286)
    fid_header['nblocks'], = struct.unpack(
        endian + 'i', f.read(int32size))
    #  nblocks is the number of data blocks present in the file. ('int32')
    fid_header['ntraces'], = struct.unpack(
        endian + 'i', f.read(int32size))
    # ntraces is the number of traces in each block. ('int32')
    fid_header['np'], = struct.unpack(
        endian + 'i', f.read(int32size))
    # np is the number of simple elements (16-bit integers, 32-bit
    # integers, or 32-bit floating point numbers) in one trace. It is
    # equal to twice the number of complex data points. ('int32')
    fid_header['ebytes'], = struct.unpack(
        endian + 'i', f.read(int32size))
    #  ebytes is the number of bytes in one element, either 2 (for
    #  16-bit integers in single precision FIDs) or 4 (for all
    #  others). ('int32')
    fid_header['tbytes'], = struct.unpack(
        endian + 'i', f.read(int32size))
    # tbytes is set to (np*ebytes).('int32')
    fid_header['bbytes'], = struct.unpack(
        endian + 'i', f.read(int32size))
    # bbytes is set to (ntraces*tbytes + nbheaders*sizeof(struct
    # datablockhead)). The size of the datablockhead structure is 28
    # bytes. ('int32')

    fid_header['vers_id'], = struct.unpack(
        endian + 'h', f.read(int16size))  # vers_id is the version identification of present VNMR. ('int16')
    status = f.read(int16size)
    fid_header['status'], = struct.unpack(
        endian + 'h', status)  # see below ('int16')
    fid_header['nbheaders'], = struct.unpack(
        endian + 'i', f.read(int32size))  # nbheaders is the number of block headers per data block. ('int32')
    if args.verbose:
        print 'status : ', fid_header['status'], type(fid_header['status']), type(status)

  #       status is bits as defined below with their hexadecimal values.
  # All other bits must be zero.
  # Bits 0-6: file header and block header status bits (bit 6 is unused):
  #  0         S_DATA                   0x1          0 = no data, 1 = data
  #  1         S_SPEC                   0x2          0 = FID, 1 = spectrum
  #  2         S_32                     0x4          *
  #  3         S_FLOAT                  0x8          0 = integer, 1 = floating point
  #  4         S_COMPLEX                0x10         0 = real, 1 = complex
  #  5         S_HYPERCOMPLEX           0x20         1 = hypercomplex
  #
  #* If S_FLOAT=0, S_32=0 for 16-bit integer, or S_32=1 for 32-bit integer.
  #  If S_FLOAT=1, S_32 is ignored.

    fid_header['s_data'] = int(get_bit(fid_header['status'], 1))
    fid_header['s_spec'] = int(get_bit(fid_header['status'], 2))
    fid_header['s_32'] = int(get_bit(fid_header['status'], 3))
    fid_header['s_float'] = int(get_bit(fid_header['status'], 4))
    fid_header['s_complex'] = int(get_bit(fid_header['status'], 5))
    fid_header['s_hyper'] = int(get_bit(fid_header['status'], 6))

  #  Bits 7-14: file header status bits (bits 10 and 15 are unused):
  # 7        S_ACQPAR            0x80       0 = not Acqpar, 1 = Acqpar
  # 8        S_SECND             0x100      0 = first FT, 1 = second FT
  # 9        S_TRANSF            0x200      0 = regular, 1 = transposed
  # 11       S_NP                0x800      1 = np dimension is active
  # 12       S_NF                0x1000     1 = nf dimension is active
  # 13       S_NI                0x2000     1 = ni dimension is active
  # 14       S_NI2               0x4000     1 = ni2 dimension is active
    fid_header['s_acqpar'] = int(get_bit(fid_header['status'], 8))
    fid_header['s_secnd'] = int(get_bit(fid_header['status'], 9))
    fid_header['s_transf'] = int(get_bit(fid_header['status'], 10))
    fid_header['s_np'] = int(get_bit(fid_header['status'], 12))
    fid_header['s_nf'] = int(get_bit(fid_header['status'], 13))
    fid_header['s_ni'] = int(get_bit(fid_header['status'], 14))
    fid_header['s_ni2'] = int(get_bit(fid_header['status'], 15))

    if args.verbose:
        # ['s_data'], fid_header['s_spec'], fid_header['s_32'], fid_header['s_float'], fid_header['s_complex'], fid_header['s_hyper']
        print fid_header

    if fid_header['s_float'] == 1:
        # data = fread(fid, fid_header.np*fid_header.ntraces, '*float32')
        dtype_str = np.dtype(endian + 'f4')  # 'float32'
        if args.verbose:
            print 'reading 32bit floats '
    elif fid_header['s_32'] == 1:
        # data = fread(fid, fid_header.np*fid_header.ntraces, '*int32')
        dtype_str = 'int32'
        if args.verbose:
            print 'reading 32bit int'
    else:
        # data = fread(fid, fid_header.np*fid_header.ntraces, '*int16')
        dtype_str = 'int16'
        if args.verbose:
            print 'reading 16bit int'

    dims = [0, 0, 0]
    # validate dimensions
    if procpar['nD'] == 2:
        if 'diff' in procpar.keys() and procpar['diff'] == 'y':
            dims[0] = procpar['fn'] / 2,
            dims[1] = procpar['fn1'] / 2,
            dims[2] = procpar['ns']
        else:
            dims[0] = procpar['np'] / 2        # num phase encode lines / 2
            # num frequency lines acquired / # echoes
            dims[1] = procpar['nf'] / procpar['ns']
            dims[2] = procpar['ns']          # if 2D, num slices, else ni2
        # fse3d sequence has nD == 2, but is a 3d acquisition???
        if procpar['ni2'] > 1:
            dims[2] = procpar['ni2']
    elif procpar['nD'] == 3:
        dims[0] = procpar['np'] / 2        # NUM phase encode lines / 2
        # num frequency lines acquired / # echoes
        dims[1] = procpar['nf'] / procpar['ne']
        dims[2] = procpar['ni2']
    else:
        raise ValueError(
            "Can only handle 2D or 3D files (based on procpar field nD)")
    if args.verbose:
        print 'Dimensions: ', dims, fid_header['dims']

    if fid_header['np'] != int(procpar['np']) or \
       fid_header['ntraces'] != int(procpar['nf']) or \
       fid_header['nblocks'] != int(procpar['arraydim']):
        if args.verbose:
            print 'NP ', fid_header['np'], procpar['np'], \
                ' NF ', fid_header['ntraces'], procpar['nf'], \
                ' Blocks ', fid_header['nblocks'], procpar['arraydim']
        raise ValueError(
            "Cannot resolve fid header with procpar. We're probably not interpreting the procpar correctly.")
    if fid_header['nChannels'] != int(fid_header['nblocks'] / procpar['acqcycles']):
        print "ReadFID error nChannels %d %f " % (fid_header['nChannels'], float(fid_header['nblocks'] / procpar['acqcycles']))
        fid_header['nChannels'] = int(
            fid_header['nblocks'] / procpar['acqcycles'])
    fid_header['nPhaseEncodes'] = fid_header['ntraces'] / fid_header['nEchoes']
    fid_header['rank'] = procpar['nD']

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
    if args.verbose:
        print "Shifts of slices along z axis ", procpar['pss'], procpar['pss0']
    # Create FDF-like header variables
    fid_header['location'] = [-procpar['pro'], procpar['ppe'], procpar['pss0']]
    fid_header['span'] = [procpar['lro'], procpar['lpe'], procpar['lpe2']]
    fid_header['roi'] = [procpar['lro'], procpar['lpe'], procpar['lpe2']]
    fid_header['origin'] = np.array(
        fid_header['location']) - np.array(fid_header['span']) / 2.0
    # TODO orientation alignment still needs fixing - what does this tag do to
    # the data?
    if procpar['orient'] == "sag":
        # TODO fid_header['orientation']= [0,0,1,1,0,0,0,1,0]
        fid_header['orientation'] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
    elif procpar['orient'] == "oblique":
        fid_header['orientation'] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
    else:
        fid_header['orientation'] = [1, 0, 0, 0, 1, 0, 0, 0, 1]

    # reset output structures
    ksp_data_real = np.empty([fid_header['np'] / 2, fid_header['ntraces'],
                              fid_header['nblocks']], dtype=np.float32)
    ksp_data_imag = np.empty([fid_header['np'] / 2, fid_header['ntraces'],
                              fid_header['nblocks']], dtype=np.float32)
    if args.verbose:
        print "Raw data shape:", ksp_data_real.shape
    # We have to read data every time in order to increment file pointer
    nchar = 0
    iblock = 0
    for iblock in xrange(0, fid_header['nblocks']):
        # fprintf(1, repmat('\b', 1, nchar))
        if args.verbose:
            print 'reading block ', iblock + 1, ' of ', fid_header['nblocks']
        # Read a block header
        header = dict()
        header['scale'], = struct.unpack(
            endian + 'h', f.read(int16size))  # scaling factor fid,1,'int16')
        header['bstatus'], = struct.unpack(
            endian + 'h', f.read(int16size))  # status of data in block fid,1,'int16')
        header['index'], = struct.unpack(
            endian + 'h', f.read(int16size))  # block index  fid,1,'int16')
        header['mode'], = struct.unpack(
            endian + 'h', f.read(int16size))  # mode of data in block fid,1,'int16')
        header['ctcount'], = struct.unpack(
            endian + 'i', f.read(int32size))  # ct value for FID fid,1,'int32')
        header['lpval'], = struct.unpack(
            endian + 'f', f.read(int32size))  # f2 (2D-f1) left phase in phasefile fid,1,'float32')
        header['rpval'], = struct.unpack(
            endian + 'f', f.read(int32size))  # f2 (2D-f1) right phase in phasefile fid,1,'float32')
        header['lvl'], = struct.unpack(
            endian + 'f', f.read(int32size))  # level drift correction fid,1,'float32')
        header['tlt'], = struct.unpack(
            endian + 'f', f.read(int32size))  # tilt drift correction  fid,1,'float32')

        header['s_data'] = int(get_bit(header['bstatus'], 1))
        header['s_spec'] = int(get_bit(header['bstatus'], 2))
        header['s_32'] = int(get_bit(header['bstatus'], 3))
        header['s_float'] = int(get_bit(header['bstatus'], 4))
        header['s_complex'] = int(get_bit(header['bstatus'], 5))
        header['s_hyper'] = int(get_bit(header['bstatus'], 6))


# status is bits 0 to 6 defined the same as for file header status. Bits 7-11 are defined
# below (all other bits must be zero):
#  7                                   0x80       0 = absent, 1 = present
#           MORE_BLOCKS
#  8                                   0x100      0 = real, 1 = complex
#           NP_CMPLX
#  9                                   0x200      0 = real, 1 = complex
#           NF_CMPLX
#  10                                  0x400      0 = real, 1 = complex
#           NI_CMPLX
#  11                                  0x800      0 = real, 1 = complex
#           NI2_CMPLX
        header['s_more'] = int(get_bit(header['bstatus'], 8))
        header['s_np_cmplx'] = int(get_bit(header['bstatus'], 9))
        header['s_nf_cmplx'] = int(get_bit(header['bstatus'], 10))
        header['s_ni_cmplx'] = int(get_bit(header['bstatus'], 11))
        header['s_ni2_cpmlx'] = int(get_bit(header['bstatus'], 12))
        if header['s_hyper']:
            # short s_spare1; /* short word: spare */
            header['hyper_spare1'], = struct.unpack(
                endian + 'h', f.read(int16size))
            # short status;   /* status word for block header */
            header['hyper_status'], = struct.unpack(
                endian + 'h', f.read(int16size))
            # short s_spare2; /* short word: spare */
            header['hyper_spare2'], = struct.unpack(
                endian + 'h', f.read(int16size))
            # short s_spare3; /* short word: spare */
            header['hyper_spare3'], = struct.unpack(
                endian + 'h', f.read(int16size))
            # long l_spare1;  /* long word: spare */
            header['hyper_lspare1'], = struct.unpack(
                endian + 'i', f.read(int32size))
            # float lpval1;   /* 2D-f2 left phase */
            header['hyper_lpval1'], = struct.unpack(
                endian + 'f', f.read(int32size))
            # float rpval1;   /* 2D-f2 right phase */
            header['hyper_rpval1'], = struct.unpack(
                endian + 'f', f.read(int32size))
            # float f_spare1; /* float word: spare */
            header['hyper_f_spare1'], = struct.unpack(
                endian + 'f', f.read(int32size))
            # float f_spare2; /* float word: spare */
            header['hyper_f_spare2'], = struct.unpack(
                endian + 'f', f.read(int32size))

# Main data block header mode bits 0-15:
#    Bits 0-3: bit 3 is currently unused
#      0        NP_PHMODE                           0x1    1 = ph mode
#      1        NP_AVMODE                           0x2    1 = av mode
#      2        NP_PWRMODE                           0x4    1 = pwr mode
#    Bits 4-7: bit 7 is currently unused
#      4        NF_PHMODE                           0x10   1 = ph mode
#      5        NF_AVMODE                           0x20   1 = av mode
#      6        NF_PWRMODE                           0x40   1 = pwr mode
#    Bits 8-11: bit 11 is currently unused
#      8        NI_PHMODE                           0x100  1 = ph mode
#      9        NI_AVMODE                           0x200  1 = av mode
#      10       NI_PWRMODE                           0x400  1 = pwr mode
#    Bits 12-15: bit 15 is currently unused
#      12       NI2_PHMODE                           0x8    1 = ph mode
#      13       NI2_AVMODE                           0x100  1 = av mode
#      14       NI2_PWRMODE                           0x2000 1 = pwr mode
        header['np_phmode'] = int(get_bit(header['mode'], 1))
        header['np_avmode'] = int(get_bit(header['mode'], 2))
        header['np_pwrmode'] = int(get_bit(header['mode'], 3))
        header['nf_phmode'] = int(get_bit(header['mode'], 5))
        header['nf_avmode'] = int(get_bit(header['mode'], 6))
        header['nf_pwrmode'] = int(get_bit(header['mode'], 7))
        header['ni_phmode'] = int(get_bit(header['mode'], 9))
        header['ni_avmode'] = int(get_bit(header['mode'], 10))
        header['ni_pwrmode'] = int(get_bit(header['mode'], 11))
        header['ni2_phmode'] = int(get_bit(header['mode'], 13))
        header['ni2_avmode'] = int(get_bit(header['mode'], 14))
        header['ni2_pwrmode'] = int(get_bit(header['mode'], 15))

        if args.verbose:
            print header
        data = np.fromfile(
            f, count=fid_header['np'] * fid_header['ntraces'], dtype=dtype_str)
        if args.verbose:
            print "Raw FID data: Dim and shape: ", data.ndim, data.shape
        data = np.reshape(data, [fid_header['ntraces'], fid_header['np']])
        if args.verbose:
            print "Reshaped raw FID data: Dim and shape: ", data.ndim, data.shape, " max np ", fid_header['np']
        # fid_header['np'] #[::2, :] #
        ksp_data_real[:, :, iblock] = np.matrix(data[:, :fid_header['np']:2]).T
        # fid_header['np'] #[1::2, :]      #
        ksp_data_imag[:, :, iblock] = np.matrix(
            data[:, 1:fid_header['np']:2]).T
        # break
    if f.tell() != os.fstat(f.fileno()).st_size:
        print "ReadFID: OS Error, fid file completed without finishing to end of file."
    f.close()
    # print iblock
    if iblock == 0:
        if args.verbose:
            print "Reshaping single block data ", dims
        ksp_data_real = np.reshape(ksp_data_real, dims)
        ksp_data_imag = np.reshape(ksp_data_imag, dims)
    # fid_header.procpar = procpar
    if args.verbose:
        print "Data Row 1:   %.5g %.5g %.5g %.5g %.5g  %.5g ... %.5g %.5g %.5g %.5g" % (data[0, 0], data[0, 1], data[0, 2], data[0, 3], data[0, 4], data[0, 5], data[0, -4], data[0, -3], data[0, -2], data[0, -1])
        print "Data Col 1:   %.5g %.5g %.5g %.5g %.5g %.5g  ... %.5g %.5g %.5g %.5g" % (data[0, 0], data[1, 0], data[2, 0], data[3, 0], data[4, 0], data[5, 0], data[-4, 0], data[-3, 0], data[-2, 0], data[-1, 0])
        print "Data Row -1:   %.5g %.5g %.5g %.5g %.5g %.5g ... %.5g %.5g %.5g %.5g" % (data[-1, 0], data[-1, 1], data[-1, 2], data[-1, 3], data[-1, 4], data[-1, 5], data[-1, -4], data[-1, -3], data[-1, -2], data[-1, -1])
        print "Data Col -1:   %.5g %.5g %.5g %.5g %.5g %.5g ... %.5g %.5g %.5g %.5g" % (data[0, -1], data[1, -1], data[2, -1], data[3, -1], data[4, -1], data[5, -1], data[-4, -1], data[-3, -1], data[-2, -1], data[-1, -1])
        print "ksp_data_real : %.5g  %.5g %.5g | %.5g %.5g | %.5g %.5g |%.5g" % (ksp_data_real[0, 0, iblock], ksp_data_real[0, 1, iblock], ksp_data_real[0, 2, iblock], ksp_data_real[1, 0, iblock], ksp_data_real[2, 0, iblock], ksp_data_real[0, -1, iblock], ksp_data_real[-1, 0, iblock], ksp_data_real[-1, -1, iblock])
        print "ksp_data_imag : %.5g  %.5g %.5g | %.5g %.5g | %.5g %.5g |%.5g" % (ksp_data_imag[0, 0, iblock], ksp_data_imag[0, 1, iblock], ksp_data_imag[0, 2, iblock], ksp_data_imag[1, 0, iblock], ksp_data_imag[2, 0, iblock], ksp_data_imag[0, -1, iblock], ksp_data_imag[-1, 0, iblock], ksp_data_imag[-1, -1, iblock])
        print "Final dims and shape of ksp_data_real: ", dims, ksp_data_real.shape
        return procpar, fid_header, dims, ksp_data_real, ksp_data_imag
# end readfid


def Cardiac_ASL_recon(procpar, fid_header, dims, ksp_data_real, ksp_data_imag,  args):
    """ASL_recon modified from seg_recon_ms_complex.m script in UCL groups Cardiac ASL script


    """
    #procpar, fid_header, dims, ksp_data_real, ksp_data_imag = readfid(
    #    fidfolder, procpar, args)
    FID = np.complex(ksp_data_real, ksp_data_imag)
    ns = procpar['ns']
    frames = procpar['noframes']
    ro = procpar['lro'] * frames  # or np/2 of nv
    pe = procpar['lpe'] * frames
    seg = procpar['etl']  # locetl
    ASL = 2
    numb = ASL * frames
    # pelist = [-63,  -31,    1,   33,   -62,  -30,    2,    34,   -61,   -29,     3,  35,  -60, -28,    4,   36,  -59, -27,     5,    37,  -58,   -26,   6,	38, - 57,  -25,   7,  39, -56,  -24, 8,    40,   -55,   -23,   9,   41, -54,   -22,   10,  42,  -53,  -21,   11,   43,  -52,  -20,	12, 44,  -51, -19, 13,    45, -50, -18,    14,   46,  -49,  -17,  15,    47,   -48,   -16,              16, 48, -47, -15,  17,    49,   -46,	-14,  18,   50,   -45, -13,    19,   51, -44,  -12,  20,   52,  -43,  -11,   21,   53,   -42,   -10,  22,    54, -41,  -9,  23,  55, -40,  -8,   24,   56,  -39,   -7,    25,   57,  -38,   -6,   26, 58, -37,   -5, 27,  59, -36,   -4,    28,    60, -35,   -3, 29,	61,  -34,    -2,   30,   62,  -33,   -1,   31,  63,  -32,   0,   32,  64]
    for step in xrange(0,pe/4.0):
        #pelist[(step*4+1):(step*4+4)]= step-63:pe/4.0:pe/2.0
        pelist[(step*4 +1):(step*4+5)]= np.arange((step-int(pe/2.0-1)),int(pe/2.0)+1,int(pe/4.0))
    #Check pelist with procpar
    if 'pelist' in procpar.keys() and len(procpar['pelist']) == pe:
        pediff = np.sum(np.array(pelist)-np.array(procpar['pelist']))
        if pediff != 0:
            print 'Procpar pelist does not match ASL recon list'
            print 'Using calculated version instead of procpar.'
        
    fidtmp = np.reshape(FID, [ro, seg, ns, pe / seg * numb])
    ksp = np.zeros(ro, pe, numb, ns)
    for ll in xrange(1, ns):
        fidtmp2 = np.squeeze(fidtmp[:, :, ll, :])
        for j in xrange(1, numb):
            fidtmp3[:, :, j] = np.reshape(fidtmp2[:, :, j::numb], [ro, pe])
        ksp[:, :, :, ll] = fidtmp3[:, pelist + 63, :]

#        for jj in xrange(1,numb):
    img = fftshift(ifftn(fftshift(ksp)))
    return ksp, img, procpar, fid_header


def recon(procpar, dims, fid_header, ksp_data_real, ksp_data_imag, args):
    """recon Reconstruct k-space image data into N-D image
    :param procpar:   procpar dictionary
    :param dims: dimension array
    :param fid_header: header info in fid
    :param ksp_data_real: real component of image k-space
    :param ksp_data_imag: imaginary component of image k-space
    """
    if args.verbose:
        print 'Reconstructing image'
        print dims[0], dims[1], dims[2], fid_header['nChannels'], fid_header['nEchoes']
        print 'nchannels ', fid_header['nblocks'], 'acqcycles ', procpar['acqcycles']
        print 'Image shape ', ksp_data_real.shape

    if np.product(dims) != np.product(ksp_data_real.shape):
        print "ksp not arranged properly"

    if fid_header['nChannels'] == 1 and fid_header['nEchoes'] == 1:
        ksp = np.empty(
            [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32
        img = np.empty(
            [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32
    else:  # two float32
        ksp = np.empty([dims[0], dims[1], dims[2], fid_header['nChannels'],
                        fid_header['nEchoes']], dtype=np.complex64)  # float32
        img = np.empty([dims[0], dims[1], dims[2], fid_header['nChannels'],
                        fid_header['nEchoes']], dtype=np.complex64)  # float32
    # Setup pelist
    if 'pelist' in procpar.keys():
        pelist = np.array(procpar['pelist']).astype(
            int) - int(min(procpar['pelist']))
    if procpar['nD'] == 2 and procpar['ni2'] == 1:
        print 'Reconstructing mode 1: 2D slices'
        if fid_header['nEchoes'] == 1 and fid_header['nChannels'] == 1:
            # two float32
            ksp = np.empty(
                [dims[0], dims[1], dims[2]], dtype=np.complex64)
            # two float32
            img = np.empty(
                [dims[0], dims[1], dims[2]], dtype=np.complex64)
            # [:, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
            ksp[:, :, :].real = ksp_data_real
            # [:, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
            ksp[:, :, :].imag = ksp_data_imag
            for islice in xrange(0, int(dims[2])):
                # 'pelist'' in procpar.keys() and len(procpar['pelist']) == int(dims[2]):
                if 'pelist' in locals():
                    if args.verbose:
                        print procpar['pelist']
                    img[:, :, islice] = fftshift(ifftn(ifftshift(ksp[:,
                                                                     pelist,
                                                                     islice])))
                else:
                    img[:, :, islice] = fftshift(
                        ifftn(ifftshift(ksp[:, :, islice])))
        elif fid_header['nChannels'] == 1:
            print 'Reconstructing mode 2: 2D slices with ni2=' + procpar['ni2']
            for echo in xrange(0, int(fid_header['nEchoes'])):
                for n in xrange(0, int(fid_header['nChannels'])):
                    for islice in xrange(0, int(dims[2])):
                        # ksp(:, :, islice, n, echo) = complex(ksp_data_real(:,
                        # echo:fid_header['nEchoes']:end,
                        # n:fid_header['nChannels']:end), ksp_data_imag(:,
                        # echo:fid_header['nEchoes']:end,
                        # n:fid_header['nChannels']:end))
                        ksp[:, :, islice, n, echo].real = ksp_data_real[:, echo::fid_header['nEchoes'],
                                                                        n::fid_header['nChannels']]
                        ksp[:, :, islice, n, echo].imag = ksp_data_imag[:, echo::fid_header['nEchoes'],
                                                                        n::fid_header['nChannels']]
                        img[:, :, islice, n, echo] = fftshift(ifftn(ifftshift(
                            ksp[:, pelist, islice, n, echo])))
        elif fid_header['nEchoes'] == 1:
            print "Reconstructing mode 3: 2D Multi channel."

            #ksp.real = np.reshape(ksp_data_real,[dims[0], dims[1], dims[2], fid_header['nChannels'], 1])
            # ksp.imag = np.reshape(ksp_data_imag,[dims[0], dims[1], dims[2], fid_header['nChannels'], 1])

            for n in xrange(0, int(fid_header['nChannels'])):
                if args.verbose:
                    print "Processing  channel ", n

                for islice in xrange(0, int(dims[2])):
                    ksp[:, :, islice, n, echo].real = ksp_data_real[
                        :, islice * dims[2] + 1:(islice + 1) * dims[2], n]  # islice::dims[2], n]
                    ksp[:, :, islice, n, echo].imag = ksp_data_imag[
                        :, islice * dims[2] + 1:(islice + 1) * dims[2], n]  # islice::dims[2], n]
                    if 'pelist' in locals():
                        img[:, :, islice, n, 0] = fftshift(
                            ifftn(ifftshift(ksp[:, pelist, islice, n, 0])))
                    else:
                        img[:, :, islice, n, 0] = fftshift(
                            ifftn(ifftshift(ksp[:, :, islice, n, 0])))
    else:  # if procpar.nD == 3
        if fid_header['nEchoes'] == 1 and fid_header['nChannels'] == 1:
            print 'Reconstructing mode: 3D basic'
            ksp = np.empty(
                [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32
            img = np.empty(
                [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32
            # [:, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
            ksp[:, :, :].real = ksp_data_real
            # [:, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
            ksp[:, :, :].imag = ksp_data_imag
            if 'pelist' in locals():
                img[:, :, :] = fftshift(ifftn(ifftshift(ksp[:, pelist, :])))
            else:
                img[:, :, :] = fftshift(ifftn(ifftshift(ksp[:, :, :])))
        elif fid_header['nChannels'] == 1:
            print 'Reconstructing mode 4: 3D Multi echo '
            # if len(ksp_data_real.shape) > 3:
            #     if fid_header['nEchoes'] != ksp_data_real.shape(3):
            #         ksp_data_real = np.reshape(ksp_data_real,[fid_header['dims'], fid_header['nChannels'], fid_header['nEchoes']])
            #         ksp_data_imag = np.reshape(ksp_data_imag,[fid_header['dims'], fid_header['nEchoes'], fid_header['nChannels']])
            #     elif fid_header['nChannels'] != ksp_data_real.shape(3):
            #         ksp_data_real = np.reshape(ksp_data_real,[fid_header['dims'], fid_header['nChannels'], fid_header['nEchoes']])
            #         ksp_data_imag = np.reshape(ksp_data_imag,[fid_header['dims'], fid_header['nChannels'], fid_header['nEchoes']])

            for echo in xrange(0, int(fid_header['nEchoes'])):
                for n in xrange(0, int(fid_header['nChannels'])):
                    if args.verbose:
                        print "Processing echo ", echo, " channel ", n
                    ksp[:, :, :, n, echo].real = ksp_data_real[:, echo::fid_header['nEchoes'],
                                                               n::fid_header['nChannels']]
                    ksp[:, :, :, n, echo].imag = ksp_data_imag[:, echo::fid_header['nEchoes'],
                                                               n::fid_header['nChannels']]
                    if 'pelist' in locals():
                        img[:, :, :, n, echo] = fftshift(
                            ifftn(ifftshift(ksp[:, pelist, :, n, echo])))
                    else:
                        img[:, :, :, n, echo] = fftshift(
                            ifftn(ifftshift(ksp[:, :, :, n, echo])))
        else:
            print "Reconstructing mode 4: 3D Multi channel."

            # ksp.real = np.reshape(ksp_data_real,
            #                       [dims[0], dims[1], dims[2],
            #                        fid_header['nChannels'], 1])
            # ksp.imag = np.reshape(ksp_data_imag,
            #                       [dims[0], dims[1], dims[2],
            #                        fid_header['nChannels'], 1])
            print 'KSP shape ', ksp.shape
            for n in xrange(0, int(fid_header['nChannels'])):
                for islice in xrange(0, int(dims[2])):
                    ksp[:, :, :, n, 0].real = ksp_data_real[
                        :, n::fid_header['nChannels'], n]
                    ksp[:, :, :, n, 0].imag = ksp_data_imag[
                        :, n::fid_header['nChannels'], n]

                if args.verbose:
                    print "Processing echo 0 channel ", n
                if 'pelist' in locals():
                    img[:, :, :, n, 0] = fftshift(
                        ifftn(ifftshift(ksp[:, pelist, :, n, 0])))
                else:
                    img[:, :, :, n, 0] = fftshift(
                        ifftn(ifftshift(ksp[:, :, :, n, 0])))
    return img, ksp
# end recon


def convksp(procpar, dims, fid_header, ksp_data_real, ksp_data_imag, args):
    """convksp Convert k-space image data into N-D k-space
    :param procpar:   procpar dictionary
    :param dims: dimension array
    :param fid_header: header info in fid
    :param ksp_data_real: real component of image k-space
    :param ksp_data_imag: imaginary component of image k-space
    """
    if args.verbose:
        print 'Reordering k-space'
        print dims[0], dims[1], dims[2], fid_header['nChannels'], fid_header['nEchoes']
        print 'nchannels ', fid_header['nblocks'], 'acqcycles ', procpar['acqcycles']
        print 'Image shape ', ksp_data_real.shape

    if np.product(dims) != np.product(ksp_data_real.shape):
        print "ksp not arranged properly"

    if fid_header['nChannels'] == 1 and fid_header['nEchoes'] == 1:
        ksp = np.empty(
            [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32
    else:
        ksp = np.empty([dims[0], dims[1], dims[2],
                        fid_header['nChannels'], fid_header['nEchoes']],
                       dtype=np.complex64)
        # float32
    if 'pelist' in procpar.keys():
        pelist = np.array(procpar['pelist']).astype(
            int) - int(min(procpar['pelist']))

    if procpar['nD'] == 2 and procpar['ni2'] == 1:
        if fid_header['nEchoes'] == 1 and fid_header['nChannels'] == 1:
            # two float32
            ksp = np.empty(
                [dims[0], dims[1], dims[2]], dtype=np.complex64)
            # [:, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
            ksp[:, :, :].real = ksp_data_real
            # [:, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
            ksp[:, :, :].imag = ksp_data_imag
            if 'pelist' in locals():
                ksp[:, :, :] = ksp[:, pelist, :]

        elif fid_header['nChannels'] == 1:
            for echo in xrange(0, int(fid_header['nEchoes'])):
                for n in xrange(0, int(fid_header['nChannels'])):
                    for islice in xrange(0, int(dims[2])):
                        ksp[:, :, islice, n, echo].real = ksp_data_real[
                            :, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
                        ksp[:, :, islice, n, echo].imag = ksp_data_imag[
                            :, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
                        if 'pelist' in locals():
                            ksp[:, :, islice, n, echo] = ksp[
                                :, pelist, islice, n, echo]
        elif fid_header['nEchoes'] == 1:
            print "Reconstructing mode 3: 2D Multi channel."
            for n in xrange(0, int(fid_header['nChannels'])):
                if args.verbose:
                    print "Processing  channel ", n
                for islice in xrange(0, int(dims[2])):
                    ksp[:, :, islice, n, echo].real = ksp_data_real[
                        :, islice * dims[2] + 1:(islice + 1) * dims[2], n]  # islice::dims[2], n]
                    ksp[:, :, islice, n, echo].imag = ksp_data_imag[
                        :, islice * dims[2] + 1:(islice + 1) * dims[2], n]  # islice::dims[2], n]
                    if 'pelist' in locals():
                        ksp[:, :, islice, n, 0] = ksp[:, pelist, islice, n, 0]

    else:  # if procpar.nD == 3
        if fid_header['nEchoes'] == 1 and fid_header['nChannels'] == 1:
            print 'Reorder kspace mode: 3D basic'
            ksp = np.empty(
                [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32
            # [:, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
            ksp[:, :, :].real = ksp_data_real
            # [:, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
            ksp[:, :, :].imag = ksp_data_imag
            if 'pelist' in locals():
                ksp[:, :, :] = ksp[:, pelist, :]
        elif fid_header['nChannels'] == 1:
            print 'Reconstructing mode 4: 3D Multi echo '
            print "Processing multi-echo or multi-channel 3D image"
            for echo in xrange(0, int(fid_header['nEchoes'])):
                for n in xrange(0, int(fid_header['nChannels'])):
                    if args.verbose:
                        print "Processing echo ", echo, " channel ", n
                    ksp[:, :, :, n, echo].real = ksp_data_real[:, echo::fid_header['nEchoes'],
                                                               n::fid_header['nChannels']]
                    ksp[:, :, :, n, echo].imag = ksp_data_imag[
                        :, echo::fid_header['nEchoes'], n::fid_header['nChannels']]
                    if 'pelist' in locals():
                        ksp[:, :, islice, n, echo] = ksp[
                            :, pelist, islice, n, echo]
        else:
            print "Reconstructing mode: 3D Multi channel."
            for n in xrange(0, int(fid_header['nChannels'])):
                for islice in xrange(0, int(dims[2])):
                    ksp[:, :, :, n, 0].real = ksp_data_real[
                        :, n::fid_header['nChannels'], n]
                    ksp[:, :, :, n, 0].imag = ksp_data_imag[
                        :, n::fid_header['nChannels'], n]
                    if 'pelist' in locals():
                        ksp[:, :, islice, n, echo] = ksp[
                            :, pelist, islice, n, echo]
    return ksp
# end convksp


def reconksp(procpar, dims, fid_header, ksp, args, pelistflag=0):
    """recon Reconstruct k-space image data into N-D image
    :param procpar:   procpar dictionary
    :param dims: dimension array
    :param fid_header: header info in fid
    :param ksp: k-space
    """
    if args.verbose:
        print 'Reconstructing image'
        print dims[0], dims[1], dims[2], fid_header['nChannels'], fid_header['nEchoes']
        print 'nchannels ', fid_header['nblocks'], 'acqcycles ', procpar['acqcycles']
        print 'Image shape ', ksp.shape

    if np.product(dims) != np.product(ksp.shape):
        print "ksp not arranged properly"

    if fid_header['nChannels'] == 1 and fid_header['nEchoes'] == 1:
        img = np.empty(
            [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32
    else:  # two float32
        img = np.empty([dims[0], dims[1], dims[2], fid_header['nChannels'],
                        fid_header['nEchoes']], dtype=np.complex64)  # float32
    # Setup pelist
    if 'pelist' in procpar.keys() and pelistflag == 0:
        pelist = np.array(procpar['pelist']).astype(
            int) - int(min(procpar['pelist']))
    if procpar['nD'] == 2 and procpar['ni2'] == 1:
        print 'Reconstructing mode 1: 2D slices'
        if fid_header['nEchoes'] == 1 and fid_header['nChannels'] == 1:
            for islice in xrange(0, int(dims[2])):
                # 'pelist'' in procpar.keys() and len(procpar['pelist']) == int(dims[2]):
                if 'pelist' in locals():
                    if args.verbose:
                        print procpar['pelist']
                    img[:, :, islice] = fftshift(ifftn(ifftshift(ksp[:,
                                                                     pelist,
                                                                     islice])))
                else:
                    img[:, :, islice] = fftshift(ifftn(ifftshift(ksp[:,
                                                                     :,
                                                                     islice])))
        else:  # if fid_header['nChannels'] == 1:
            print 'Reconstructing mode 2: 2D slices with ni2=' + procpar['ni2']
            print 'Echoes ', fid_header['nEchoes'], 'Channels ', fid_header['nChannels']
            for echo in xrange(0, int(fid_header['nEchoes'])):
                for n in xrange(0, int(fid_header['nChannels'])):
                    for islice in xrange(0, int(dims[2])):
                        if 'pelist' in locals():
                            img[:, :, islice, n, echo] = fftshift(ifftn(ifftshift(
                                ksp[:, pelist, islice, n, echo])))
                        else:
                            img[:, :, islice, n, echo] = fftshift(ifftn(ifftshift(
                                ksp[:, :, islice, n, echo])))
        # elif fid_header['nEchoes'] == 1:
        #     print "Reconstructing mode 3: 2D Multi channel."

        #     for n in xrange(0, int(fid_header['nChannels'])):
        #         if args.verbose:
        #             print "Processing  channel ", n
        #         for islice in xrange(0, int(dims[2])):

        #             if 'pelist' in locals():
        #                 img[:, :, islice, n, 0] = fftshift(
        #                     ifftn(ifftshift(ksp[:, pelist, islice, n, 0])))
        #             else:
        #                 img[:, :, islice, n, 0] = fftshift(
        #                     ifftn(ifftshift(ksp[:, :, islice, n, 0])))
    else:  # if procpar.nD == 3
        if fid_header['nEchoes'] == 1 and fid_header['nChannels'] == 1:
            print 'Reconstructing mode: 3D basic'

            img = np.empty(
                [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32

            if 'pelist' in locals():
                img[:, :, :] = fftshift(ifftn(ifftshift(ksp[:, pelist, :])))
            else:
                img[:, :, :] = fftshift(ifftn(ifftshift(ksp[:, :, :])))
        elif fid_header['nChannels'] == 1:
            print 'Reconstructing mode 4: 3D Multi echo '

            for echo in xrange(0, int(fid_header['nEchoes'])):
                for n in xrange(0, int(fid_header['nChannels'])):
                    if args.verbose:
                        print "Processing echo ", echo, " channel ", n
                    if 'pelist' in locals():
                        img[:, :, :, n, echo] = fftshift(
                            ifftn(ifftshift(ksp[:, pelist, :, n, echo])))
                    else:
                        img[:, :, :, n, echo] = fftshift(
                            ifftn(ifftshift(ksp[:, :, :, n, echo])))
        else:
            print "Reconstructing mode 4: 3D Multi channel."

            # ksp.real = np.reshape(ksp_data_real,
            #                       [dims[0], dims[1], dims[2],
            #                        fid_header['nChannels'], 1])
            # ksp.imag = np.reshape(ksp_data_imag,
            #                       [dims[0], dims[1], dims[2],
            #                        fid_header['nChannels'], 1])
            print 'KSP shape ', ksp.shape
            for n in xrange(0, int(fid_header['nChannels'])):
                if args.verbose:
                    print "Processing echo 0 channel ", n
                if 'pelist' in locals():
                    img[:, :, :, n, 0] = fftshift(
                        ifftn(ifftshift(ksp[:, pelist, :, n, 0])))
                else:
                    img[:, :, :, n, 0] = fftshift(
                        ifftn(ifftshift(ksp[:, :, :, n, 0])))
    return img
# end reconksp


def sliceifft(procpar,  ksp, args):
    """sliceifft Reconstruct k-space image data into 3 or 5-D image by slice in 3rd dim
    :param procpar:   procpar dictionary
    :param ksp: k-space
    :return img:
    """
    dims = ksp.shape
    if args.verbose:
        print 'Reconstructing image'
        print 'Image shape ', ksp.shape

    if len(dims) == 3:
        img = np.empty(
            [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32
        for islice in xrange(0, dims[2]):
            img[:, :, islice] = fftshift(
                ifftn(ifftshift(ksp[:, :, islice])))
    else:  # two float32
        img = np.empty(
            [dims[0], dims[1], dims[2], dims[3], dims[4]], dtype=np.complex64)  # float32
        if args.verbose:
            print 'Reconstructing mode 2: 2D slices with ni2=' + procpar['ni2']
            print 'Echoes ', dims[4], 'Channels ', dims[3]
        for echo in xrange(0, dims[4]):
            for n in xrange(0, dims[3]):
                for islice in xrange(0, dims[2]):
                    img[:, :, islice, n, echo] = fftshift(ifftn(ifftshift(
                        ksp[:, :, islice, n, echo])))
    return img
# end sliceifft


def simpleifft(procpar, dims, fid_header, ksp, args):
    """recon Reconstruct k-space image data into N-D image
    :param procpar:   procpar dictionary
    :param dims: dimension array
    :param fid_header: header info in fid
    :param ksp: k-space data
    :return img: reconstructed image
    """
    if args.verbose:
        print 'Reconstructing image'
        print dims[0], dims[1], dims[2], fid_header['nChannels'], fid_header['nEchoes']
        print 'nchannels ', fid_header['nblocks'], 'acqcycles ', procpar['acqcycles']
        print 'Image shape ', ksp.shape

    if np.product(dims) != np.product(ksp.shape):
        print "ksp not arranged properly ", dims, ksp.shape

    if len(ksp.shape) == 3 or (fid_header['nChannels'] == 1 and fid_header['nEchoes'] == 1):
        img = np.empty(
            [dims[0], dims[1], dims[2]], dtype=np.complex64)  # float32
    else:  # two float32
        img = np.empty([dims[0], dims[1], dims[2], fid_header['nChannels'],
                        fid_header['nEchoes']], dtype=np.complex64)  # float32

    if procpar['nD'] == 2 and procpar['ni2'] == 1:
        if args.verbose:
            print 'Reconstructing mode 1: 2D slices'
        if len(ksp.shape) == 3 or (fid_header['nEchoes'] == 1 and fid_header['nChannels'] == 1):
            for islice in xrange(0, int(dims[2])):
                img[:, :, islice] = fftshift(
                    ifftn(ifftshift(ksp[:, :, islice])))
        else:  # if fid_header['nChannels'] == 1:
            if args.verbose:
                print 'Reconstructing mode 2: 2D slices with ni2=' + procpar['ni2']
                print 'Echoes ', fid_header['nEchoes'], 'Channels ', fid_header['nChannels']
            for echo in xrange(0, int(fid_header['nEchoes'])):
                for n in xrange(0, int(fid_header['nChannels'])):
                    for islice in xrange(0, int(dims[2])):
                        img[:, :, islice, n, echo] = fftshift(ifftn(ifftshift(
                            ksp[:, :, islice, n, echo])))
    else:
        # if len(ksp.shape) == 3:
        if len(ksp.shape) == 3 or (fid_header['nEchoes'] == 1 and fid_header['nChannels'] == 1):
            if args.verbose:
                print 'Reconstructing mode: 3D basic'
            img = fftshift(ifftn(ifftshift(ksp)))
        else:
            if args.verbose:
                print 'Reconstructing mode 4: 3D Multi echo Multi channel '
            for echo in xrange(0, int(fid_header['nEchoes'])):
                for n in xrange(0, int(fid_header['nChannels'])):
                    img[:, :, :, n, echo] = fftshift(
                        ifftn(ifftshift(ksp[:, :, :, n, echo])))
    return img
# end simpleifft


def save_as_nifti(image, basename):
    """save_nifti
    Save image as NIFTI
    """
    import nibabel as nib
    affine = np.eye(4)
    if image.ndim == 5:
        for echo in xrange(0, image.shape[4]):
            for channel in xrange(0, image.shape[3]):
                new_image = nib.Nifti1Image(
                    np.abs(image[:, :, :, channel, echo]), affine)
                new_image.set_data_dtype(np.float32)
                nib.save(new_image, basename + '_' + str(channel) + str(echo)
                         + '.nii.gz')
    else:
        new_image = nib.Nifti1Image(np.abs(image), affine)
        new_image.set_data_dtype(np.float32)
        nib.save(new_image, basename + '.nii.gz')


def cuda_ifft(cmplx_array):
    from reikna.cluda import dtypes, any_api
    from reikna.fft import FFT
    from reikna.core import Annotation, Type, Transformation, Parameter

    # Pick the first available GPGPU API and make a Thread on it.
    api = any_api()
    thr = api.Thread.create()
    data_dev = thr.to_device(cmplx_array)
    # Create the FFT computation and attach the transformation above to its
    # input.
    # (A shortcut: using the array type saved in the transformation)
    ifft = FFT(data_dev)
    # fft.parameter.input.connect(trf, trf.output, new_input=trf.input)
    cifft = ifft.compile(thr)
    # Run the computation
    cifft(data_dev, data_dev, inverse=0)
    result = data_dev.get()
    return result


def double_resolution2(ksp, basename, procpar, hdr, args):
    """double_resolution creates double resolution image from k-space data
    based on super-resolution methodsfor multiple averages, this just expands the
    kspace data and reconstructs image.  Equivalent to interpolation in image space.
    """
    print "Double res and " + basename + " filter"
    # two 32-bit float

    szmin = np.array(ksp.shape[0:3]) / 2 - 1
    szmax = np.array(ksp.shape[0:3]) + szmin
    if len(ksp.shape) == 3:
        ksplarge = np.zeros(np.array(ksp.shape[:3]) * 2, dtype=np.complex64)
        ksplarge[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]] = ksp
        print "Double resolution k-space created. Starting reconstruction ...(may take some time)"

        image_filtered = simpleifft(
            procpar, ksplarge.shape, hdr, ksplarge, args)
        print 'Double res recon ' + basename
        del ksplarge
        print "Saving double resolution magnitude image: " + basename + " filtered"
        save_as_nifti(np.abs(image_filtered), basename + '-double')
        if args.phase:
            print "Saving double-resolution phase image: " + basename + " filtered"
            save_as_nifti(np.angle(image_filtered), basename + '-double-pha')
    else:
        ksplarge = np.zeros(np.array(ksp.shape[:3]) * 2, dtype=np.complex64)
        for echo in xrange(0, int(hdr['nEchoes'])):
            for n in xrange(0, int(hdr['nChannels'])):
                ksplarge[szmin[0]:szmax[0], szmin[1]:szmax[1], szmin[2]:szmax[2]] = ksp[
                    :, :, :, n, echo]
                print "Double resolution k-space created. Starting reconstruction ...(may take some time)"

                image_filtered = simpleifft(
                    procpar, ksplarge.shape, hdr, ksplarge, args)
                print 'Double res recon ' + basename
                print "Saving double resolution magnitude image: " + basename + " filtered"
                save_as_nifti(
                    np.abs(image_filtered), basename + '-double_0' + str(n) + '_0' + str(echo))
                if args.phase:
                    print "Saving double-resolution phase image: " + basename + " filtered"
                    save_as_nifti(
                        np.angle(image_filtered), basename + '-double-pha_0' + str(n) + '_0' + str(echo))
    return image_filtered


def RescaleFIDImage(ds, image_data, args):
    """ RescaleFIDImage
    """
    # Read in data from all files to determine scaling
    datamin = float("inf")
    datamax = float("-inf")

    # RescaleSlope of phase imgs set to [-pi, pi]
    if ds.ImageType[2] == "PHASE MAP" or \
            (hasattr(ds, 'ComplexImageComponent') and ds.ComplexImageComponent == 'PHASE'):
        # this implies either args.phase is on or procpar['imPH']=='y'
        datamin = -math.pi
        datamax = math.pi
    else:
        datamin = np.min([datamin, image_data.min()])
        datamax = np.max([datamax, image_data.max()])

    RescaleIntercept = datamin
    # Np recast to int16 with range (-32768 or 32767)
    if ds.PixelRepresentation == 0:
        slope_factor = RescaleFDF.UInt16MaxRange
    else:
        slope_factor = RescaleFDF.Int16MaxRange
    RescaleSlope = (datamax - datamin) / slope_factor

    if args.verbose:
        print "Rescale data to uint16"
        print "Intercept: ", RescaleIntercept, "  Slope: ", RescaleSlope
        print "Current data min: ", image_data.min(), " max ", image_data.max()
    image_data = (image_data - RescaleIntercept) / RescaleSlope
    image_data_int16 = image_data.astype(np.int16)

    # Adjusting Dicom parameters for rescaling
    ds.RescaleIntercept = RescaleFDF.ShortenFloatString(
        RescaleIntercept, "RescaleIntercept")  # (0028,1052) Rescale Intercept
    ds.RescaleSlope = RescaleFDF.ShortenFloatString(
        RescaleSlope, "RescaleSlope")  # (0028,1053) Rescale Slope

    return ds, image_data_int16
# end RescaleFIDImage


def flipdims(image):
    """flipdims flip dimensions on each axis 3D image
    """
    from scipy.signal import _arraytools as ar
    for n in xrange(0, 3):
        print "Reversing axis ", n
        image = ar.axis_reverse(image, n)
    return image


def RearrangeImage(image, axis_order, args):
    """Rearrange Image
    """
    from scipy.signal import _arraytools as ar
    axes = re.findall(r'[-]*\d', axis_order)
    raxes = np.array(axes)
    iaxes = np.abs(raxes.astype(int))
    dims = image.shape
    if len(axes) != 3:
        print " Axis order not compatible.  Must be arragment of 0, 1, 2 with no spaces and a comma in between."
    else:

        if image.ndim == 5:
            image_treal = np.empty([dims[iaxes[0]], dims[iaxes[1]],
                                    dims[iaxes[2]], dims[3],
                                    dims[4]], dtype=np.float32)
            image_timag = np.empty([dims[iaxes[0]], dims[iaxes[1]],
                                    dims[iaxes[2]], dims[3],
                                    dims[4]], dtype=np.float32)
            if args.verbose:
                print "Transposing 5D image to ", axis_order
                print image.shape, ' -> ', image_treal.shape

            for echo in xrange(0, image.shape[4]):
                # Transpose image dimensions
                image_treal[:, :, :, 0, echo] = np.transpose(
                    (image.real)[:, :, :, 0, echo], (iaxes[0], iaxes[1], iaxes[2]))
                image_timag[:, :, :, 0, echo] = np.transpose(
                    (image.imag)[:, :, :, 0, echo], (iaxes[0], iaxes[1], iaxes[2]))
                # Then do reversing second
                for n in xrange(0, 3):
                    if re.search('-', axes[n]):
                        if args.verbose:
                            print "Reversing axis ", n
                        image_treal[:, :, :, 0, echo] = ar.axis_reverse(
                            image_treal[:, :, :, 0, echo], n)
                        image_timag[:, :, :, 0, echo] = ar.axis_reverse(
                            image_timag[:, :, :, 0, echo], n)
                image.reshape([dims[iaxes[0]], dims[iaxes[1]],
                               dims[iaxes[2]], dims[3], dims[4]])
#                print image.shape
        else:
            if args.verbose:
                print "Transposing 3D image to ", axis_order
            for n in xrange(0, 3):
                if re.search('-', axes[n]):
                    if args.verbose:
                        print "Reversing axis ", n
                    image_treal = ar.axis_reverse(image_treal, n)
                    image_timag = ar.axis_reverse(image_timag, n)
            image_treal = np.transpose(
                image.real, (iaxes[0], iaxes[1], iaxes[2]))
            image_timag = np.transpose(
                image.imag, (iaxes[0], iaxes[1], iaxes[2]))
            image.reshape([dims[iaxes[0]], dims[iaxes[1]], dims[iaxes[2]]])
            # print image.shape
        image = np.empty([dims[iaxes[0]], dims[iaxes[1]],
                          dims[iaxes[2]], dims[3], dims[4]],
                         dtype=np.complex64)
        image.real = image_treal
        image.imag = image_timag
        if args.verbose:
            print "Final rearranged image shape: ", image.shape, image_treal.shape
    return image
# end RearrangeImage


def AssertImplementation(testval, fidfilename, comment, assumption):
    """ASSERTksp_data_imagPLEMENTATION - Check FID properties match up with interpretation of procpar
  Due to lack of documentation on the use of the procpar file, the interpretation
 implemented in this script is based on various documents and scripts and may have errors.
 This function seeks to double check some of the interpretations against the fid properties.
    """
    if testval:
        if len(fidfilename) > 0:
            FIDStr = "fid file: " + fidfilename + "\n"
        else:
            FIDStr = ""

        print "\nImplementation check error:\n" + FIDStr + comment + '\nAssumption:' + assumption + '\n'
        # sys.exit(1)


def ParseDiffusionFID(ds, procpar, diffusion_idx, args):
    """ParseDiffusionFID is a variation of ParseDiffusionFDF

    :param ds: Dicom dataset
    :param procpar: Procpar dictionary tag/value pairs
    :param fid_properties: Tag/value pairs of local fid file

    :param args: Input arguments
    :returns: Dicom struct

    VNMRJ 4.0:

SGLDiffusion Module
    The DIFFUSION_T structure is defined to have the following structure members in
sglCommon.h
typedef struct
{
      DIFF_TYPE_T     type;                  /* type of diffusion encoding (GE, SE or STE)*/
      GENERIC_GRADIENT_T    *grad;           /* diffusion gradient */
      double gdiff;          /* diffusion gradient amplitude */
      double delta;         /* diffusion delta (gradient duration) */
      double DELTA;         /* diffusion DELTA */
      double tadd;          /* additional diffusion encoding time */
      double te;            /* echo time for diffusion */
      char     minte;       /* minimum echo time flag for diffusion */
      double tau1;          /* duration of events in first half of te for SE diffusion */
      double tau2;          /* duration of events in second half of te for SE diffusion */
      double tm;            /* mixing time for STE diffusion */
      char     mintm;       /* minimum mixing time flag for STE diffusion */
      double d1;            /* d1,d2,d3,d4 are delays around the diffusion gradients. */
      double d2;            /* The scheme is [d1 - Gdiff - d2] - taudiff - [d3 - Gdiff - d4] */
      double d3;
      double d4;
      double dm;            /* diffusion mixing time delay */
      double *dro;          /* multiplier for diffusion gradient in readout */
      double *dpe;          /* multiplier for diffusion gradient in phase encode */
      double *dsl;          /* multiplier for diffusion gradient in slice */
      int      nbval;       /* total number of bvalues*directions */
      int      nbro;        /* number of bvalues*directions along readout */
      int      nbpe;        /* number of bvalues*directions along phase encode */
      int      nbsl;        /* number of bvalues*directions along slice */
      double *bro;          /* b value along readout */
      double *bpe;          /* b value along phase encode */
      double *bsl;          /* b value along slice */
      double *brp;          /* b value cross-term (readout - phase encode) */
      double *brs;          /* b value cross-term (readout - slice) */
      double *bsp;          /* b value cross-term (slice - phase encode) */
      double *btrace;       /* trace */
      double max_bval;      /* the maximum trace */
      double Time;          /* the duration of diffusion components */
} DIFFUSION_T;
typedef enum {
      DIFF_NULL = 0,
      DIFF_GE = 1,
      DIFF_SE = 2,
      DIFF_STE = 3,
} DIFF_TYPE_T;


    """
    if args.verbose:
        print 'Processing diffusion image'

    # Get procpar diffusion parameters
    Bvalue = procpar['bvalue']  # 64 element array for 31 directions
    # 10 for 3 directions
    BValueSortIdx = np.argsort(Bvalue)
    BvalSave = procpar['bvalSave']
    if 'bvalvs' in procpar.keys():
        BvalVS = procpar['bvalvs']  # excluded in external recons by vnmrj
    BvalueRS = procpar['bvalrs']  # 64
    BvalueRR = procpar['bvalrr']  # 64
    BvalueRP = procpar['bvalrp']  # 64
    BvaluePP = procpar['bvalpp']  # 64
    BvalueSP = procpar['bvalsp']  # 64
    BvalueSS = procpar['bvalss']  # 64

    # Assume all images were reconned by this program
    # if procpar['recon'] == 'external':
    # diffusion_idx=0
    # while True:
    #    if (math.fabs(Bvalue[ diffusion_idx ] - fdf_properties['bvalue']) < 0.005):
    #        break
    #    diffusion_idx += 1
    #diffusion_idx = fdf_properties['array_index'] - 1
    # else:
    #    diffusion_idx = fdf_properties['array_index']*2

    # if diffusion_idx > len(Bvalue):
    # print 'Procpar Bvalue does not contain enough values determined by
    # fdf_properties array_index'

    # if args.verbose:
    # print 'Diffusion index ', diffusion_idx, ' arrary index ',
    # fdf_properties['array_index']

    # Sort diffusion based on sorted index of Bvalue instead of
    if ds.MRAcquisitionType == '2D':
        fdf_properties['array_index']
        ds.FrameAcquisitionNumber = fdf_properties['array_index']
        loc = numpy.array(fdf_properties['location'], dtype='|S9')
        ds.ImagePositionPatient = [loc[0], loc[1],loc[2]]
        orient = numpy.array(fdf_properties['orientation'], dtype='|S9')
        ds.ImageOrientationPatient = [orient[0], orient[1], orient[2],
                                      orient[3], orient[4], orient[5],
                                      orient[6], orient[7], orient[8]]    
    else:
        ds.AcquisitionNumber = BValueSortIdx[diffusion_idx]

#    if (math.fabs(Bvalue[ diffusion_idx ] - fid_properties['bvalue']) > 0.005):
# print 'Procpar and fdf B-value mismatch: procpar value ', Bvalue[
# diffusion_idx ], ' and  local fdf value ', fdf_properties['bvalue'], '
# array idx ', fdf_properties['array_index']

    # MR Diffusion Sequence (0018,9117) see DiffusionMacro.txt
    # B0 scan does not need the MR Diffusion Gradient Direction Sequence macro and its directionality should be set to NONE
    # the remaining scans relate to particular directions hence need the
    # direction macro

    diffusionseq = Dataset()
    if Bvalue[diffusion_idx] < 20:
        if Bvalue[diffusion_idx] != 0:
            print "Diffusion reference not quite zero, setting to zero instead!"
        diffusionseq.DiffusionBValue = 0
        diffusionseq.DiffusionDirectionality = 'NONE'
    else:
        diffusionseq.DiffusionBValue = int(Bvalue[diffusion_idx])
        # TODO  One of: DIRECTIONAL, BMATRIX, ISOTROPIC, NONE
        diffusionseq.DiffusionDirectionality = 'BMATRIX'

    # Diffusion Gradient Direction Sequence (0018,9076)
        diffusiongraddirseq = Dataset()
        # Diffusion Gradient Orientation (0018,9089)
        # diffusiongraddirseq.add_new((0x0018,0x9089), 'FD',[
        #fdf_properties['dro'], fdf_properties['dpe'],
        # fdf_properties['dsl']])
        diffusiongraddirseq.DiffusionGradientOrientation = [
            procpar['dro'][diffusion_idx], procpar['dpe'][diffusion_idx],
            procpar['dsl'][0]]
        diffusionseq.DiffusionGradientDirectionSequence = Sequence(
            [diffusiongraddirseq])
        #diffusionseq.add_new((0x0018,0x9076), 'SQ',Sequence([diffusiongraddirseq]))

    # Diffusion b-matrix Sequence (0018,9601)
        diffbmatseq = Dataset()
        diffbmatseq.DiffusionBValueXX = BvalueRR[diffusion_idx]
        diffbmatseq.DiffusionBValueXY = BvalueRP[diffusion_idx]
        diffbmatseq.DiffusionBValueXZ = BvalueRS[diffusion_idx]
        diffbmatseq.DiffusionBValueYY = BvaluePP[diffusion_idx]
        diffbmatseq.DiffusionBValueYZ = BvalueSP[diffusion_idx]
        diffbmatseq.DiffusionBValueZZ = BvalueSS[diffusion_idx]
        diffusionseq.DiffusionBMatrixSequence = Sequence([diffbmatseq])

    # TODO  One of: FRACTIONAL, RELATIVE, VOLUME_RATIO
    diffusionseq.DiffusionAnisotropyType = 'FRACTIONAL'
    ds.MRDiffusionSequence = Sequence([diffusionseq])

    MRImageFrameType = Dataset()
    MRImageFrameType.FrameType = [
        "ORIGINAL", "PRksp_data_imagARY", "DIFFUSION", "NONE"]  # same as ds.ImageType
    MRImageFrameType.PixelPresentation = ["MONOCHROME"]
    MRImageFrameType.VolumetrixProperties = ["VOLUME"]
    MRImageFrameType.VolumeBasedCalculationTechnique = ["NONE"]
    MRImageFrameType.ComplexImageComponent = ["MAGNITUDE"]
    MRImageFrameType.AcquisitionContrast = ["DIFFUSION"]
    ds.MRImageFrameTypeSequence = Sequence([MRImageFrameType])

    return ds


def ParseASL(ds, procpar, fid_properties):
    """ParseASL
    FIXME- this is not implemented for FID
    """
    # (0018,9257)	    1C	The purpose of the Arterial Spin Labeling.
    #   Enumerated Values:
    #	       LABEL
    #              CONTROL
    #		M_ZERO_SCAN
    #	Required if Frame Type (0008,9007) is
    #	ORIGINAL. May be present otherwise.
    #	See C.8.13.5.14.1 for further
    #	explanation.
    # if fdf_properties["asltag"] == 1:
    #     ds.MRArterialSpinLabeling[0].ASLContext = 'LABEL'
    # elif fdf_properties["asltag"] == -1:
    #     ds.MRArterialSpinLabeling[0].ASLContext = 'CONTROL'
    # else:
    #     ds.MRArterialSpinLabeling[0].ASLContext = 'M_ZERO_SCAN'

    # # FIX ME : this could be either array_index or slice_no
    # ds.MRArterialSpinLabeling[0].ASLSlabSequence[
    #     0].ASLSlabNumber = fdf_properties["array_index"]

    # # ASL Mid slab position

    # # The Image Plane Attributes, in conjunction with the Pixel Spacing
    # # Attribute, describe the position and orientation of the image slices
    # # relative to the patient-based coordinate system. In each image frame the
    # # Image Position (Patient) (0020,0032) specifies the origin of the image
    # # with respect to the patient-based coordinate system. RCS and the Image
    # # Orientation (Patient) (0020,0037) attribute values specify the
    # # orientation of the image frame rows and columns. The mapping of pixel
    # # location i, j to the RCS is calculated as follows:
    # #					 X x i Yx j 0 S x
    # #				Px				  i	   i
    # #					 X y i Yy j 0 S y
    # #				Py				  j	   j
    # #								    =M
    # #					 X z i Yz j 0 S z	  0	   0
    # #				Pz
    # #				1	   0	   0	  01	  1	   1
    # # Where:
    # #	 Pxyz The coordinates of the voxel (i,j) in the frame's
    # # image plane in units of mm.
    # #	 Sxyz The three values of the Image Position (Patient) (0020,0032)
    # # attributes. It is the
    # #	       location in mm from the origin of the RCS.
    # #	 Xxyz The values from the row (X) direction cosine of the Image
    # # Orientation (Patient)
    # #	       (0020,0037) attribute.
    # #	 Yxyz The values from the column (Y) direction cosine of the Image
    # # Orientation (Patient)
    # #	       (0020,0037) attribute.
    # #	 i     Column index to the image plane. The first column is index zero.
    # #	   i Column pixel resolution of the Pixel Spacing (0028,0030)
    # # attribute in units of mm.
    # #	 j     Row index to the image plane. The first row index is zero.
    # # j Row pixel resolution of the Pixel Spacing (0028,0030) attribute in
    # # units of mm.

    # # ds.MRArterialSpinLabelingSequence.ASLSlabSequence[0].ASLMidSlabPosition
    # # = [str(ImagePositionPatient[0]), str(ImagePositionPatient[1]),
    # # str(ImagePositionPatient[2] + (islice-1)*SliceThickness))]

    # #            print ImagePositionPatient
    # #           M = np.matrix([[PixelSpacing[0] *
    # #           ImageOrientationPatient[0], PixelSpacing[1] *
    # #           ImageOrientationPatient[1], SliceThinkness *
    # #           ImageOrientationPatient[2] ImagePositionPatient[0]],
    # #           [PixelSpacing[0] * ImageOrientationPatient[3], PixelSpacing[1]
    # #           * ImageOrientationPatient[4], SliceThinkness *
    # #           ImageOrientationPatient[5] ImagePositionPatient[1]],
    # #           [PixelSpacing[0] * ImageOrientationPatient[6], PixelSpacing[1]
    # #           * ImageOrientationPatient[8], SliceThinkness *
    # #           ImageOrientationPatient[8] ImagePositionPatient[2]], [0, 0, 0,
    # #           1]])
    # #           pos = np.matrix([[ceil(ds.Rows / 2)],[ ceil(ds.Columns /
    # #           2],[fdf_properties['slice_no'],[1]])
    # #           Pxyz = M * pos

    # # ds.MRArterialSpinLabelingSequence.ASLSlabSequence[0].ASLMidSlabPosition
    # # = [str(Pxyz[0,0]),str(Pxyz[1,0]),str(Pxyz[2,0]))]
    return ds


def ParseFID(ds, fid_properties, procpar, args):
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
    #     fid_tmp = fid_properties['roi']
    #     fid_properties['roi'][0:1] = fid_tmp[1:2]
    #     fid_properties['roi'][2] = fid_tmp[0]
    #     fid_tmp = fid_properties['dims']
    #     fid_properties['dims'][0:1] = fid_tmp[1:2]
    #     fid_properties['dims'][2] = fid_tmp[0]

    #----------------------------------------------------------
    # General implementation checks
    filename = os.path.join(args.inputdir, 'fid')  # fid_properties['filename']

    # File dimensionality or Rank fields
    # rank is a positive integer value (1, 2, 3, 4,...) giving the
    # number of dimensions in the data file (e.g., int rank=2;).
    fidrank = fid_properties['rank']
    acqndims = procpar['acqdim']
    CommentStr = '''Acquisition dimensionality (ie 2D or 3D) does not match between fid and procpar'''
    AssumptionStr = '''Procpar nv2 > 0 indicates 3D acquisition and fid rank property indicates dimensionality.\n
        Using local FID value %d instead of procpar value %d.''' % (fidrank, acqndims)
    if args.verbose:
        print 'Acqdim (type): ' + ds.MRAcquisitionType + " acqndims " + str(acqndims)

    AssertImplementation(
        acqndims != fidrank, filename, CommentStr, AssumptionStr)

    # matrix is a set of rank integers giving the number of data
    # points in each dimension (e.g., for rank=2, float
    # matrix[]={256,256};)
    if fidrank == 3:
        fid_size_matrix = fid_properties['dims'][0:3]
    else:
        fid_size_matrix = fid_properties['dims'][0:2]
    if args.verbose:
        print "FID size matrix ", fid_size_matrix, type(fid_size_matrix)
    # fid_size_matrix = np.array(fid_matrix)

    # spatial_rank is a string ("none", "voxel", "1dfov", "2dfov",
    # "3dfov") for the type of data (e.g., char
    # *spatial_rank="2dfov";).
    # spatial_rank = fid_properties['spatial_rank']

    #  0018,0023 MR Acquisition Type (optional)
    # Identification of spatial data encoding scheme.
    # Defined Terms: 1D 2D 3D
    fid_MRAcquisitionType = '2D'
    if fidrank == 3:
        fid_MRAcquisitionType = '3D'
    CommentStr = 'MR Acquisition type does not match between fid and procpar'
    AssumptionStr = '''In fid, MR Acquisition type defined by spatial_rank and matrix.
        For 2D, spatial_rank="2dfov" and matrix has two elements eg. {256,256}.\n
        For 3D, spatial_rank="3dfov" and matrix has three elements.\n
        In procpar, MR Acquisition type is defined by nv2 > 0 or lpe2 > 0.\n
        Using local FID value %s
        instead of procpar value %s. ''' % (fid_MRAcquisitionType, fid_MRAcquisitionType)
    # AssertImplementation( ds.MRAcquisitionType != fid_MRAcquisitionType, filename, CommentStr, AssumptionStr)
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
    if fidrank == 3:
        roi = fid_properties['FOVcm'][0:3]
    else:
        roi = fid_properties['FOVcm'][0:2]
    if args.verbose:
        print "FID roi ", roi, type(roi)
    # roi = np.array(roi_text)

    # PixelSpacing - 0028,0030 Pixel Spacing (mandatory)
    PixelSpacing = fid_properties['voxelmm']
    if PixelSpacing[0] != float(ds.PixelSpacing[0]) or PixelSpacing[1] != float(ds.PixelSpacing[1]):
        print "Pixel spacing mismatch, procpar ", ds.PixelSpacing, " fid spacing ", str(PixelSpacing[0]), ', ', str(PixelSpacing[1])
    if args.verbose:
        print "Pixel Spacing : Procpar   ", ds.PixelSpacing
        print "Pixel Spacing : FID props ", PixelSpacing
    # (0028,0030) Pixel Spacing
    ds.PixelSpacing = [str(PixelSpacing[0]), str(PixelSpacing[1])]

    # FID slice thickness
    if fidrank == 3:
        fidthk = fid_properties['voxelmm'][2]
    else:
        fidthk = fid_properties['voxelmm'][2]
    #    fidthk = procpar['thk']*10

    CommentStr = 'Slice thickness does not match between fid and procpar'
    AssumptionStr = '''In fid, slice thickness defined by roi[2] for 2D or roi[2]/matrix[2].\n
        In procpar, slice thickness defined by thk (2D) or lpe2*10/(fn2/2) or lpe2*10/nv2.\n
        Using local FID value %s
        instead of procpar value %s .''' % (str(fidthk), str(ds.SliceThickness))
    if args.verbose:
        print 'fidthk : ' + str(fidthk)
        print 'SliceThinkness: ' + str(ds.SliceThickness)

    SliceThickness = float(ds.SliceThickness)

    # fix me Quick hack to avoid assert errors for diffusion and 3D magnitude images
    # if not ('diff' in procpar.keys() and procpar["diff"] == 'y'):
    #	 if MRAcquisitionType == '3D':
    #	     print 'Not testing slicethickness in diffusion and 3D MR FIDs'
    #	else:
    #AssertImplementation(SliceThickness != fidthk, filename, CommentStr, AssumptionStr)

    # Slice Thickness 0018,0050 Slice Thickness (optional)
    if fidrank == 3:
        if len(PixelSpacing) != 3:
            print "Slice thickness: 3D procpar spacing not available, fidthk ", fidthk
    else:
        if PixelSpacing[2] != ds.SliceThickness:
            print "Slice Thickness mismatch, procpar ", ds.SliceThickness, " fid spacing ", PixelSpacing[2], fidthk

    # Force slice thickness to be from fid props
    ds.SliceThickness = str(fidthk)
    SliceThickness = fidthk

    # -------------------------------------------------------------------------
    # GROUP 0020: Relationship

    ds.ImageComments = A2D.FID2DCM_Image_Comments + '\n' +\
        ', '.join("%s=%r" % (key, val)
                  for (key, val) in fid_properties.iteritems())
    # str(fid_properties) #['filetext']
    if args.verbose:
        print 'Calculating orientation and span'
        print fid_properties['orientation']
        print np.array(fid_properties['location']) * 10
        print fid_properties['span']
    orientation = np.matrix(fid_properties['orientation']).reshape(3, 3)
    location = np.matrix(np.array(fid_properties['location']) * 10.0)
    if len(fid_properties['span']) == 2:
        span = np.matrix(np.array(np.append(fid_properties['span'], 0) * 10.0))
    else:
        span = np.matrix(np.array(fid_properties['span']) * 10.0)
    if args.verbose:
        print "Span: ", span, span.shape
        print "Location: ", location, location.shape
        print orientation.transpose()

    ds, tmatrix = ProcparToDicomMap.CalcTransMatrix(ds, orientation,
                                                    location, span, fidrank,
                                                    PixelSpacing,
                                                    SliceThickness)
    if args.verbose:
        print 'Transformation matrix'
        print tmatrix
    # Nuclear Data Fields

    # Data fields may contain data generated by interactions between
    # more than one nucleus (e.g., a 2D chemical shift correlation map
    # between protons and carbon). Such data requires interpreting the
    # term ppm for the specific nucleus, if ppm to frequency
    # conversions are necessary, and properly labeling axes arising
    # from different nuclei. To properly interpret ppm and label axes,
    # the identity of the nucleus in question and the corresponding
    # nuclear resonance frequency are needed. These fields are related
    # to the abscissa values "ppm1", "ppm2", and "ppm3" in that the 1,
    # 2, and 3 are indices into the nucleus and nucfreq fields. That
    # is, the nucleus for the axis with abscissa string "ppm1" is the
    # first entry in the nucleus field.
    #   - nucleus is one entry ("H1", "F19", same as VNMR tn
    #       parameter) for each rf channel (e.g., char
    #       *nucleus[]={"H1","H1"};).
    #    - nucfreq is the nuclear
    #       frequency (floating point) used for each rf channel (e.g.,
    #       float nucfreq[]={200.067,200.067};).

#    if fid_properties['nucleus'][0] != ds.ImagedNucleus:
#        print 'Imaged nucleus mismatch: ', fid_properties['nucleus'], ds.ImagedNucleus
#    if math.fabs(fid_properties['nucfreq'][0] - float(ds.ImagingFrequency)) > 0.01:
# print 'Imaging frequency mismatch: ', fid_properties['nucfreq'],
# ds.ImagingFrequency

    # Change patient position and orientation in
    # if procpar['recon'] == 'external' and fid_properties['rank'] == '3':

    #-------------------------------------------------------------------------
    # GROUP 0028: Image Presentation
    # A good short description of this section can be found here:
    # http://dicomiseasy.blogspot.com.au/2012/08/chapter-12-pixel-data.html

    # Implementation check
    CommentStr = 'Number of rows does not match between fid and procpar'
    AssumptionStr = '''In FID, number of rows is defined by matrix[1]. \n
        In procpar, for 3D datasets number of rows is either fn1/2 or nv (%s,%s).\n
        For 2D datasets, number of rows is fn/2.0 or np (%s,%s).\n
        Using local FID value %s
        instead of procpar value %s.''' % (str(procpar['fn1'] / 2.0),
                                           str(procpar['nv']),
                                           str(procpar['fn'] / 2.0),
                                           str(procpar['np']),
                                           str(fid_properties['dims'][1]),
                                           str(ds.Rows))
    AssertImplementation(int(float(ds.Rows)) != int(float(fid_properties['dims'][1])),
                         filename, CommentStr, AssumptionStr)
    if args.verbose:
       # print 'Rows ', procpar['fn']/2.0, procpar['fn1']/2.0, procpar['nv'],
       # procpar['np']/2.0
        print '   Procpar: rows ', ds.Rows, '   FID prop rows ', fid_properties['dims'][1]
    ds.Rows = fid_properties['dims'][1]  # (0028,0010) Rows

    # Implementation check
    CommentStr = 'Number of columns does not match between fid and procpar'
    AssumptionStr = '''In FID, number of columns is defined by matrix[0]. \n
        In procpar, for 3D datasets number of columns is either fn/2 or np (%s,%s).\n
        For 2D datasets, number of rows is fn1/2.0 or nv (%s,%s).\n
        Using local FID value %s
         instead of procpar value %s.''' % (str(procpar['fn'] / 2.0),
                                            str(procpar['np']),
                                            str(procpar['fn1'] / 2.0),
                                            str(procpar['nv']),
                                            str(fid_properties['dims'][0]),
                                            str(ds.Columns))
    AssertImplementation(int(float(ds.Columns)) != int(float(fid_properties['dims'][0])),
                         filename, CommentStr, AssumptionStr)
    if args.verbose:
       # print 'Columns ', procpar['fn']/2.0, procpar['fn1']/2.0,
       # procpar['nv'], procpar['np']/2.0, fid_properties['dims'][0]
        print '   Procpar: Cols ', ds.Rows, '   FID prop Cols ', fid_properties['dims'][0]
    ds.Columns = fid_properties['dims'][0]  # (0028,0011) Columns

    #-------------------------------------------------------------------------
    # Number of frames
    # DICOMFID_HEADER code:
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
    #	    if $DEBUG then write('alpha',' new value = "%s"',$value) endif

    # if fidrank == 3:
    #	 ds.NumberOfFrames = fid_properties['dims'][2]

    # dicom3tool uses frames to create enhanced MR
    # ds.NumberOfFrames = fid_properties['slices']
    # ds.FrameAcquisitionNumber = fid_properties['slice_no']

#    if 'ne' in procpar.keys() and procpar['ne'] > 1:
#	 print 'Processing multi-echo sequence image'


#    if 'echo_no' in fid_properties.keys():
#        volume = fid_properties['echo_no']

    if ds.ImageType[2] == "MULTIECHO":  # and fid_properties['echoes'] > 1:
        if args.verbose:
            print 'Multi-echo sequence ...'
        # TE 0018,0081 Echo Time (in ms) (optional)
        if 'TE' in fid_properties.keys():
            if ds.EchoTime != str(fid_properties['TE'] * 1000):
                print "Echo Time mismatch: ", ds.EchoTime, fid_properties['TE'] * 1000
            ds.EchoTime = str(fid_properties['TE'] * 1000)
        # 0018,0086 Echo Number (optional)
    if 'echo_no' in fid_properties.keys():
        if ds.EchoNumber != fid_properties['echo_no']:
            print "Echo Number mismatch: ", ds.EchoNumber, fid_properties['echo_no']
        ds.EchoNumber = fid_properties['echo_no']

    if ds.ImageType[2] == "ASL":
        ds = ParseASL(ds, procpar, fid_properties)

    # if 'echoes' in fid_properties.keys() and fid_properties['echoes'] > 1 and fid_properties['array_dim'] == 1:
    #    ds.AcquisitionNumber = fid_properties['echo_no']
    #    ds.ImagesInAcquisition = fid_properties['echoes']
    # else:

    ds.AcquisitionNumber = '1'  # fid_properties['array_index']
    ds.ImagesInAcquisition = '1'  # fid_properties['array_dim']

    if ds.ImageType[2] == 'DIFFUSION':
        print 'Parsing diffusion module'
        ds = ParseDiffusionFID(ds, procpar, fid_properties, args)

    # Multi dimension Organisation and Index module
    DimOrgSeq = Dataset()
    #ds.add_new((0x0020,0x9164), 'UI', DimensionOrganizationUID)

    if (ds.ImageType[2] == "MULTIECHO") or (ds.ImageType[2] == "DIFFUSION" and ds.AcquisitionNumber == 1):  # or SEQUENCE == "Diffusion":
        DimensionOrganizationUID = [ProcparToDicomMap.CreateUID(A2D.UID_Type_DimensionIndex1,
                                                                [], [], args.verbose),
                                    ProcparToDicomMap.CreateUID(A2D.UID_Type_DimensionIndex2,
                                                                [], [], args.verbose)]
        DimOrgSeq.add_new((0x0020, 0x9164), 'UI', DimensionOrganizationUID)
        ds.DimensionOrganizationType = '3D_TEMPORAL'  # or 3D_TEMPORAL
    else:
        DimensionOrganizationUID = ProcparToDicomMap.CreateUID(
            A2D.UID_Type_DimensionIndex1, [], [], args.verbose)
        # if args.verbose:
        #    print "DimUID", DimensionOrganizationUID
        DimOrgSeq.add_new((0x0020, 0x9164), 'UI', [DimensionOrganizationUID])
        ds.DimensionOrganizationType = '3D'  # or 3D_TEMPORAL

    ds.DimensionOrganizationSequence = Sequence([DimOrgSeq])

    if ds.ImageType[2] == 'MULTIECHO':
        DimIndexSeq1 = Dataset()
        # Image position patient 20,32 or 20,12
        DimIndexSeq1.DimensionIndexPointer = (0x0020, 0x0032)

        # #DimIndexSeq1.DimensionIndexPrivateCreator=
        # #DimIndexSeq1.FunctionalGrouppointer=
        # #DimIndexSeq1.FunctionalGroupprivateCreator=
        DimIndexSeq1.add_new(
            (0x0020, 0x9164), 'UI', DimOrgSeq.DimensionOrganizationUID[0])
        DimIndexSeq1.DimensionDescriptionLabel = 'Third Spatial dimension'

        DimIndexSeq2 = Dataset()
        DimIndexSeq2.DimensionIndexPointer = (0x0018, 0x0081)  # Echo Time
        # DimIndexSeq2.DimensionIndexPrivateCreator=
        # DimIndexSeq2.FunctionalGrouppointer=
        # DimIndexSeq2.FunctionalGroupprivateCreator=
        DimIndexSeq2.add_new(
            (0x0020, 0x9164), 'UI', DimOrgSeq.DimensionOrganizationUID[1])
        DimIndexSeq2.DimensionDescriptionLabel = 'Fourth dimension (multiecho)'
        ds.DimensionIndexSequence = Sequence([DimIndexSeq2, DimIndexSeq1])
           
    elif (ds.ImageType[2] == "DIFFUSION" and ds.AcquisitionNumber == 1):
        DimIndexSeq1 = Dataset()
        # Image position patient 20,32 or 20,12
        DimIndexSeq1.DimensionIndexPointer = (0x0020, 0x0032)

        # #DimIndexSeq1.DimensionIndexPrivateCreator=
        # #DimIndexSeq1.FunctionalGroupPointer=
        # #DimIndexSeq1.FunctionalGroupPrivateCreator=
        DimIndexSeq1.add_new(
            (0x0020, 0x9164), 'UI', DimOrgSeq.DimensionOrganizationUID[0])
        DimIndexSeq1.DimensionDescriptionLabel = 'Third Spatial dimension'

        DimIndexSeq2 = Dataset()
        DimIndexSeq2.DimensionIndexPointer = (0x0018, 0x9087) # Diffusion b-value
        # DimIndexSeq2.DimensionIndexPrivateCreator=
        # DimIndexSeq2.FunctionalGroupPointer=
        # DimIndexSeq2.FunctionalGroupPrivateCreator=
        DimIndexSeq2.add_new(
            (0x0020, 0x9164), 'UI', DimOrgSeq.DimensionOrganizationUID[1])
        DimIndexSeq2.DimensionDescriptionLabel = 'Fourth dimension (diffusion b value)'
        ds.DimensionIndexSequence = Sequence([DimIndexSeq2, DimIndexSeq1])
    else:
        DimIndexSeq1 = Dataset()
        # Image position patient 20,32 or 20,12
        DimIndexSeq1.DimensionIndexPointer = (0x0020, 0x0032)
        # #DimIndexSeq1.DimensionIndexPrivateCreator=
        # #DimIndexSeq1.FunctionalGrouppointer=
        # #DimIndexSeq1.FunctionalGroupprivateCreator=
        DimIndexSeq1.add_new(
            (0x0020, 0x9164), 'UI', [DimensionOrganizationUID])
        DimIndexSeq1.DimensionDescriptionLabel = 'Third Spatial dimension'
        ds.DimensionIndexSequence = Sequence([DimIndexSeq1])

        # Module: Image Pixel (mandatory)
        # Reference: DICOM Part 3: Information Object Definitions C.7.6.3
        # ds.Rows                                                              # 0028,0010 Rows (mandatory)
        # ds.Columns                                                           # 0028,0011 Columns (mandatory)
        # ds.BitsStored                                                        # 0028,0101 (mandatory)
        # ds.HighBit                                                           # 0028,0102 (mandatory)
        # ds.PixelRepresentation                                               # 0028,0103 Pixel Representation (mandatory)
        # ds.PixelData                                                        #
        # 7fe0,0010 Pixel Data (mandatory)

    FrameContentSequence = Dataset()
    # FrameContentSequence.FrameAcquisitionNumber = '1' #fid_properties['slice_no']
    # FrameContentSequence.FrameReferenceDateTime
    # FrameContentSequence.FrameAcquisitionDateTime
    # FrameContentSequence.FrameAcquisitionDuration
    # FrameContentSequence.CardiacCyclePosition
    # FrameContentSequence.RespiratoryCyclePosition
    # FrameContentSequence.DimensionIndexValues = 1 #islice # fid_properties['array_no']
    #FrameContentSequence.TemporalPositionIndex = 1
    FrameContentSequence.StackID = [str(1)]  # fourthdimid
    FrameContentSequence.InStackPositionNumber = [int(1)]  # fourthdimindex
    FrameContentSequence.FrameComments = ''  # fid_properties['filetext']
    FrameContentSequence.FrameLabel = 'DimX'
    ds.FrameContentSequence = Sequence([FrameContentSequence])

    return ds, fid_size_matrix, tmatrix
# end ParseFID


import errno


def mkdir_or_cleandir(path, args):
    """
    https://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
For Python >= 3.2, os.makedirs has an optional third argument exist_ok that, when true, enables the mkdir -p functionality  - unless mode is provided and the existing directory has different permissions than the intended ones; in that case, OSError is raised as previously.
    """
    try:
        if os.path.isdir(path):
            if os.listdir(path) == []:
                if args.verbose:
                    print "Output path exists and is empty"
            else:
                if args.verbose:
                    print "Cleaning output path"
                import shutil
                shutil.rmtree(path, ignore_errors=True)
                os.makedirs(path)
        else:
            if args.verbose:
                print "Creating dir ", path
            os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            print "Cannot create dir ", path
            raise


def SaveFIDtoDicom(ds, procpar, image_data, fid_properties, M, args, outdir):
    """
    Multi-dimension (3D) export of FID format image data and metadata to DICOM
    Wrapper for Save3dFIDtoDicom

    :param ds:         Baseline pydicom dataset
    :param image_data: Pixel data of image
    :param procpar:    General header info
    :param fid_properties: Local fid header info
    :param M:          Image transformation matrix
    :param args:       argparse struct
    :param outdir:     output directory
    :return ds:        Return the pydicom dataset
    """

    if args.verbose:
        print "FID export of complex image to DICOM ", outdir

    if args.magnitude:
        if args.verbose:
            print "Exporting magnitude ... POP!! POP!!"
        orig_imagetype = ds.ImageType[2]
        if ds.ImageType[2] == "PHASE MAP":
            ds.ImageType[2] = "MAGNITUDE"
        ds.ComplexImageComponent = 'MAGNITUDE'
        voldata = np.abs(image_data)  # Magnitude
        outdir1 = outdir  # re.sub('.dcm', '-mag.dcm', outdir)
        if not os.path.isdir(outdir1):
            mkdir_or_cleandir(outdir1, args)
        Save3dFIDtoDicom(ds, procpar, voldata, fid_properties, M,
                         args, outdir1)
        ds.ImageType[2] = orig_imagetype

    if args.phase:
        if args.verbose:
            print "Exporting phase ..."
        voldata = np.angle(image_data)  # Phase
        orig_imagetype = ds.ImageType[2]
        ds.ImageType[2] = "PHASE MAP"
        ds.ComplexImageComponent = 'PHASE'
        voldata = np.angle(image_data)
        outdir1 = re.sub('.dcm', '-pha.dcm', outdir)
        if not os.path.isdir(outdir1):
            mkdir_or_cleandir(outdir1, args)
        Save3dFIDtoDicom(ds, procpar, voldata, fid_properties, M,
                         args, outdir1)
        ds.ImageType[2] = orig_imagetype
        ds.ComplexImageComponent = 'MAGNITUDE'
    if args.realimag:
        if args.verbose:
            print "Exporting real and imag ..."
        voldata = np.real(image_data)
        ds.ComplexImageComponent = 'REAL'
        outdir1 = re.sub('.dcm', '-real.dcm', outdir)
        if not os.path.isdir(outdir1):
            mkdir_or_cleandir(outdir1, args)
        Save3dFIDtoDicom(
            ds, procpar, voldata, fid_properties, M, args, outdir1)
        voldata = np.imag(image_data)
        ds.ComplexImageComponent = 'IMAGINARY'
        outdir1 = re.sub('.dcm', '-imag.dcm', outdir)
        if not os.path.isdir(outdir1):
            mkdir_or_cleandir(outdir1, args)
        Save3dFIDtoDicom(
            ds, procpar, voldata, fid_properties, M, args, outdir1)
        ds.ComplexImageComponent = 'MAGNITUDE'

# end SaveFIDtoDicom


def Save3dFIDtoDicom(ds, procpar, voldata, fid_properties, M, args, outdir):
    """
    Multi-dimension (3D) export of FID format image data and metadata to DICOM

    :param ds:  Baseline pydicom dataset
    :param image_data: Pixel data of image
    :param procpar: General header info
    :param fid_properties: Local fid header info
    :param M: Image transformation matrix
    :param args: argparse struct
    :param outdir: output directory

    :return ds: Return the pydicom dataset
    """
    if np.iscomplexobj(voldata):
        print "Error: Volume data is complex prior to saving as DICOM. Forcing magnitude as output image."
        voldata = np.abs(voldata)
    # Calculate the max and min throughout all dataset values;
    # calculate the intercept and slope for casting to UInt16
    ds, voldata = RescaleFIDImage(ds, voldata, args)


#    if not os.path.isdir(outdir):
#        if args.verbose:
#            print "Save3dFIDtoDicom output path has not been created."
    mkdir_or_cleandir(outdir, args)

    # if procpar['recon'] == 'external':
    #
    #        pdb.set_trace()
    # if procpar['recon'] == 'external' and fid_properties['rank'] == 3:
    #   if procpar['seqfil'] == "epip":
    #       print "Transposing external recon 3D"
    #       voldata = np.transpose(voldata,(1,2,0)) # 1,2,0
    #   if procpar['seqfil'] == "fse3d":
    #       print "Transposing external recon 3D"
    #       voldata = np.transpose(voldata,(2,0,1)) # 0,2,1 works
    # readprocpar.m procpar('nD') == 3
    #        acq.FOVcm = [pps.lro pps.lpe pps.lpe2];
    #        acq.dims = [pps.nf pps.np/2 pps.nv2];
    #        acq.voxelmm = acq.FOVcm./acq.dims*10;

    if args.verbose:
        print "Image data shape: ", str(voldata.shape)
        print "Vol data shape: ", voldata.shape
        print "fid properties matrix: ", fid_properties['dims']
        print "Slice points: ", fid_properties['dims'][0] * fid_properties['dims'][1]
    #  slice_data = np.zeros_like(np.squeeze(image_data[:,:,1]))
    #   if 'ne' in procpar.keys():

    range_max = voldata.shape[2]  # fid_properties['dims'][2]
    # fid_properties['dims'][0]*fid_properties['dims'][1]
    num_slicepts = voldata.shape[0] * voldata.shape[1]
    # if procpar['recon'] == 'external' and procpar['seqfil'] == 'fse3d' and fid_properties['rank'] == 3:
    #    range_max = fid_properties['dims'][1]
    #    num_slicepts = fid_properties['dims'][0]*fid_properties['dims'][2]
    #    ds.Columns = fid_properties['dims'][2]
    #    ds.Rows = fid_properties['dims'][0]
    # FIXME FSE3d still producing bad dicoms

    if args.verbose:
        print "Columns ", ds.Columns, " Rows ", ds.Rows
        print "Range max and no slice points: ", range_max, num_slicepts

    # Export dicom to file
#    if MRAcquisitionType == '3D':

#    else:


#    if ds.ImageType[2]=="MULTIECHO":
#        print ds.FrameContentSequence[0].StackID, ds.FrameContentSequence[0].StackID[0]
#        print type(ds.FrameContentSequence[0].StackID), type(ds.FrameContentSequence[0].StackID[0])
#        ds.FrameContentSequence[0].StackID = str(int(ds.FrameContentSequence[0].StackID[0]) + 1)
#        if args.verbose:
# print "Incrementing volume StackID ", ds.FrameContentSequence[0].StackID

    ds.FrameContentSequence[0].StackID = [str(1)]
    if voldata.ndim == 5:
        if args.verbose:
            print "Saving 5D image "
        for echo in xrange(0, voldata.shape[4]):
            ds.EchoNumber = str(echo + 1)
            for n in xrange(0, voldata.shape[3]):
                # Indexing in np matrix begins at 0, fid/dicom filenames
                # begin at 1
                for islice in xrange(0, range_max):
                    if args.verbose:
                        print "Reshape volume slice to 1D array ", voldata.shape[0] * voldata.shape[1], num_slicepts

                    slice_data = np.reshape(np.squeeze(voldata[:, :, islice,
                                                               n, echo]),
                                            (num_slicepts, 1))
                    # Convert Pixel data to string
                    ds.PixelData = slice_data.tostring()
                    new_filename = "slice%03dimage%03decho%03d.dcm" % (
                        islice + 1, n + 1, echo + 1)

# if procpar['recon'] == 'external' and fid_properties['rank'] == 3 and
# procpar['seqfil'] == 'fse3d':
                    pos = np.matrix([[0], [0], [islice], [1]])
#        else:
#            pos = np.matrix([[0],[0],[islice],[1]])

                    # Get slice's position
                    Pxyz = M * pos
                    ds.ImagePositionPatient = [
                        str(Pxyz[0, 0]), str(Pxyz[1, 0]), str(Pxyz[2, 0])]

                    ds.FrameContentSequence[0].InStackPositionNumber = [
                        int(islice)]  # THIRDdimindex
                    ds.FrameContentSequence[
                        0].TemporalPositionIndex = int(ds.EchoNumber)

                    # Save DICOM
                    if args.verbose:
                        print "Saving 2D slice in 5D DICOM slice %d, image %d, echo %d" % (islice, n, echo)

                    ds.save_as(os.path.join(outdir, new_filename))
            # Increment StackID
            if args.verbose:
                print "Incrementing volume StackID (current) ", ds.FrameContentSequence[0].StackID
            ds.FrameContentSequence[0].StackID = str(
                int(ds.FrameContentSequence[0].StackID[0]) + 1)
    else:
        if args.verbose:
            print "Saving 3D image "
        for islice in xrange(0, range_max):
            # Reshape volume slice to 1D array
            slice_data = np.reshape(np.squeeze(voldata[:, :, islice]),
                                    (num_slicepts, 1))
            # Convert Pixel data to string
            ds.PixelData = slice_data.tostring()
            new_filename = "slice%03dimage001echo001.dcm" % (islice + 1)

# if procpar['recon'] == 'external' and fid_properties['rank'] == 3 and
# procpar['seqfil'] == 'fse3d':
            pos = np.matrix([[0], [0], [islice], [1]])
#        else:
#            pos = np.matrix([[0],[0],[islice],[1]])

            # Get slice's position
            Pxyz = M * pos
            ds.ImagePositionPatient = [
                str(Pxyz[0, 0]), str(Pxyz[1, 0]), str(Pxyz[2, 0])]
            ds.FrameContentSequence[0].InStackPositionNumber = [int(islice)]
            # THIRDdimindex
            ds.FrameContentSequence[0].TemporalPositionIndex = int(0)

            # Save DICOM
            if args.verbose:
                print "Saving 2D slice in 3d DICOM slice %d " % (islice)
            ds.save_as(os.path.join(outdir, new_filename))

    return ds
# end Save3dFIDtoDicom


def Save2dFIDtoDicom(ds, image_data, fid_properties, outdir, filename):
    """
    Export 2D image and metadata to DICOM

    :param ds: dicom dataset struct
    :param image_data: image data
    :param fid_properties: FID header info

    """

    # Common export format
    ds.PixelData = image_data.tostring()         # (7fe0,0010) Pixel Data

    if ds.ImageType[2] == "ASL":
        image_number = fid_properties['array_index']
        if fid_properties["asltag"] == 1:               # Labelled
            new_filename = "slice%03dimage%03decho%03d.dcm" % (
                fid_properties['slice_no'], image_number, 1)
        elif fid_properties["asltag"] == -1:  # Control
            new_filename = "slice%03dimage%03decho%03d.dcm" % (
                fid_properties['slice_no'], image_number, 2)
        else:                                            # Unknown
            new_filename = "slice%03dimage%03decho%03d.dcm" % (
                fid_properties['slice_no'], image_number, 3)
        dcmfilename = os.path.join(outdir, new_filename)
    else:
        dcmfilename = os.path.join(
            outdir, os.path.splitext(filename)[0] + '.dcm')
    # Save DICOM
    ds.save_as(dcmfilename)

    return 1
# end SaveSave2dFIDtoDicom


def SaveKspace(ksp, args):
    """SaveKspace Save the k-space data to MAT file

    """
    if not args.outputdir:
        outdir = os.path.splitext(args.inputdir)[0]
        if not outdir.find('.img'):
            outfile = outdir + '-ksp.mat'
        else:
            (dirName, imgdir) = os.path.split(outdir)
            while imgdir == '':
                (dirName, imgdir) = os.path.split(dirName)

            (ImgBaseName, ImgExtension) = os.path.splitext(imgdir)
            outfile = os.path.join(dirName, ImgBaseName + '-ksp.mat')
    else:
        if not outdir.find('.dcm'):
            outfile = args.outputdir + '-ksp.mat'
        else:
            (dirName, imgdir) = os.path.split(outdir)
            while imgdir == '':
                (dirName, imgdir) = os.path.split(dirName)

            (ImgBaseName, ImgExtension) = os.path.splitext(imgdir)
            outdir = os.path.join(dirName, ImgBaseName + '-ksp.mat')
    # Force ksp type to be two 32-bit floats
    ksp = ksp.astype(np.complex64)
    scipy.io.savemat(outfile, mdict={'ksp': ksp}, do_compression=True)
# end SaveSaveKspace

if __name__ == "__main__":

    import nibabel as nib

    parser = argparse.ArgumentParser(usage=' ReadFID.py -i "Input FID directory"',
                                     description='''agilent2dicom is an FID to
                                     Enhanced MR DICOM converter from MBI. ReadFID
                                     takes header info from fid files and adjusts
                                     the dicom dataset *ds* then rescales the image
                                     data.''')
    parser.add_argument(
        '-i', '--inputdir', help='''Input directory name. Must be an Agilent FID
        image directory containing procpar and fid files''', required=True)
    parser.add_argument(
        '-o', '--outputdir', help='Output directory name for DICOM files.')
    parser.add_argument(
        '-m', '--magnitude', help='Magnitude component flag.', action="store_true")
    parser.add_argument(
        '-p', '--phase', help='Phase component flag.', action="store_true")
    parser.add_argument(
        '-s', '--sequence', help='Sequence type (one of Multiecho, Diffusion, ASL).')
    parser.add_argument('-a', '--axis_order', help='Axis order eg 1,0,2.')
    parser.add_argument(
        '-v', '--verbose', help='Verbose comments.', action="store_true")
    args = parser.parse_args()

    # import ProcparToDicomMap as ptd

    procpar, procpartext = ReadProcpar.ReadProcpar(
        os.path.join(args.inputdir, 'procpar'))
    ds, MRAcq_type = ProcparToDicomMap.ProcparToDicomMap(procpar, args)
    print "Rows: ", ds.Rows, " Columns: ", ds.Columns

    files = os.listdir(args.inputdir)
    fidfiles = [f for f in files if f.endswith('fid')]
    if args.verbose:
        print "Number of FID files ", len(fidfiles)

    # for filename in fidfiles:
    if args.verbose:
        print "Reading FID"
    filename = fidfiles[len(fidfiles) - 1]
    procpar, fid_header, dims, data_real, data_imag = readfid(
        args.inputdir, procpar, args)
    if args.verbose:
        print "Dims: ", dims, " Echoes: ", fid_header['nEchoes'], " Channels: ", fid_header['nChannels']
    affine = np.eye(4)
    # # affine[:3, :3]= np.arange(9).reshape((3, 3))
    # raw_data = nib.Nifti1Image(image_data_real, affine)
    # nib.save(raw_data, 'raw_data.nii.gz')

    # Change dicom for specific FID header info
    ds, matsize, ImageTransformationMatrix = ParseFID(
        ds, fid_header, procpar, args)

    if os.path.exists('raw_image_00.nii.gz'):
        if args.verbose:
            print "Getting Original image (reconstruction)"
        nii = nib.load('raw_image_00.nii.gz')

        image = np.empty(
            [nii.shape[0], nii.shape[1], nii.shape[2], 1, 3], dtype=np.complex64)
        image[:, :, :, 0, 0] = nii.get_data()
        nii = nib.load('raw_image_01.nii.gz')
        image[:, :, :, 0, 1] = nii.get_data()
        nii = nib.load('raw_image_02.nii.gz')
        image[:, :, :, 0, 2] = nii.get_data()
    else:
        if args.verbose:
            print "Computing Original image (reconstruction)"
        image, ksp = recon(
            procpar, dims, fid_header, data_real, data_imag, args)

        if args.verbose:
            print "Saving raw image"
        if image.ndim == 5:
            for i in xrange(0, image.shape[4]):
                raw_image = nib.Nifti1Image(
                    np.abs(image[:, :, :, 0, i]), affine)
                nib.save(raw_image, 'raw_image_0' + str(i) + '.nii.gz')
        else:
            raw_image = nib.Nifti1Image(np.abs(image), affine)
            nib.save(raw_image, 'raw_image.nii.gz')
        #    raw_ksp = nib.Nifti1Image(np.abs(ksp), affine)
        #    nib.save(raw_ksp, 'raw_ksp.nii.gz')

    if args.axis_order:
        timage = RearrangeImage(image, args.axis_order, args)
        tr_image = nib.Nifti1Image(np.abs(timage[:, :, :, 0, 0]), affine)
        nib.save(tr_image, 'tr_image.nii.gz')

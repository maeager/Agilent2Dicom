#!/usr/bin/env python
"""ReadProcpar is used to read Agilent FDF config files

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
import os
import sys
import re
import argparse



def ReadProcpar(procparfilename):
    """
READPROCPAR - Read procpar file and return procpar dictionary and text

Procpar element format

First line: name subtype basictype maxvalue minvalue stepsize Ggroup Dgroup protection active intptr
  name: string
  subtype: 0 (undefined), 1 (real), 2 (string), 3 (delay), 4 (flag), 5 (frequency), 6 (pulse), 7 (integer).
  basictype: 0 (undefined), 1 (real), 2 (string).
  maxvalue: the maximum value that the parameter can contain, or an index to a maximum 
            value in the parameter parmax (found in /vnmr/conpar). Applies to both
            string and real types of parameters.
  minvalue: the minimum value that the parameter can contain or an index to a minimum 
            value in the parameter parmin (found in /vnmr/conpar). Applies to real types
            of parameters only.
  stepsize: a real number for the step size in which parameters can be entered or index
            to a step size in the parameter parstep (found in /vnmr/conpar). If stepsize
            is 0, it is ignored. Applies to real types only.
  Ggroup: 0 (ALL), 1 (SAMPLE), 2 (ACQUISITION), 3 (PROCESSING), 4 (DISPLAY), 5 (SPIN).
  Dgroup: The specific application determines the usage of this integer.
  protection: a 32-bit word made up of the following bit masks, which are summed to form
              the full mask:
                  0  1    Cannot array the parameter
                  1  2    Cannot change active/not active status
                  2  4    Cannot change the parameter value
                  3  8    Causes _parameter macro to be executed (e.g., if parameter is named sw, the macro _sw is executed when sw is changed)
                  4  16   Avoids automatic redisplay
                  5  32   Cannot delete parameter
                  6  64   System parameter for spectrometer or data station
                  7  128  Cannot copy parameter from tree to tree
                  8  256  Cannot set array parameter
                  9  512  Cannot set parameter enumeral values
                  10 1024 Cannot change the parameter's group
                  11 2048 Cannot change protection bits
                  12 4096 Cannot change the display group
                  13 8192 Take max, min, step from /vnmr/conpar parameters parmax, parmin, parstep.
  active: 0 (not active), 1 (active).
  intptr: not used (generally set to 64).

if basictype=1,
  Second line: numvalues value1 [value2] [value3] ...

if basictype=2,
  Second line: numvalues value1 
  [Third line: value2]
  [Fourth line: value3]
  ...

Last line: 0, or if subtype=4 (flag), an array of possible flag values formatted a:
   numvalues flag1 [flag2] [flag3]

Notes:
1. All strings are enclosed in double quotes.
2. Floating point values have been found associated with the subtype 'integer',   therefore it is advised to read these values as floats.


    """

    f = open(procparfilename, 'r')
    line = f.readline()
    procpar = {}
    while line != '':
        # parse first line in property
        tokens = line.split()
        propname = tokens[0]
        propsubtype = tokens[1]
        proptype = tokens[2]
        # parse second line in property: [number of values] [first value] ...
        line = f.readline()
        tokens = line.strip().split(None, 1)
        propnumvalues = int(tokens[0])
        # handle property values
        if proptype == '1': # real number
            if propnumvalues == 1:
                propvalue = float(tokens[1])
            else:
                propvalue = map(float, tokens[1].split())
        elif proptype == '2': # string
            if propnumvalues == 1:
                propvalue = tokens[1].strip('"')
            if propnumvalues > 1:
                propvalue = [tokens[1].strip('"')]
                for i in range(2, propnumvalues + 1):
                    propvalue.append(f.readline().strip('"\n"'))
        line = f.readline() # last line in property
        line = f.readline() # next property        
        lastprop = propvalue
        procpar[propname] = propvalue
    f.seek(0)
    procpartext = f.readlines()
    return (procpar, procpartext)
# ReadProcpar


def ProcparInfo(procpar):
    """ProcparInfo create basic information from procpar for displaying

       This method is primarily for the Agilent2DicomAppQt GUI.
    """
    header = dict()
    header['Protocol_Name'] = procpar['seqfil']
    header['StudyID'] = procpar['name']
    header['SeriesDescription'] = procpar['comment']
    header['Mode'] = '%dD' % procpar['nD']
    if procpar['nD'] == 2:
        header['FOV_cm'] = [procpar['lro'], procpar['lpe']]
        header['Dimensions'] = [procpar['nf']/procpar['ns'], procpar['np']/2, procpar['ns']]
        #if len(procpar['thk']) > 1:
        #    print "procpar thk size greater than 1"
        header['Voxel_Res_mm'] = [procpar['lro']*10/header['Dimensions'][0],
                                  procpar['lpe']*10/header['Dimensions'][1], procpar['thk']*10]
    elif procpar['nD'] == 3:
        header['FOV_cm'] = [procpar['lro'], procpar['lpe'], procpar['lpe2'], procpar['lpe3']]
        header['Dimensions'] = [procpar['nf'], procpar['np']/2, procpar['nv2']]
        header['Voxel_Resolution_mm'] = [procpar['lro']*10/header['Dimensions'][0],
                                         procpar['lpe']*10/header['Dimensions'][1],
                                         procpar['lpe2']*10/header['Dimensions'][2]]
    header['Volumes'] = procpar['volumes']
    header['NumberOfEchoes'] = procpar['ne']
    header['FOV_ORIENTATION'] = [procpar['orient'], procpar['psi'],
                                 procpar['phi'], procpar['theta']]
    header['GRADIENTS'] = [procpar['gcoil'], procpar['gro'],
                           procpar['gpe']]  #, procpar['gss'], procpar['gspoil'], procpar['rewind']]
    header['ACQ_CONTROL'] = 'seqcon:' + procpar['seqcon'] + 'nD:' + procpar['nD']
    + 'ni:' + procpar['ni'] + 'nf:' + procpar['nf'] + 'cf:' + procpar['cf']
    + 'ne:' + procpar[', ne'] + 'ns:' + procpar['ns'] + 'np:' + procpar[', np']
    + 'nv:' + procpar['nv'] + 'nv2:' + procpar['nv2'] + 'nv3:' + procpar['nv3']
    rcvrs = re.findall('y', procpar['rcvrs'])
    if rcvrs:
        header['NumberOfChannels'] = len(rcvrs);
    else:
        header['NumberOfChannels'] = 1
                            
    
    dims=[1, 1, 1, 1, 1]
    # validate dimensions
    if procpar['nD'] == 2:
        dims[0] = procpar['np']/2        # num phase encode lines / 2
        dims[1] = procpar['nf']/procpar['ns'] # num frequency lines acquired / # echoes
        dims[2] = procpar['ns']          # if 2D, num slices, else ni2
        if procpar['ni2'] > 1:           # fse3d sequence has nD == 2, but is a 3d acquisition???
            dims[2] = procpar['ni2']
    elif procpar['nD'] == 3:
        dims[0] = procpar['np']/2        # NUM phase encode lines / 2
        dims[1] = procpar['nf']/procpar['ne'] # num frequency lines acquired / # echoes
        dims[2] = procpar['ni2']
    dims[3] = header['NumberOfChannels']
    dims[4] = header['NumberOfEchoes']
    header['Dims2'] = dims
    header['DimXYZ'] = [procpar[procpar['dimX']], procpar[procpar['dimY']],
                        procpar[procpar['dimZ']]]
    header['PosXYZ'] = [procpar[procpar['posX']], procpar[procpar['posY']],
                        procpar[procpar['posZ']]]
    header['VoxelSize'] = [procpar['vox1'], procpar['vox2'], procpar['vox3']]
    header['TR'] = str(procpar['tr']*1000.0)
    if 'te' in procpar.keys():
        header['TE'] = str(procpar['te']*1000.0)
    return header

if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage='ReadProcpar -i "Input FDF directory"',
                                     description='ReadProcpar is part of agilent2dicom, an FDF to Enhanced MR DICOM converter from MBI.')
    parser.add_argument('-i', '--inputdir', help='Input directory name. Must be an Agilent FDF directory', required=True);
    parser.add_argument('-s', '--show', help='Show procpar info.', action="store_true");
    parser.add_argument('-v', '--verbose', help='Verbose.', action="store_true");
    args = parser.parse_args()


    procpar, procpartext = ReadProcpar(os.path.join(args.inputdir, 'procpar'))
    if not args.show:
        print procpar
        print "Basic info:"

    p = ProcparInfo(procpar)
    #print '\n'.join(p)
    print '\n'.join("%s:  %r" % (key, val) for (key, val) in p.iteritems())

#!/usr/bin/env python

"""ReadProcpar is used to read Agilent FDF config files

   (c) 2014  Michael Eager (michael.eager@monash.edu)
"""
import os,sys
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

if basictype = 1, 
  Second line: numvalues value1 [value2] [value3] ... 

if basictype = 2, 
  Second line: numvalues value1  
  [Third line: value2]
  [Fourth line: value3]
  ...

Last line: 0, or if subtype = 4 (flag), an array of possible flag values formatted as
   numvalues flag1 [flag2] [flag3]

Notes:
1. All strings are enclosed in double quotes.
2. Floating point values have been found associated with the subtype 'integer', 
   therefore it is advised to read these values as floats.


    """

    f = open(procparfilename,'r')
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
        tokens = line.strip().split(None,1)
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
                for i in range(2,propnumvalues+1):
                    propvalue.append(f.readline().strip('"\n"'))
        line = f.readline() # last line in property
        line = f.readline() # next property        
        lastprop = propvalue
        procpar[propname] = propvalue
    f.seek(0)
    procpartext = f.readlines()
    return (procpar, procpartext)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage=' ReadProcpar -i "Input FDF directory" ',description='ReadProcpar is part of agilent2dicom, an FDF to Enhanced MR DICOM converter from MBI.')
    parser.add_argument('-i','--inputdir', help='Input directory name. Must be an Agilent FDF directory',required=True);
    args = parser.parse_args()


    procpar, procpartext = ReadProcpar(args.inputdir + '/procpar')
    print procpar
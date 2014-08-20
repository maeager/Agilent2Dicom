#!/usr/bin/env python

"""ParseFDF is used to parse Agilent FDF image files and correct errors in procpar

   (c) 2014  Michael Eager (michael.eager@monash.edu)
"""

import os,sys,math
import numpy
import argparse

# Numpy recast to int16 with range (-32768 or 32767)
UInt16MaxRange = 65533


def RescaleFDF(procpar,args):
    """RescaleFDF
     Calculate the max and min throughout all fdf iles in dataset;
     calculate the intercept and slope for casting to UInt16

    :param procpar: FDF procpar dictionary
    :param args:    Input arguments
    :param slope_factor: Default divisor for data slope is the uint16 range maximum
    :returns: RescaleIntercept Minimum data value in image. RescaleSlope data range divided by uint16 range.
    """

    # Read in data from all files to determine scaling
    datamin = float("inf")
    datamax = float("-inf")

    # RescaleSlope of phase imgs set to [-pi,pi]
    if args.phase and 'imPH' in procpar.keys() and procpar["imPH"] == 'y':
        datamin = -math.pi
        datamax = math.pi
    else:
        for filename in fdffiles:
            fdf_properties, data = ReadFDF(args.inputdir + '/' + filename)
            datamin = numpy.min([datamin,data.min()])
            datamax = numpy.max([datamax,data.max()])


    RescaleIntercept = datamin
    # Numpy recast to int16 with range (-32768 or 32767)
    if ds.PixelRepresentation == 0:
        slope_factor=UInt16MaxRange
    else:
        slope_factor=32767
    RescaleSlope = (datamax - datamin) / slope_factor # 65533 # / 32767    

    return RescaleIntercept,RescaleSlope


if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage=' procpartodicommapping -i "Input FDF directory"',description='agilent2dicom is an FDF to Enhanced MR DICOM converter from MBI. Version ' + VersionNumber)
    parser.add_argument('-i','--inputdir', help='Input directory name. Must be an Agilent FDF image directory containing procpar and *.fdf files',required=True);
    parser.add_argument('-m','--magnitude', help='Magnitude component flag.',action="store_true");
    parser.add_argument('-p','--phase', help='Phase component flag.',action="store_true");
    parser.add_argument('-v','--verbose', help='Verbose.',action="store_true");
    args = parser.parse_args()
    
    from ReadFDF import ReadFDF
    from ReadProcpar import ReadProcpar
#    from ProcparToDicomMap import ProcparToDicomMap
#    from ProcparToDicomMap import CreateUID

    procpar, procpartext = ReadProcpar(args.inputdir+'/procpar')
#    ds,MRAcq_type,SEQUENCE,BValue = ProcparToDicomMap(procpar, args)

    RescaleIntercept,RescaleSlope = RescaleFDF(procpar,args)

    print "Intercept: ", RescaleIntercept
    print "Slope: ", RescaleSlope

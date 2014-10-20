#!/usr/bin/env python

"""ParseFDF is used to parse Agilent FDF image files and correct errors in procpar

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

import os,sys,math
import re
import numpy
import argparse

# Numpy recast to int16 with range (-32768 or 32767)
UInt16MaxRange = 65533
Int16MaxRange = 32767

from ReadProcpar import *
import ReadFDF


def ShortenFloatString(val,origin):
    """Shorten float to 16 element string
    """
    if numpy.iscomplex(val):
        print "Error: ShortenFloatString given complex number, Source: %s" % origin
        val=numpy.abs(val)
    if val < 0:
        print "Negative val in shorten string"
    if len(str(val))>=15:
        stringval=str(val)
        if re.search('e',stringval): #scientific notation
            pos = stringval.index('e')
            stripped_val = stringval[:pos - (len(stringval)-15)]+stringval[pos:]
        else: # normal float 
            stripped_val = stringval[:15]
        print "Cropping float string from ", str(val), " to ", stripped_val
        return stripped_val[:15]
    else:
        return str(val)



def FindScale(fdffiles,ds,procpar,args):
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
    if (len(ds.ImageType)>=3 and ds.ImageType[2] == "PHASE MAP") or \
            (hasattr(ds,'ComplexImageComponent') and ds.ComplexImageComponent == 'PHASE'):
        # this implies either args.phase is on or procpar['imPH']=='y'
        datamin = -math.pi
        datamax = math.pi
    else:
        for filename in fdffiles:
            fdf_properties, data = ReadFDF.ReadFDF(os.path.join(args.inputdir, filename))
            datamin = numpy.min([datamin,data.min()])
            datamax = numpy.max([datamax,data.max()])


    RescaleIntercept = datamin
    # Numpy recast to int16 with range (-32768 or 32767)
    if ds.PixelRepresentation == 0:
        slope_factor = UInt16MaxRange
    else:
        slope_factor = Int16MaxRange
    RescaleSlope = (datamax - datamin) / slope_factor 

    return RescaleIntercept,RescaleSlope

    
    
def RescaleImage(ds,image_data,RescaleIntercept,RescaleSlope,args):

    if args.verbose:
        print "Rescale data to uint16"
        print "Intercept: ", RescaleIntercept, "  Slope: ", RescaleSlope
        print "Current data min: ", image_data.min(), " max ", image_data.max()
    image_data = (image_data - RescaleIntercept) / RescaleSlope
    image_data_int16 = image_data.astype(numpy.int16)

    ## Adjusting Dicom parameters for rescaling
    # Rescale intercept string must not be longer than 16
    ds.RescaleIntercept = ShortenFloatString(RescaleIntercept,"RescaleIntercept")  #(0028,1052) Rescale Intercept
    # Rescale slope string must not be longer than 16
    ds.RescaleSlope = ShortenFloatString(RescaleSlope,"RescaleSlope")  #(0028,1053) Rescale Slope

    return ds,image_data_int16



if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage=' procpartodicommapping -i "Input FDF directory"',description='agilent2dicom is an FDF to Enhanced MR DICOM converter from MBI. Version ' + VersionNumber)
    parser.add_argument('-i','--inputdir', help='Input directory name. Must be an Agilent FDF image directory containing procpar and *.fdf files',required=True);
    parser.add_argument('-m','--magnitude', help='Magnitude component flag.',action="store_true");
    parser.add_argument('-p','--phase', help='Phase component flag.',action="store_true");
    parser.add_argument('-v','--verbose', help='Verbose.',action="store_true");
    args = parser.parse_args()
    

    import ReadProcpar as rp
    import ProcparToDicomMap as ptd
#    from ProcparToDicomMap import CreateUID

    procpar, procpartext = rp.ReadProcpar(os.path.join(args.inputdir,'procpar'))
    ds,MRAcq_type = ptd.ProcparToDicomMap(procpar, args)

    files = os.listdir(args.inputdir)
    fdffiles = [ f for f in files if f.endswith('.fdf') ]

    RescaleIntercept,RescaleSlope = FindScale(fdffiles,ds,procpar,args)

    print "Intercept: ", RescaleIntercept
    print "Slope: ", RescaleSlope

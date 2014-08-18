#!/usr/bin/env python

"""ReadFDF is used to read Agilent FDF image files

   (c) 2014  Michael Eager (michael.eager@monash.edu)
"""

import os,sys
import numpy
import argparse

def ReadFDF(fdffilename):
    """
    READFDF - Read FDF file and return properties derived from fdf header and image data
    
    Infomation on header parameters can be found in the "Agilent VNMRJ 3.2 User Programming User Guide"
    :param fdffilename: Name string of FDF file
    :return fdf_properties: Label/value dictionary of FDF header properties
    :return image data: 1D float array of pixel data 
    """
    f = open(fdffilename,'r')
    #print fdffilename
    # Read in dataset properties
    fdftext = ''
    fdf_properties = dict()
    line = f.readline()
    while line[0] != '\x0c': # or line[0] != '\x00':
        fdftext = fdftext + line
        if line[0] == '#':
            line = f.readline()
            continue
        # print line
        if line.find("=") == -1 and line[0] != '\n':
            # print 'Unknown header line in fdf.'
            continue
        tokens = line.strip(' ;\n').split(' ',1)
        tokens = tokens[1].strip().split('=')
        tokens[0] = tokens[0].strip(' *[]')
        tokens[1] = tokens[1].strip()
        if tokens[1][0] == '{' and tokens[1][-1] == '}':
            tokens[1] = '[' + tokens[-1].strip('{}') + ']'
        exec('fdf_properties["' + tokens[0] + '"] = ' + tokens[1])
        line = f.readline()
    fdf_properties['filename'] = fdffilename
    fdf_properties['filetext'] = fdftext
    
    # Find NULL indicating start of image
    c = f.read(1)
    while c != '\x00':
        c = f.read(1)
        
    # read in data
    if fdf_properties['storage'] == "integer":
        dt = "int"
    elif fdf_properties['storage'] == "float":
        dt = "float"
    else:
        print "Error: unrecognised fdf header storage value"
        sys.exit(1)
        
    dt = dt + str(fdf_properties['bits'])
    
    data = numpy.fromfile(f,dtype=dt)
    
    f.close()
    
    return (fdf_properties, data)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage=' ReadFDF -i "Input FDF file" ',description='ReadFDF is part of agilent2dicom, an FDF to Enhanced MR DICOM converter from MBI.')
    parser.add_argument('-i','--inputfile', help='Input file name. Must be an Agilent FDF file',required=True);
    args = parser.parse_args()

    fdf_properties, image_data = ReadFDF(args.inputfile)

    print 'Image_data shape:', str(image_data.shape)
    print "Origin: ", fdf_properties['origin']
    print "ROI: ",fdf_properties['roi']
    print "Orientation: ", fdf_properties['orientation']
    print "StudyID: ", fdf_properties['studyid']
    print "Slice ", fdf_properties['slice_no'], " of ", fdf_properties['slices']
    print "Matrix size", fdf_properties['matrix']

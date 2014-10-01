# README #

  Copyright 2014 Michael Eager  (michael.eager@monash.edu).

  This file is part of the Agilent2Dicom package.

  The Agilent2Dicom package is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.  

  The Agilent2Dicom package is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Agilent2Dicom source code.  If not, see
  <http://www.gnu.org/licenses/>.




## Quick summary ##

The *Agilent2Dicom* package is a series of bash and python scripts to
convert FDF images from the Agilent 9.4T MR scanner at Monash
Biomedical Imaging (MBI) into enhanced MR DICOM images.

Homepage: [MBI's Confluence homepage](https://confluence-vre.its.monash.edu.au/display/MBI/Agilent+FDF+to+Dicom+converter)

Source: [https://bitbucket.org/mbi-image/agilent2dicom](https://bitbucket.org/mbi-image/agilent2dicom)


The internal python script agilent2dicom.py is at version 0.6.

The main shell script fdf2dcm.sh and Agilent2Dicom is at version 1.1.


## Basic overview ##

The script fdf2dcm.sh is a wrapping bash script that cleans up the arguments and
calls a python script agilent2dicom to do the conversion of FDF files to basic
2D DICOM files, then calls dicom3tools' dcmulti to compress the basic MR format
to an enhanced MR format. It requires an input directory containing FDF images,
with optional arguments for the output directory that defaults to a subdrectory
called .dcm in the input directory. By default the output dicom file is named
0001.dcm and any mutiple dicom files are incremented from 1.


```
#!bash

[eagerm@m2101 Agilent2Dicom]$ ./fdf2dcm.sh -h
 usage: ./fdf2dcm.sh -i inputdir [-o outputdir] [-v] [-m] [-p]
 -i <inputdir> FDF source directory
 -o <outputdir> Optional destination DICOM directory. Default is input_dir/.dcm.
 -v             verbose output.
 -m,-p          Enable magnitude and phase subdirectory conversion. These flags are
passed to agilent2dicom and should only be used from within fdf2dcm or with knowledge
of input fdf data.
 -s <SEQ>            Sequence type (one of MULTIECHO,DIFFUSION,ASL).
 -x             Debug mode.
 -v		Verbose output
 -h             this help
```

The crux of the script is performed within the agilent2dicom python script. It's usage is listed below:

```
#!bash

[eagerm@m2101 Agilent2Dicom]$ ./agilent2dicom.py -h
usage: agilent2dicom -i "Input FDF directory" [-o "Output directory"] [-m] [-p] [-v] 
agilent2dicom is an FDF to Enhanced MR DICOM converter from MBI.
optional arguments:
  -h, --help            show this help message and exit
  -i INPUTDIR, --inputdir INPUTDIR
                        Input directory name. Must be an Agilent .img
                        directory containing procpar and fdf files
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Output directory name for DICOM files.
  -m, --magnitude       Magnitude component flag.
  -p, --phase           Phase component flag.
  -s SEQUENCE, --sequence SEQUENCE
                        Sequence type (one of Multiecho, Diffusion, ASL.
  -v, --verbose         Verbose.) 
```

A simple GUI interface *fdf2dicom*, was developed to be used by experimental scientists on the Agilent RedHat console. This requires python 2.6 and Tkinter.


## Current bugs and faults ##

For bugs and faults see [debugging page](https://confluence-vre.its.monash.edu.au/display/MBI/FDF2DCM+debugging) on Confluence.

## How do I get set up? ##

See INSTALL.txt

### Dependencies ###

 * python (2.6 or greater)
   - python-dicom
   - python-numpy
   - tkinter (for python 2.6 GUI fdf2dicom on Agilent console)
   - pyqt4 (for FDF2DicomQt and FDF2DicomApp)
 * dicom3tools
 * dcmtk
 * mrtrix (0.2.12 or mrtrix3)


## How to run tests ##

FDF examples are available on MASSIVE working repository, these
include standard 2d, multi-echo 3D, multiecho 2D complex, fast spin
echo 3D, diffusion (EPI and FSE), ASL, and CINE images. Most FDF
images are reconstructed from Vnmrj internally, but some can be
externally reconstructed and these are slightly different and need to
be tested.

The Makefile contains examples to be used on MASSIVE. These include
*standard2d*, *me3d*, *me2d*, *diffusion*, *asl*,
*cine* and other specific protocol examples. These will execute:

* *run_{}* runs the FDF agilent to dicom process and conversion of DICOMs to NIFTI;
* *check_{}*, check the dimensions of the Dicom and Nifti files; 
* *test_{}*, runs the dciodvfy on the first dicom to validate the conversion to enhanced MR;
* *view_{}*, view the DICOM image using mrview; and
* *viewn_{}*, view the NIFTI formatted image.
 
To make all the examples and test them use:

```
#!bash

make all
make test_all
```

### Python script testing ###

* FDF scripts

** ReadProcpar.py contains ReadProcpar method. ReadProcpar is used to read Agilent FDF config files

** ProcparToDicomMap.py contains the ProcparToDicomMap, getColumns and CreatUID methods. ProcparToDicomMap is used to map Agilent FDF header file (procpar) to DICOM dataset format.

Test the Procpar methods for multi-echo 3D:
```
#!bash

python ./ProcparToDicomMap.py -v -i  ~/Monash016/amanda/ExampleAgilentData/multiecho3d_magonly/
```

** ReadFDF.py contains the ReadFDF method.
** RescaleFDF.py contains the Rescale and FindScale methods.
** ParseFDF.py contains the ParseFDF, ParseDiffusionFDF and ParseASL methods.

Testing the ParseFDF methods for multi-echo 3D, gradient-echo 3D and diffusion examples:
```
#!bash

python ./ParseFDF.py -v -i  ~/Monash016/amanda/ExampleAgilentData/multiecho3d_magonly/
python ./ParseFDF.py -v -i ~/Monash016/RatKidney/Agilent/20120522/kidney512iso_01.fid
python ./ParseFDF.py -v -i ~/Monash016/amanda/ExampleAgilentData/diffusion/
```



* Complex filtering and reconstruction of FID data

cplxfilter.py includes the complex filtering methods from
scipy.ndimage.filters gaussian_filter, median_filter, gaussian_laplace
and the scipy.signal method wiener_filter.

ReadFID.py contains the readfid and recon methods for FID data conversion.

Python-scipy not supported in python/2.7.1-gcc on MASSIVE, but python-dicom not supported in python/2.7.3.
By including the site-package path of both versions allows testing on the cplxfilter methods.
```
#!bash

module unload python 
module load python/2.7.3-gcc
export PYTHONPATH=$PYTHONPATH:/usr/local/python/2.7.1-gcc/lib/python2.7/site-packages 
python2.7 ./cplxfilter.py -i ~/Monash016/RatKidney/Agilent/20120522/kidney512iso_01.fid
fslview raw_image.nii.gz new_image.nii.gz median_image.nii.gz
```


## Contribution guidelines ##

*Writing tests*

New example FDF image type testing procedures should have:

* a 'run_<type>' routine that converts the FDF folder to DICOM and converts the DICOMs to Nifti;
* a check_<type>' routine that displays the dimensions of the dicom and nifti files;
* a test_<type> routine that runs dciodvfy on the Dicom files;
* [optional] a view_ routine for mrview and fslview to show the DICOM or NIFTI images

See examples in Makefile.


*Code review*


*Other guidelines*

## Who do I talk to? ##

* Michael Eager (michael.eager@monash.edu) or someone in the Imaging Team at MBI

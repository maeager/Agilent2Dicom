# README #

## Quick summary ##

The *Agilent2Dicom* package is a series of bash and python scripts to
convert FDF and FID images from the Agilent 9.4T MR scanner at Monash
Biomedical Imaging (MBI) into enhanced MR DICOM images. 

Homepage: [MBI's Confluence homepage](https://confluence-vre.its.monash.edu.au/display/MBI/Agilent+FDF+to+Dicom+converter)

Source code: [https://bitbucket.org/mbi-image/agilent2dicom](https://bitbucket.org/mbi-image/agilent2dicom)

Agilent2Dicom is at version 1.3.0.


## Basic overview ##

The Agilent2Dicom package is for converting high-field MR images to
enhanced DICOM formats.  The idea stemmed from poor reconstruction and
implementation of the DICOM format on the Vnmrj system used by the
Agilent 9.4 MR scanner. Amanda Ng initiated the project with a simple
fdf converter. Michael Eager took over and expanded the features to
include enhanced DICOM conversion, full implementation for all
sequences, FID complex filtering, and pyqt4 GUIs.

FDF reconstructed images are converted to enhanced DICOM by the
FDF2DCM suite (*fdf2dcm.sh*, *agilent2dicom.py*, *ReadFDF.py* and
*ParseFDF.py*).  fdf2dcm.sh is a wrapping bash script that cleans up
the input arguments and calls a python script agilent2dicom to do the
conversion of FDF files to basic 2D DICOM files. fdf2dcm.sh then calls
dicom3tools' dcmulti to compress the basic MR format to an enhanced MR
format. It requires an input directory containing FDF images, with
optional arguments for the output directory that defaults to a
subdrectory called .dcm in the input directory. By default the output
dicom file is named 0001.dcm and any mutiple dicom files are
incremented from 1.

FID file formats contain the raw k-space data from the scanner.  The
FID2DCM suite (*fid2dcm.sh*, *fid2dicom.py*, *ReadFID.py*) enable the
reconstruction and convertions of images to DICOM.  A feature of the
reconstruction enables adjustable filtering in complex space and
storage of phase information.  *cplxfilter.py* allows filtering of
real and imaginary reconstructed images using gaussian, laplace-
gaussian, median and wiener filters.


## Usage ##

### FDF conversion to DICOM ##

Commandline usage of the FDF2DCM suite is as follows:
```
#!bash

[eagerm@m2101 Agilent2Dicom]$ ./fdf2dcm.sh -v -i  ../ExampleAgilentData/standard2d -o ../output_data/standard2d
```

The verbose argument (-v) gives more information and prompts the user
if they would like to delete any existing intermediate or final DICOM
files.

The crux of the script is performed within the agilent2dicom python
script. For conversion to standard 2D dicoms without combination
using DCMULTI, use agilent2dicom.py as follows:
```
#!bash

[eagerm@m2101 Agilent2Dicom]$ ./agilent2dicom.py -v -i ../ExampleAgilentData/standard2d -o ../output_data/standard2d
```


### FID filtering and conversion to DICOM ###

A script for processing raw FID data was developed and enables
filtering and bespoke reconstruction of MR data.  *fid2dicom* is the
main wrapper script that uses *ReadFID.py* and *cplxfilter.py* as its
core. The arguments to *fid2dicom* are similar to *fdf2dcm* in that it
requires an input directory, and an optional output directory.
*fid2dicom* has additional arguments for complex filtering including
Gaussian, Gaussian Laplace, Median, and Wiener filters.

Basic FID2DCM usage from the commandline is as follows:
```
#!bash

[eagerm@m2102 Agilent2Dicom]$ ./fid2dcm.sh -i inputdir.fid -o outputdir.dcm -v -m -p -k -g 1.0 
```

The enhanced DICOM files of the raw reconstruction will be stored in
"outputdir.dcm".  A Gaussian filter with sigma=1.0 is used on the real
and imaginary components of the reconstructed image. The magnitude and
phase of the Gaussian filtered image will be stored in
<outputdir>-gaussian-mag.dcm and <outputdir>-gaussian-pha.dcm,
respectively.  The -k argument saves the k-space data to a MATLAB mat
file (<outputdir>-ksp.mat).
 

## Advanced usage ##

### FDF ###
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

agilent2dicom.py has short and long arguments:
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


### FID  ###

Advanced FID2DCM usage from the commandline is as follows:
```
#!bash

[eagerm@m2102 Agilent2Dicom]$ ./fid2dcm.sh -h

 usage: ./fid2dcm.sh -i inputdir [-o outputdir] [-v] [-m] [-p] [-r] [-k] [[-g 1.0 -j 0 -e wrap] [-l 1.0] [-n 5] [ -w 5 -z 0.001]]
. To export components use magnitude (-m), phase (-p), real and imag (-r) or k-space (-k). Filtering is available for Gaussian filter (-g sigma), Laplace Gaussian filter (-l sigma), median filter (-n window_size), or Wiener filter (-w window_size).
 -i <inputdir>  FID source directory
 -o <outputdir> Optional destination DICOM directory. Default is input_dir/.dcm.
 -d             Disable dcmodify fixes to DICOMs.
 -m,-p,          Save magnitude and phase components.  These flags are passed to fid2dicom and should only be used from within fid2dcm or with knowledge of input fid data.
 -r             Save real and imaginary components of filtered image.
 -k             Save Kspace data.
 -g <sigma>     Gaussian filter smoothing of reconstructed RE and IM components. Sigma variable argument, default 1/srqt(2).
 -j <order>     Gaussian filter order variable argument, default 0.
 -e {'wrap','reflect','nearest','mirror'}  Gaussian filter mode variable argument, default=nearest.
 -l <simga>     Gaussian Laplace filter smoothing of reconstructed RE and IM components. Sigma variable argument, default 1/srqt(2).
 -n <wsize>     Median filter smoothing of reconstructed RE and IM components. Size variable argument, default 5.
 -w <wsize>     Wiener filter smoothing of reconstructed RE and IM components. Size variable argument, default 5.
 -z <noise>     Wiener filter noise variable, default 0 (none=default variance calculated in local region).
 -x             Debug mode.
 -v             Verbose output.
 -h             this help
```


*fid2dicom.py* has slighly different arguments to *fid2dcm.sh*
```
#!bash

[eagerm@m2108 Agilent2Dicom]$ python ./fid2dicom.py -h
usage:  fid2dicom.py -i "Input FID directory" [-o "Output directory"] [-m] [-p] [-v] [[-g -s 1.0 [-go 0 -gm wrap]] [-l -s 0.707] [-d -n 5] [-w -n 5]]

fid2dicom is an FID to Enhanced MR DICOM converter from MBI. Complex filtering
enabled with -g, -l, -n or -w arguments. FID2DICOM Version 1.3.0

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTDIR, --inputdir INPUTDIR
                        Input directory name. Must be an Agilent FID image
                        directory containing procpar and fid files
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Output directory name for DICOM files.
  -m, --magnitude       Save Magnitude component. Default output of filtered
                        image outputput
  -p, --phase           Save Phase component.
  -k, --kspace          Save Kspace data in outputdir-ksp.mat file
  -r, --realimag        Save real and imaginary data in outputdir-real and
                        outputdir-imag.
  -g, --gaussian_filter
                        Gaussian filter smoothing of reconstructed RE and IM
                        components.
  -l, --gaussian_laplace
                        Gaussian Laplace filter smoothing of reconstructed RE
                        and IM components. --sigma variable must be declared.
  -s SIGMA, --sigma SIGMA
                        Gaussian and Laplace-Gaussian sigma variable. Default
                        1/srqt(2).
  -go {0,1,2,3}, --gaussian_order {0,1,2,3}
                        Gaussian and Laplace-Gaussian order variable. Default
                        0.
  -gm {wrap,nearest,constant,reflect,mirror}, --gaussian_mode {wrap,nearest,constant,reflect,mirror}
                        Gaussian and Laplace-Gaussian mode variable. Default
                        nearest.
  -d, --median_filter   Median filter smoothing of reconstructed RE and IM
                        components.
  -w, --wiener_filter   Wiener filter smoothing of reconstructed RE and IM
                        components.
  -n WINDOW_SIZE, --window_size WINDOW_SIZE
                        Window size of Wiener and Median filters. Default 5.
  -wn WIENER_NOISE, --wiener_noise WIENER_NOISE
                        Wiener filter noise. Default None.
  -v, --verbose         Verbose.
```

# Agilent Console GUI #

A first attempt at a GUI interface, *fdf2dicom*, was developed to be
used by instrument scientists on the Agilent RedHat console for FDF
conversion only. This requires python 2.6 and Tkinter.
```
#!bash

[vnmr1@vnmrj1 s_2014062002]$ fdf2dicom
```

An updated GUI, *Agilent2DicomApp*, now uses PyQt4 and enables more
features including FID conversion.  With the Agilent2Dicom path
included in the PATH variable, the GUI can be loaded from the terminal
from any folder:
```
#!bash

[vnmr1@vnmrj1 s_2014062002]$ Agilent2DicomApp
```

On MASSIVE, enable the scipy and pyqt4 packages using the following commands:
```
#!bash

module load python/2.7.3-gcc pyqt4
export PYTHONPATH=/usr/local/pyqt4/4.11/lib/python2.7/site-packages/:/usr/local/python/2.7.3-gcc/lib/python2.7/site-packages:/usr/local/python/2.7.1-gcc/lib/python2.7/site-packages 
./Agilent2DicomApp.py
```



## Current bugs and faults ##

For bugs and faults see [debugging page](https://confluence-vre.its.monash.edu.au/display/MBI/FDF2DCM+debugging).

## How do I get set up? ##

See INSTALL.txt

### Dependencies ###

 * python (2.6 or greater)
   - python-dicom
   - python-numpy
   - tkinter (for python 2.6 GUI fdf2dicom on Agilent console)
   - pyqt4 (for Agilent2DicomQt and Agilent2DicomApp)
 * dicom3tools
 * dcmtk
 * mrtrix (0.2.12 or mrtrix3)


## How to run tests ##

FDF and FID examples are available on MASSIVE working repository
ExampleAgilentData, these include standard 2d, multi-echo 3D,
multiecho 2D complex, fast spin echo 3D, diffusion (EPI and FSE), ASL,
and CINE images. Most FDF images are reconstructed from Vnmrj
internally, but some can be externally reconstructed and these are
slightly different and need to be tested.

The Makefile contains examples to be used on MASSIVE. These include
*standard2d*, *me3d*, *me2d*, *diffusion*, *asl*, *cine* and other
specific protocol examples. These will execute:

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

* Procpar scripts

*ReadProcpar.py* contains ReadProcpar method. ReadProcpar is used to
 read Agilent FDF config files

*ProcparToDicomMap.py* contains the ProcparToDicomMap, getColumns and
 CreatUID methods. ProcparToDicomMap is used to map Agilent FDF header
 file (procpar) to DICOM dataset format.

Test the Procpar methods for multi-echo 3D:
```
#!bash

python ./ProcparToDicomMap.py -v -i  ../ExampleAgilentData/multiecho3d_magonly/
```

* FDF scripts

*ReadFDF.py* contains the ReadFDF method.
*RescaleFDF.py* contains the Rescale and FindScale methods.
*ParseFDF.py* contains the ParseFDF, ParseDiffusionFDF and ParseASL methods.

Testing the ParseFDF methods for multi-echo 3D, gradient-echo 3D and
diffusion examples:
```
#!bash

python ./ParseFDF.py -v -i ../ExampleAgilentData/multiecho3d_magonly/
python ./ParseFDF.py -v -i ../ExampleAgilentData//kidney512iso_01.fid
python ./ParseFDF.py -v -i ../ExampleAgilentData/diffusion/
```



### Testing the complex filtering and reconstruction of FID data feature ###

*ReadFID.py* contains the readfid and recon methods for FID data conversion.
*cplxfilter.py* includes the complex filtering methods from
scipy.ndimage.filters gaussian_filter, median_filter, gaussian_laplace
and the scipy.signal method wiener_filter.  Testing requires the
python-nibabel package for faster exporting and viewing of medical
images in NIFTI format.


Python-scipy is not supported in python/2.7.1-gcc on MASSIVE, but python-dicom not supported in python/2.7.3.
By including the site-package path of both versions allows testing on the cplxfilter methods.

*cplxfilter.py* will export simple nigti files for raw, gaussian, etc. given an input path:
```
#!bash

module unload python 
module load python/2.7.3-gcc
export PYTHONPATH=$PYTHONPATH:/usr/local/python/2.7.1-gcc/lib/python2.7/site-packages 
python2.7 ./cplxfilter.py -i ../ExampleAgilentData/kidney512iso_01.fid
fslview raw_image.nii.gz new_image.nii.gz median_image.nii.gz
```

*cplxfilter.py* also allows rearranging the image based on a three
 element argument (eg. -a 2,-0,-1).  The negative sign will reverse
 the slices along the axis. 


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

Mercurial keyword expansion is used to set some important variables. The hgrc file in the development .hg folder should be the following:

```
#!bash

[paths]
default = ssh://hg@bitbucket.org/mbi-image/agilent2dicom
[extensions]
keyword =
[keyword]
agilentFDF2dicom.py =
fdf2dcm.sh =
agilent2dicom_globalvars.py =
fid2dicom.py =
fid2dcm.sh =
Agilent2DicomAppQt.py =
[keywordmaps]
Author = {author|user}
Date = {date|utcdate}
Header = {root}/{file},v {node|short} {date|utcdate} {author|user}
Id = {file|basename},v {node|short} {date|utcdate} {author|user}
Revision = {node|short}
```

## Who do I talk to? ##

* Dr. Michael Eager (michael.eager@monash.edu) or someone in the Imaging Team at MBI


## Licence ##

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



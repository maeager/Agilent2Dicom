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


## What is this repository for? ##

### Quick summary ###

The *Agilent2Dicom* package is a series of bash and python scripts to
convert FDF images from the Agilent 9.4T MR scanner at Monash
Biomedical Imaging (MBI) into enhanced MR DICOM images.

Homepage: [MBI's Confluence homepage](https://confluence-vre.its.monash.edu.au/display/MBI/Agilent+FDF+to+Dicom+converter)

Source: [https://bitbucket.org/mbi-image/agilent2dicom](https://bitbucket.org/mbi-image/agilent2dicom)


### Version ###

The internal python script agilent2dicom.py is at version 0.5.

The main shell script fdf2dcm.sh and Agilent2Dicom is at version 1.1.

### Current bugs and faults ###

For bugs and faults see [debugging page](https://confluence-vre.its.monash.edu.au/display/MBI/FDF2DCM+debugging) on Confluence.

## How do I get set up? ##

### Summary of set up ###
### Configuration ###


Create a suitable directory (e.g. Agilent2DicomProject) in your source or projects folder.

```
#!bash

cd ~/src                   # on MASSIVE this is ~/Monash016/eagerm
mkdir Agilent2DicomProject
cd Agilent2DicomProject
```


### Dependencies ###

 * python (2.6 or greater)
   - python-dicom
   - python-numpy
   - tkinter (for python 2.6 GUI fdf2dicom on Agilent console)
 * dicom3tools
 * dcmtk
 * mrtrix (0.2.12 or mrtrix3)


### Setup and install ###

Pull the development code from MASSIVE using:

```
#!bash

hg clone username@m2.massive.org.au://gpfs/M2Home/projects/Monash016/eagerm/Agilent2Dicom/Agilent2Dicom
```

or a stable version from the MBI bitbucket repository, using:

```
#!bash

hg clone ssh://hg@bitbucket.org/mbi-image/agilent2dicom
```

An obsolete subversion repository is still available from the MeRC SVN service.
```
#!bash

svn checkout https://svn-vre.its.monash.edu.au/mbi/trunk/Sandbox/meager/Agilent2Dicom
```

### Install Dicom3tools ###

Install dicom3tools from http://www.dclunie.com/dicom3tools.html into the
Agilent2Dicom folder.  The default settings for MBI's Agilent 9.4T MR scanner
are used below.  The default UID root value '1.3.6.1.4.1' is unique to the MBI
Agilent 9.4T MR scanner and should be changed for other devices.  The dependencies for 
compiling dicom3tools are imake, gcc and binutil-essentials.

On MASSIVE use: 
```
#!bash

module load imake gcc.  
```

On the Agilent console, use: 
```
#!bash

sudo yum install imake gcc binutil-essentials
```

Once the dependencies have been installed, begin by downloading the
latest dicom3tools [bz2 tar ball](http://www.dclunie.com/dicom3tools/workinprogress/), adjust the
_config/site.p-def_ file and then follow the compile instructions
below:

```
#!bash

cd ~/src/Agilent2DicomProject    # On MASSIVE use cd ~/Monash016/eagerm/Agilent2Dicom
wget   http://www.dclunie.com/dicom3tools/workinprogress/dicom3tools_1.00.snapshot.20140306142442.tar.bz2
tar jxvf dicom3tools_1.00.snapshot.20140306142442.tar.bz2
cd dicom3tools_1.00.snapshot.20140306142442

sed -i 's/CLUNIE/MBIAGILENT/' config/site.p-def
./Configure
##            setenv IMAKEINCLUDE -I./config                              # only needed for tcsh
imake -I./config -DInstallInTopDir -DUseMBIAGILENT1ID -DDefaultUIDRoot=1.3.6.1.4.1
make World
make install                          # into ./bin
make install.man                      # into ./man
 
```


### How to run tests ###

FDF examples are available on MASSIVE working repository, these
include standard 2d, multi-echo 3D, multiecho 2D complex, fast spin
echo 3D, diffusion (EPI and FSE), ASL, and CINE images. Most FDF
images are reconstructed from Vnmrj internally, but some can be
externally reconstructed and these are slightly different and need to
be tested.

The Makefile contains examples to be used on MASSIVE. These include
*standard2d*, *me3d*, *me2d*, *diffusion*, *asl*,
*cine* and other specific protocol examples. These will execute:
       - *run_{}* runs the FDF agilent to dicom process and conversion of DICOMs to NIFTI;
       - *check_{}*, check the dimensions of the Dicom and Nifti files; and
       - *test_{}*, runs the dciodvfy on the first dicom to validate the conversion to enhanced MR.
 
To make all the examples and test them use:

```
#!bash

make all
make test_all
```

To view dicoms use *make view_<>*.

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
 -h             this help
```

The crux of the script is performed within the agilent2dicom python script. It's usage is listed below:

```
#!bash

[eagerm@m2101 Agilent2Dicom]$ ./agilent2dicom.py -h
usage: agilent2dicom -i "Input FDF directory" [-o Output directory] [-m] [-p] [-v]
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
 
```

A simple GUI interface *fdf2dicom*, was developed to be used by experimental scientists on the Agilent RedHat console. This requires python 2.6 and Tkinter.



## Contribution guidelines ##

* Writing tests
* Code review
* Other guidelines

## Who do I talk to? ##

* Michael Eager (michael.eager@monash.edu) or someone in the Imaging Team at MBI
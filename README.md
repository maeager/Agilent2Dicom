# README #

This small set of scripts aims to accurately convert FDF images from the Agilent 9.4T MR scanner at MBI into enhanced MR DICOM images.

[MBI intranet homepage](https://confluence-vre.its.monash.edu.au/display/MBI/Agilent+FDF+to+Dicom+converter)

### What is this repository for? ###

* Quick summary
* Version

Agilent2Dicom is at version 0.4

fdf2dcm is at 1.0

* Current bugs and faults

See [current issues](https://confluence-vre.its.monash.edu.au/display/MBI/FDF2DCM+debugging)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* How to run tests
* Deployment instructions

Create a suitable directory (Agilent2Dicom) in your source or projects folder (SOURCEDIR).

```
#!bash

cd ~/src                   # on MASSIVE this is ~/Monash016/eagerm
mkdir Agilent2Dicom
cd Agilent2Dicom
```

Pull the development code from MASSIVE using:

```
#!bash

hg clone username@m2.massive.org.au://gpfs/M2Home/projects/Monash016/eagerm/Agilent2Dicom/Agilent2Dicom
```

or a stable version from the MBI bitbucket repository, using:  (Note this is a stable repo without limited version history)

```
#!bash

hg clone ssh://hg@bitbucket.org/mbi-image/agilent2dicom

## Obsolete MeRC subversion reop
## svn checkout https://svn-vre.its.monash.edu.au/mbi/trunk/Sandbox/meager/Agilent2Dicom
```


Install dicom3tools from http://www.dclunie.com/dicom3tools.html into the Agilent2Dicom folder.  The default settings for MBI's Agilent 9.4T MR scanner are used below.  The default UID root value '1.3.6.1.4.1' is unique to the MBI Agilent 9.4T MR scanner and should be changed for other devices.
Dependencies: imake, gcc and other unix binutil essentials.  On MASSIVE use: module load imake gcc.  On the Agilent console, use: sudo yum install imake gcc binutil-essential

Install Dicom3tools


```
#!bash

cd ~/src/Agilent2Dicom
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

### Basic overview ###

The script fdf2dcm.sh is a wrapping bash script that cleans up the arguments and calls a python script agilent2dicom to do the conversion
of FDF files to basic 2D DICOM files, then calls dicom3tools' dcmulti to compress the basic MR format to an enhanced MR format. It requires
an input directory containing FDF images, with optional arguments for the output directory that defaults to a subdrectory called .dcm in the input directory. By default the
output dicom file is named 0001.dcm and any mutiple dicom files are incremented from 1.


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

[eagerm@m2101 Agilent2Dicom]$ ./agilent2dicom -h
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

A simple GUI interface *fdf2dicom.sh*, was developed to be used by experimental scientists on the Agilent RedHat console. 



### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Michael Eager (michael.eager@monash.edu) or someone in the Imaging Team at MBI
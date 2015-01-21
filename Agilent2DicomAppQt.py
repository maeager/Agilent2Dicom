#!/usr/bin/env python
#
# $Header$
# $Id$
#
# Version 1.2.5: Working version on Redhat Workstation
# Version 1.3.0: Info tab panels show information from Procpar

#
# Copyright 2014 Michael Eager
#
# This file is part of the Agilent2Dicom package
# See https://bitbucket.org/mbi-image/agilent2dicom/ for documentation.
#
# Agilent2dicom is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Agilent2dicom is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Agilent2dicom. If not, see <http://www.gnu.org/licenses/>.

"""
GUI for FDF/FID Agilent to Dicom converter. Extract Ui_Form.py using:
pyuic4 --output Agilent2DicomQt.py Agilent2DicomQt.ui
""" 


import os
import subprocess
import sys
import re
from PyQt4 import Qt, QtGui, QtCore
from PyQt4.QtGui import QDialog,QFileDialog,QApplication
from Agilent2DicomQt import Ui_MainWindow
import ReadProcpar
from agilent2dicom_globalvars import *
DEBUGGING=1

#Agilent2DicomAppVersion=0.7
__author__ = "Michael Eager, Monash Biomedical Imaging"
__version__ = str(Agilent2DicomAppVersion)
__date__ = "$Date$"
__copyright__ = "Copyright 2014 Michael Eager"


Agilent2DicomAppStamp=re.sub(r'\$Id(.*)\$',r'\1',"$Id$")
cmd_header='(if test ${MASSIVE_USERNAME+defined} \n\
then \n\
echo ''On Massive'' \n\
module purge \n\
module load massive virtualgl\n\
module load python/2.7.1-gcc \n\
module load python/2.7.3-gcc \n\
module load dcmtk mrtrix \n\
#module list \n\
#export PYTHONPATH=/usr/local/python/2.7.3-gcc/lib/python2.7/site-packages:/usr/local/pyqt4/4.11/lib/python2.7/site-packages:/usr/local/python/2.7.1-gcc/lib/python2.7:/usr/local/python/2.7.1-gcc/lib/python2.7/site-packages \n\
else \n\
echo ''Not in MASSIVE'' \n\
fi \n\
echo ''Done''\n '

mrview_header='if test ${MASSIVE_USERNAME+defined} \n\
then \n\
echo ''On Massive'' \n\
module purge \n\
module load massive virtualgl\n\
module load mrtrix \n\
module list \n\
GL=vglrun\n\
else GL= \n\
fi; $GL mrview '


class Agilent2DicomWindow(QtGui.QMainWindow):

    NIFTI=0 # save to nifti flag

    def __init__(self):
        super(Agilent2DicomWindow,self).__init__()
        self.ui=Ui_MainWindow()
    # Set up the user interface from Designer.
        self.ui.setupUi(self)
        
        # Make some local modifications.
        # self.colorDepthCombo.addItem("2 colors (1 bit per pixel)")

        
        # Disable some features
        if DEBUGGING == 0:
#            self.ui.tab_diffusion.setEnabled(False)
#            self.ui.tab_multiecho.setEnabled(False)
            self.ui.pushButton_check.setEnabled(False)
            self.ui.pushButton_view.setEnabled(False)
            self.ui.pushButton_send2daris.setEnabled(False)
            self.ui.pushButton_check2.setEnabled(False)
            self.ui.pushButton_view2.setEnabled(False)
            self.ui.pushButton_send2daris2.setEnabled(False)
            self.ui.checkBox_median.setEnabled(False)
            self.ui.lineEdit_median_size.setEnabled(False)
            self.ui.checkBox_wiener.setEnabled(False)
            self.ui.lineEdit_wiener_size.setEnabled(False)
            self.ui.lineEdit_wiener_noise.setEnabled(False)
            self.ui.checkBox_magn.setEnabled(True)
                                             
        
        self.ui.checkBox_magn.setChecked(True)
        self.ui.checkBox_ksp.setChecked(False)
        self.ui.checkBox_reimag.setChecked(False)
        self.ui.checkBox_pha.setChecked(False)
        self.ui.checkBox_nodcmulti.setChecked(False)
        self.ui.checkBox_debugging.setChecked(False)
        self.ui.checkBox_nifti.setChecked(False)
        self.ui.checkBox_gaussian3D.setChecked(False)
        self.ui.checkBox_gaussian2D.setChecked(False)
        self.ui.checkBox_median.setChecked(False)
        self.ui.checkBox_wiener.setChecked(False)
        self.ui.checkBox_epanechnikov3D.setChecked(False)
        self.ui.checkBox_epanechnikov2D.setChecked(False)
        self.ui.checkBox_kspgaussian.setChecked(False)
        self.ui.checkBox_kspgaussshift.setChecked(False)
        self.ui.checkBox_kspgauss_super.setChecked(False)
        self.ui.checkBox_kspepa.setChecked(False)
        self.ui.checkBox_kspepashift.setChecked(False)
        self.ui.checkBox_kspepa_super.setChecked(False)
        # Connect up the buttons.
        #self.connect(self.ui.buttonBox, Qt.SIGNAL("accepted()"), self.accept)
        #self.connect(self.ui.buttonBox ,Qt.SIGNAL("rejected()"), self.reject)
 
        self.ui.pushButton_changefdf.clicked.connect(self.ChangeFDFpath)
        self.ui.pushButton_changedicom.clicked.connect(self.ChangeFDFDicomPath)
        self.ui.pushButton_changefid.clicked.connect(self.ChangeFIDpath)
        self.ui.pushButton_changedicom2.clicked.connect(self.ChangeFIDDicomPath)
        self.ui.pushButton_convert.clicked.connect(self.ConvertFDF)
        self.ui.pushButton_check.clicked.connect(self.CheckFDF)
        self.ui.pushButton_view.clicked.connect(self.ViewFDF)
        self.ui.pushButton_send2daris.clicked.connect(self.Send2Daris)
        self.ui.pushButton_convertfid.clicked.connect(self.ConvertFID)
        self.ui.pushButton_check2.clicked.connect(self.CheckFID)
        self.ui.pushButton_view2.clicked.connect(self.ViewFID)
        self.ui.pushButton_send2daris2.clicked.connect(self.Send2DarisFID)

        QtCore.QObject.connect(self.ui.actionSave_Filter_Outputs_to_Nifti, QtCore.SIGNAL(QtCore.QString.fromUtf8("triggered()")), self.toggleNifti)

        QtCore.QObject.connect(self.ui.actionAbout, QtCore.SIGNAL(QtCore.QString.fromUtf8("triggered()")), self.About)
        QtCore.QObject.connect(self.ui.actionHelp, QtCore.SIGNAL(QtCore.QString.fromUtf8("triggered()")), self.About)
#        QtCore.QMetaObject.connectSlotsByName(MainWindow)


        
    def accept(self):
        '''Execute the command in response to the OK button.'''
        self.close()

    def reject(self):
        '''Cancel.'''
        self.close() 

    def toggleNifti(self):
        self.NIFTI=(self.NIFTI + 1) % 2
        
    # @QtCore.pyqtSlot()
    def ChangeFDFpath(self):
        success=1
        try:
            newdir = str(QFileDialog.getExistingDirectory(self, "Select FDF Directory"))
            print newdir
            
            self.ui.lineEdit_fdfpath.setText(newdir)
            files = os.listdir(newdir)
            if 'procpar' not in files:
                print 'ChangeFDFPath Error: FDF folder does not contain a procpar file'
                success=0
            
            fdffiles = [ f for f in files if f.endswith('.fdf') ]
            if len(fdffiles) == 0:
                print 'Error: FDF folder does not contain any fdf files'
                success=0
            if success:
                if re.search('img',newdir):
                    out = re.sub('img','dcm',newdir)
                else:
                    out = newdir+'.dcm'

                self.ui.lineEdit_dicompath.setText(out)
                self.ui.lineEdit_darisid.setText(self.GetDarisID(newdir))
                
            self.ui.checkBox_gaussian3D.setChecked(False)
            self.ui.checkBox_gaussian2D.setChecked(False)
            self.ui.checkBox_median.setChecked(False)
            self.ui.checkBox_wiener.setChecked(False)
            self.UpdateGUI()
        except ValueError:
            pass
	      
    def ChangeFDFDicomPath(self):
        try:
            newdir = str(QFileDialog.getExistingDirectory(self, "Select DICOM Directory"))
            
            files = os.listdir(newdir)            
            dcmfiles = [ f for f in files if f.endswith('.dcm') ]
            if len(dcmfiles) == 0:
                print 'FDF output DICOM folder does not contain any dcm files'
            else:
                quit_msg = "Are you sure you want to delete existing dicoms?"
                reply = QtGui.QMessageBox.question(self, 'Message', 
                     quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

                if reply == QtGui.QMessageBox.Yes:
#                    event.accept()
                    self.ui.lineEdit_dicompath.setText(newdir)
                    self.ui.lineEdit_darisid.setText(self.GetDarisID(str(self.ui.lineEdit_fdfpath.text())))
#                else:
#                    event.ignore()
        
            self.UpdateGUI()
        except ValueError:
            pass
		  
    def ChangeFIDpath(self):
        success=1
        try:
            newdir = str(QFileDialog.getExistingDirectory(self, "Select FID Directory"))
            self.ui.lineEdit_fidpath.setText(newdir)
            
            files = os.listdir(newdir)
            if 'procpar' not in files:
                print 'ChangeFIDPath Error: FID folder does not contain a procpar file'
                success=0
            
            fidfiles = [ f for f in files if f.endswith('fid') ]
            if len(fidfiles) == 0:
                print 'Error: FID folder does not contain any fid files'
                success=0
            if success:
                if re.search('fid',newdir):
                    out = re.sub('fid','dcm',newdir)
                else:
                    out = newdir+'.dcm'
                self.ui.lineEdit_dicompath2.setText(out)
                self.ui.lineEdit_darisid2.setText(self.GetDarisID(newdir))
                
            self.ui.checkBox_gaussian3D.setChecked(False)
            self.ui.checkBox_gaussian2D.setChecked(False)
            self.ui.checkBox_median.setChecked(False)
            self.ui.checkBox_wiener.setChecked(False)
            self.UpdateGUI()

        except ValueError:
            pass
			  
    def ChangeFIDDicomPath(self):
        try:
            newdir = str(QFileDialog.getExistingDirectory(self, "Select or create new DICOM Directory"))
            files = os.listdir(newdir)            
            dcmfiles = [ f for f in files if f.endswith('.dcm') ]
            if len(dcmfiles) == 0:
                print 'FDF output DICOM folder does not contain any dcm files'
            else:
                quit_msg = "Are you sure you want to delete existing dicoms?"
                reply = QtGui.QMessageBox.question(self, 'Message', 
                     quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

                if reply == QtGui.QMessageBox.Yes:
                    self.ui.lineEdit_dicompath2.setText(newdir)
                    self.ui.lineEdit_darisid.setText(self.GetDarisID(str(self.ui.lineEdit_fidpath.text())))
            self.UpdateGUI()
        except ValueError:
            pass
			      
    def ConvertFDF(self):
        try:

            input_dir = self.ui.lineEdit_fdfpath.text()
            output_dir = self.ui.lineEdit_dicompath.text()
            thispath = str(os.path.dirname(os.path.realpath(os.path.abspath(__file__))))
            print 'fdf2dcm path: %s' % thispath
            cmd1 = str(os.path.join(thispath,'fdf2dcm.sh')) + ' -i ' + str(input_dir) + ' -o ' + str(output_dir)
            if self.ui.checkBox_nodcmulti.isChecked():
                cmd1+=' -d '
            if self.ui.checkBox_debugging.isChecked():
                cmd1+=' -v '
                
            print(cmd1)
            cmd= cmd_header + cmd1 +')'
            print subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
            self.UpdateGUI()
        except ValueError:
            pass
				  
    def GetDarisID(self,inpath):
        daris_id=''
        try:
            print "GetDarisID %s " % inpath
            procpar,pep = ReadProcpar.ReadProcpar(os.path.join(os.path.abspath(str(inpath)),'procpar'))
            if 'name' in procpar.keys():
                if re.search('DaRIS',procpar['name']):
                    daris_id = re.sub('DaRIS\^','',procpar['name'])
        except ValueError:
            pass
        return daris_id

    def CheckDicomDir(self,dpath):
        if dpath != "" and os.path.isdir(dpath):
            # Check folder contains procpar and *.fdf files
            files = os.listdir(dpath)
        
            dcmfiles = [ f for f in files if f.endswith('.dcm') ]
            if len(dcmfiles) == 0:
                print 'Error: DICOM folder does not contain any dcm files'
                return False
            else:
                return True
        else:
#            print 'Error: DICOM folder does not exist: %s ' % dpath
            return False
					  					  
    def Send2Daris(self):
        try:
            daris_ID = self.ui.lineEdit_darisid.text()
            if str(daris_ID)=="": 
                print "No DaRIS ID"
                QtGui.QMessageBox.warning(self, 'Warning',"Cannot send to DaRIS without a proper ID.")
                return
            dicom_dir = str(self.ui.lineEdit_dicompath.text())
            if dicom_dir=="" or not os.path.isdir(dicom_dir): 
                print "No Dicom directory"
                QtGui.QMessageBox.warning(self, 'Warning',"Cannot send to DaRIS without a proper Dicom path.")
                return
                         
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd = os.path.join(thispath,'dpush') + ' -c ' + str(daris_ID) + ' -s mf-erc ' + str(dicom_dir)
            send_msg = "Are you sure you want to send to DaRIS the dicom directory\n"+ \
                str(dicom_dir)+"\n using the ID: "+str(daris_ID)+"   ?"
            reply = QtGui.QMessageBox.question(self, 'Message', 
                     send_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.Yes:
                print subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
            self.UpdateGUI()
        except ValueError:
            pass
					      
    def CheckFDF(self):  #send_button):
        try:
            output_dir = str(self.ui.lineEdit_dicompath.text())
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd1 ='mrinfo '+ output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'

            print subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
            cmd1 = os.path.join(thispath,'dcheck.sh') + ' -o ' + output_dir

            print subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
            self.UpdateGUI()
        except ValueError:
            pass
						  
    def ViewFDF(self):
        try:
            output_dir = str(self.ui.lineEdit_dicompath.text())
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd1 =mrview_header + output_dir 
            #print(cmd1)
            print subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
        except ValueError:
            pass

    def comboSigmaScale(self,text):
        if text == "unit voxel":
            self.sigmafactor=1
        else:
            from ReadProcpar import ReadProcpar, ProcparInfo            
            procpar, procpartext = ReadProcpar(os.path.join(str(self.ui.lineEdit_fdfpath.text()),'procpar'))
            p = ProcparInfo(procpar)
            if text == "in mm":
                self.sigmafactor=p['Voxel_Res_mm']
            else:
                self.sigmafactor=p['Voxel_Res_mm']*1000
            
    def SigmaString(self):
        s=str(self.ui.lineEdit_gsigma.text())
        if len(self.sigmafactor)!=3:
            print "Sigma units needs to be three elements"
        try: #if s.find(','):
            x, y, z = map(float, s.split(','))
            argstr=' -g %f,%f,%f' % (x/self.sigmafactor[0],y/self.sigmafactor[1],z/self.sigmafactor[2])
        except ValueError:
            argstr=' -g %f,%f,%f' % (float(s)/self.sigmafactor[0],float(s)/self.sigmafactor[1],float(s)/self.sigmafactor[2])            
        return argstr
            
    def getCommandArgs(self):
        argstr=''
        if not (self.ui.checkBox_magn.isChecked() and self.ui.checkBox_pha.isChecked()
                and  self.ui.checkBox_ksp.isChecked() and self.ui.checkBox_reimag.isChecked()):
            print "Forcing magnitude export since no export checkboxes enabled."
            argstr=' -m'
        if self.ui.checkBox_magn.isChecked():
            argstr=' -m'
        if self.ui.checkBox_pha.isChecked():
            argstr+=' -p'
        if self.ui.checkBox_ksp.isChecked():
            argstr+=' -k'
        if self.ui.checkBox_reimag.isChecked():
            argstr+=' -r'
        if self.ui.checkBox_doubleresolution.isChecked():
            argstr+=' -D'
        
        if self.ui.checkBox_gaussian2D.isChecked():
            argstr+=' -2 -g ' % self.SigmaString()            
            argstr+=' -j %s' % str(self.ui.lineEdit_gorder.text())
            if self.ui.nearest.isChecked():
                argstr+=' -e nearest'
            elif self.ui.reflect.isChecked():
                argstr+=' -e reflect'
            elif self.ui.wrap.isChecked:
                argstr+=' -e wrap'
            else:
                argstr+=' -e nearest'
        if self.ui.checkBox_gaussian3D.isChecked():
            argstr+=' -g ' % self.SigmaString()
            argstr+=' -j %s' % str(self.ui.lineEdit_gorder.text())
            if self.ui.nearest.isChecked():
                argstr+=' -e nearest'
            elif self.ui.reflect.isChecked():
                argstr+=' -e reflect'
            elif self.ui.wrap.isChecked:
                argstr+=' -e wrap'
            else:
                argstr+=' -e nearest'
        if self.ui.checkBox_fgaussian.isChecked():
            argstr+=' -G ' % self.SigmaString()
            argstr+=' -j %s' % str(self.ui.lineEdit_gorder.text())
            if self.ui.nearest.isChecked():
                argstr+=' -e nearest'
            elif self.ui.reflect.isChecked():
                argstr+=' -e reflect'
            elif self.ui.wrap.isChecked:
                argstr+=' -e wrap'
            else:
                argstr+=' -e nearest'

        if self.ui.checkBox_median.isChecked():
            argstr+=' -n %s ' % (str(self.ui.lineEdit_median_size.text()))
        if self.ui.checkBox_wiener.isChecked():
            argstr+=' -w %s -z %s' % (str(self.ui.lineEdit_wiener_size.text()), str(self.ui.lineEdit_wiener_noise.text()))
        if self.NIFTI == 1:
            argstr+=' -f'
        return argstr
									    
    def ConvertFID(self):
        try:
            input_dir = self.ui.lineEdit_fidpath.text()
            output_dir = self.ui.lineEdit_dicompath2.text()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'fid2dcm path: %s' % thispath
            cmd1 = os.path.join(thispath,'fid2dcm.sh')
            cmd1 = cmd1+' -v '+ self.getCommandArgs()
            cmd1 = cmd1+' -i ' + str(input_dir) + ' -o ' + str(output_dir)
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            #print(cmd)
            print subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
            self.UpdateGUI()
        except ValueError:
            pass

    def CheckFID(self):  #send_button):
        try:            
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'check path: %s' % thispath
            dicom_raw_dir = str(self.ui.lineEdit_dicompath2.text())
            # Shorten string if it ends in '/'
            if dicom_raw_dir[-1] == '/': dicom_raw_dir=dicom_raw_dir[:-1]
            # Get tuple of root path and raw/unfiltered dicom path
            (dicom_dir_root, dicom_raw) = os.path.split(dicom_raw_dir)
            # Do a regex and get all the dicom paths produced by Agilent2Dicom
            rgx = re.compile(r''+re.sub('.dcm','',dicom_raw)+".*.dcm")
            for dicom_dir in filter(rgx.match,os.listdir(dicom_dir_root)):
                # os.system("ls "+dicom_dir_root+" | grep '"+re.sub('.dcm','',dicom_raw)+".*.dcm'"):
                dcmpath=os.path.join(dicom_dir_root,dicom_dir)
                if not os.path.isdir(dcmpath) or len(os.listdir(dcmpath))<=2:
                    cmd1 = os.path.join(thispath,'dcheck.sh') + ' -o ' + str(dcmpath)
                    #print(cmd1)
                    #cmd=cmd_header + cmd1 +')'
                    #print(cmd)
                    print subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
            self.UpdateGUI()
        except ValueError:
            pass
										    
    def ViewFID(self):
        try:
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            dicom_raw_dir = str(self.ui.lineEdit_dicompath2.text())
            # Shorten string if it ends in '/'
            if dicom_raw_dir[-1] == '/': dicom_raw_dir=dicom_raw_dir[:-1]
            # Get tuple of root path and raw/unfiltered dicom path
            (dicom_dir_root, dicom_raw) = os.path.split(dicom_raw_dir)
            # Do a regex and get all the dicom paths produced by Agilent2Dicom
            rgx = re.compile(r''+re.sub('.dcm','',dicom_raw)+".*.dcm")
            for dicom_dir in filter(rgx.match,os.listdir(dicom_dir_root)):
                #os.system("ls "+dicom_dir_root+" | grep '"+re.sub('.dcm','',dicom_raw)+".*.dcm'"):
                dcmpath=os.path.join(dicom_dir_root,dicom_dir)
                if not os.path.isdir(dcmpath) or len(os.listdir(dcmpath))<=2:
                    cmd1 =mrview_header + dcmpath +' & '
                    
                    print(cmd1)
                    print subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
                    
            self.UpdateGUI()
        except ValueError:
            pass
											
    def Send2DarisFID(self):
        try:
            daris_ID = self.ui.lineEdit_darisid2.text()
            if str(daris_ID)=="": 
                print "No DaRIS ID"
                QtGui.QMessageBox.warning(self, 'Warning',
                     "Cannot send to DaRIS without a proper ID.")
                return
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            
            dicom_raw_dir = str(self.ui.lineEdit_dicompath2.text())
            # Shorten string if it ends in '/'
            if dicom_raw_dir[-1] == '/': dicom_raw_dir=dicom_raw_dir[:-1]
            # Get tuple of root path and raw/unfiltered dicom path
            (dicom_dir_root, dicom_raw) = os.path.split(dicom_raw_dir)
            # Do a regex and get all the dicom paths produced by Agilent2Dicom
            rgx = re.compile(r''+re.sub('.dcm','',dicom_raw)+".*.dcm")
            for dicom_dir in filter(rgx.match,os.listdir(dicom_dir_root)):
                #os.system("ls "+dicom_dir_root+" | grep '"+re.sub('.dcm','',dicom_raw)+".*.dcm'"):
                dcmpath=os.path.join(dicom_dir_root,dicom_dir)
                if not os.path.isdir(dcmpath) or len(os.listdir(dcmpath))<=2:
                    QtGui.QMessageBox.warning(self, 'Warning',"Cannot send to DaRIS. Directory "+dicom_dir+" is empty")
                else:
                    cmd = os.path.join(thispath,'dpush') + ' -c ' + str(daris_ID) + ' -s mf-erc ' + str(dcmpath)
                    send_msg = "Are you sure you want to send to DaRIS the dicom directory\n"+ \
                        str(dcmpath)+"\n using the ID: "+str(daris_ID)+"   ?"
                    reply = QtGui.QMessageBox.question(self, 'Message', send_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
                    if reply == QtGui.QMessageBox.Yes:
                        print subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
                        
        except ValueError:
            pass
											    
    def UpdateGUI(self):
        self.ui.pushButton_check.setEnabled(False)
        self.ui.pushButton_view.setEnabled(False)
        self.ui.pushButton_send2daris.setEnabled(False)
        self.ui.pushButton_check2.setEnabled(False)
        self.ui.pushButton_view2.setEnabled(False)
        self.ui.pushButton_send2daris2.setEnabled(False)
        if self.CheckDicomDir(self.ui.lineEdit_dicompath.text()):
            self.ui.pushButton_check.setEnabled(True)
            self.ui.pushButton_view.setEnabled(True)
            self.ui.pushButton_send2daris.setEnabled(True)
        if self.CheckDicomDir(self.ui.lineEdit_dicompath2.text()):
            self.ui.pushButton_check2.setEnabled(True)
            self.ui.pushButton_view2.setEnabled(True)
            self.ui.pushButton_send2daris2.setEnabled(True)
        if self.ui.lineEdit_fdfpath.text():
            self.FDFInfoUpdate()
        if self.ui.lineEdit_fidpath.text():
            self.FIDInfoUpdate()
            
												  
    def FDFInfoUpdate(self):
        from ReadProcpar import ReadProcpar, ProcparInfo            
        procpar, procpartext = ReadProcpar(os.path.join(str(self.ui.lineEdit_fdfpath.text()),'procpar'))
        SequenceText="FDF image:"
        
        p = ProcparInfo(procpar)
        self.ui.FDFprocparInfo.setText(SequenceText+'\n'.join("%s:\t\t%r" % (key,val) for (key,val) in p.iteritems()))


        
    def FIDInfoUpdate(self):
        from ReadProcpar import ReadProcpar, ProcparInfo
        procpar, procpartext = ReadProcpar(os.path.join(str(self.ui.lineEdit_fidpath.text()),'procpar'))
        #ds,MRAcq_type = ProcparToDicomMap(procpar)
        try:
            ftext= open(os.path.join(str(self.ui.lineEdit_fidpath.text()),'text'),'r')
            SequenceText="FID image:\n%s" % ftext.readline()
            ftext.close()
        except ValueError:
            SequenceText="FID image:\n-no text-"
            pass
        # DisplayText=[ SequenceText+"\n",    
        #  "StudyID: %s\n" %  ds.StudyID,
        #  "Patient ID: %s\n" % ds.PatientID,
        #  "Protocol name: %s\n"%  ds.ProtocolName,
        #  "Series Desc: %s\n" % ds.SeriesDescription,
        #  "Image type: %s\n" % ds.ImageType,
        #  "MR Acquisition Type: %s\n" % MRAcq_type,
        #  "Acq Matrix: %s, %s\n" % (ds.Rows, ds.Columns),
        #  "Pixel Spacing: %s, %s\n" % (ds.PixelSpacing[0], ds.PixelSpacing[1]),
        #  "Slice thickness: %s\n" % ds.SliceThickness]
        # #for column, key in enumerate(self.data):
        # #    for row, item in enumerate(self.data[key]):
        # #        newitem = QTableWidgetItem(item)
        # #        self.ui.listViewsetItem(row, column, newitem)
        # print DisplayText
        # QtGui.QMessageBox.question(self, 'Info','\n'.join(DisplayText),QtGui.QMessageBox.Ok)

        p = ProcparInfo(procpar)
        self.ui.FIDprocparInfo.setText(SequenceText+'\n'.join("%s:\t\t%r" % (key,val) for (key,val) in p.iteritems()))
        
    def About(self):
        msg = "                      Agilent2Dicom       \n\n"+\
            "Agilent2Dicom converts FDF and FID images from the Agilent 9.4T MR scanner at Monash Biomedical Imaging (MBI) into enhanced MR DICOM images.\n\n"+\
            "Homepage: https://confluence-vre.its.monash.edu.au/display/MBI/Agilent+FDF+to+Dicom+converter \n\n"+\
            "Version: "+__version__+"\n"+\
            "Stamp: "+Agilent2DicomAppStamp+"\n"+\
            "Author: "+__author__+"\n"+\
            __copyright__
        reply = QtGui.QMessageBox.question(self, 'Message', msg, QtGui.QMessageBox.Ok)



if __name__ == "__main__":
    app = QApplication(sys.argv)
    main = Agilent2DicomWindow()
    main.show()
    sys.exit(app.exec_())
												  

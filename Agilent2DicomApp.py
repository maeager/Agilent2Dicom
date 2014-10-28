#!/usr/bin/env python
#
# $Header: $
# $Id: $
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
pyuic4 --output Agilent2DicomWidget.py Agilent2DicomWidget.ui
""" 


import os
import sys
import re
from PyQt4 import Qt, QtGui, QtCore
from PyQt4.QtGui import QDialog,QFileDialog,QApplication
from Agilent2DicomWidget import Ui_Form
import ReadProcpar
from agilent2dicom_globalvars import *
DEBUGGING=0


__author__ = "Michael Eager, Monash Biomedical Imaging"
__version__ = str(Agilent2DicomAppVersion)+"-$Revision$"
__date__ = "$Date$"
__copyright__ = "Copyright 2014 Michael Eager"


Agilent2DicomAppStamp="$Id:"
cmd_header='(if test ${MASSIVE_USERNAME+defined} \n\
then \n\
echo ''On Massive'' \n\
module unload python \n\
module load python/2.7.1-gcc \n\
module load python/2.7.3-gcc \n\
module load dcmtk mrtrix \n\
module list \n\
export PYTHONPATH=/usr/local/python/2.7.3-gcc/lib/python2.7/site-packages:/usr/local/pyqt4/4.11/lib/python2.7/site-packages:/usr/local/python/2.7.1-gcc/lib/python2.7:/usr/local/python/2.7.1-gcc/lib/python2.7/site-packages \n\
else \n\
echo ''Not in MASSIVE'' \n\
fi \n\
echo ''Done''\n '

class Agilent2DicomWindow(QtGui.QWidget):
    def __init__(self):
        super(Agilent2DicomWindow,self).__init__()
        self.ui=Ui_Form()
    # Set up the user interface from Designer.
        self.ui.setupUi(self)
        self.ui.setWindowTitle(_translate("Form", "MBI\'s Agilent to Dicom converter application ("+__version__+")", None))
                
        # Make some local modifications.
        # self.colorDepthCombo.addItem("2 colors (1 bit per pixel)")
        
        # Disable some features
        if DEBUGGING == 0:
            self.ui.tab_diffusion.setEnabled(False)
            self.ui.tab_multiecho.setEnabled(False)
            self.ui.pushButton_check.setEnabled(False)
            self.ui.pushButton_view.setEnabled(False)
            self.ui.pushButton_send2daris.setEnabled(False)
            self.ui.pushButton_check2.setEnabled(False)
            self.ui.pushButton_view2.setEnabled(False)
            self.ui.pushButton_send2daris2.setEnabled(False)
            self.ui.checkBox_median.setChecked(False)
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
      
      
        # Connect up the buttons.
        #self.connect(self.ui.buttonBox, Qt.SIGNAL("accepted()"), self.accept)
        #self.connect(self.ui.buttonBox ,Qt.SIGNAL("rejected()"), self.reject)
        #self.ui.buttonBox.Ok.clicked.connect(self.accept)
        #self.ui.buttonBox.cancelButton.clicked.connect(self.reject)
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
        self.ui.pushButton_send2daris2.clicked.connect(self.Send2Daris2)

    def accept(self):
        '''Execute the command in response to the OK button.'''
        #print 'The volume of drinks is {0} liters ({1} jiggers).'.format(self.getLiters(), self.getJiggers())
        #print 'The blender is running at speed "{0}"'.format(self.getSpeedName())
        self.close()

    def reject(self):
        '''Cancel.'''
        self.close() 
      
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
                
            self.ui.checkBox_gaussian.setChecked(False)
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
                
            self.ui.checkBox_gaussian.setChecked(False)
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
            print(cmd)
            os.system(cmd)
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
            dicom_dir = self.ui.lineEdit_dicompath.text()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd = os.path.join(thispath,'dpush') + ' -c ' + str(daris_ID) + ' -s mf-erc ' + str(dicom_dir)
            tmp_msg = "Sending DICOM directory "+str(dicom_dir)+" to DaRIS.\n"+cmd+"\n\nAre you sure?"
            reply = QtGui.QMessageBox.question(self, 'Message', 
                     tmp_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.Yes:
                os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
					      
    def CheckFDF(self):  #send_button):
        try:
            output_dir = str(self.ui.lineEdit_dicompath.text())
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'check path: %s' % thispath
            cmd1 = os.path.join(thispath,'dcheck.sh') + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            cmd1 ='mrinfo '+ output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
						  
    def ViewFDF(self):
        try:
            output_dir = str(self.ui.lineEdit_dicompath.text())
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            if os.path.isdir(os.path.join(output_dir,'tmp')):
                tmp_msg = "Tmp folder found in DICOM directory. You cannot view the dicoms if there are subdirectories with additional dicoms. Are you sure you want to delete the temporary 2D dicom files?"
                reply = QtGui.QMessageBox.question(self, 'Message', 
                     tmp_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

                if reply == QtGui.QMessageBox.Yes:
                    import shutil
                    shutil.rmtree(os.path.join(output_dir,'tmp')) 
                else:
                    return
            
            cmd1 ='vglrun mrview '+ output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
        except ValueError:
            pass
						      
    def getCommandArgs(self):
        argstr=''
        if not (self.ui.checkBox_magn.isChecked() and self.ui.checkBox_pha.isChecked()
                and  self.ui.checkBox_ksp.isChecked() and self.ui.checkBox_reimag.isChecked()):
            print "Forcing magnitude export since no export checkboxes enabled."
            tmp_msg = "Forcing magnitude export since no export checkboxes enabled."
            reply = QtGui.QMessageBox.question(self, 'Message', 
                    tmp_msg, QtGui.QMessageBox.Ok)
            argstr=' -m'
        if self.ui.checkBox_magn.isChecked():
            argstr=' -m'
        if self.ui.checkBox_pha.isChecked():
            argstr+=' -p'
        if self.ui.checkBox_ksp.isChecked():
            argstr+=' -k'
        if self.ui.checkBox_reimag.isChecked():
            argstr+=' -r'
        
        if self.ui.checkBox_gaussian.isChecked():
            argstr+=' -g %s -j %s' % (str(self.ui.lineEdit_gsigma.text()),str(self.ui.lineEdit_gorder.text()))
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
            print(cmd)
            os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
										
    def CheckFID(*args):  #send_button):
        try:
            
            output_dir = self.ui.lineEdit_dicompath2.text()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'check path: %s' % thispath
            cmd1 = os.path.join(thispath,'dcheck.sh') + ' -o ' + str(output_dir)
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            cmd1 ='[ -x mrview ] && mrview '+ str(output_dir)
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
										    
    def ViewFID(self):
        try:
            output_dir = self.ui.lineEdit_dicompath2.text()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            if os.path.isdir(os.path.join(output_dir,'tmp')):
                tmp_msg = "Tmp folder found in DICOM directory. You cannot view the dicoms if there are subdirectories with additional dicoms. Are you sure you want to delete the temporary 2D dicom files?"
                reply = QtGui.QMessageBox.question(self, 'Message', 
                     tmp_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

                if reply == QtGui.QMessageBox.Yes:
                    import shutil
                    shutil.rmtree(os.path.join(output_dir,'tmp')) 
                else:
                    return
            
            cmd1 = 'vglrun mrview  ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
											
    def Send2Daris2(*args):
        try:
            daris_ID = self.ui.lineEdit_darisid2.text()
            dicom_dir = self.ui.lineEdit_dicompath2.text()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd = os.path.join(thispath,'dpush') + ' -c ' + str(daris_ID) + ' -s mf-erc ' + str(dicom_dir)
            os.system(cmd)	
        except ValueError:
            pass
											    
    def UpdateGUI(self):
        self.ui.setWindowTitle(_translate("Form", "MBI\'s Agilent to Dicom converter application ("+__version__+")", None))
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
												  


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main = Agilent2DicomWindow()
    main.show()
    sys.exit(app.exec_())
												  

#!/usr/bin/env python
# AgilentImageProcessing GUI for Agilent 9.4T MR FDF/FID image processing
#
# $Header$
# $Id$
# $Date$
#
# Version 1.2.5: Working version on Redhat Workstation
# Version 1.3.0: Info tab panels show information from Procpar
# Version 1.6.0: Tabs for epanechnikov, fourier gauss, fourier epanechnikov
# Version 1.7.0:
# Version 1.8.0: Fourier domain, double resolution
# Version 1.8.1:
# Version 1.8.2: CleanUpDicom dialog
#   Name change: AgilentImageProcessing replaces Agilent2DicomQt2
# Version 2.0.0: Phase and non-local means image processing tools included in
#                new tab
# Version 2.0.1: Non-local means options and allows use of FDF/FID and nifti
# Version 2.0.2: Non-local means parameters and better source code package
# Version 2.0.3: B1 bias and coil sensitivity correction implementation
#
# Copyright 2015 Michael Eager
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
pyuic4 --output Agilent2DicomQt2.py Agilent2DicomQt2.ui
"""


import os
import subprocess
import sys
import re
import numpy as np
from PyQt4 import Qt, QtGui, QtCore
from PyQt4.QtGui import QDialog, QFileDialog, QApplication
from Agilent2DicomQt2 import Ui_MainWindow, _translate
from ReadProcpar import ReadProcpar, ProcparInfo
from agilent2dicom_globalvars import AGILENT2DICOM_APP_VERSION, AGILENT2DICOM_VERSION
import logging
from dcmcleanup import Ui_CleanUpDialog
DEBUGGING = 1
# Agilent2DicomAppVersion=0.7
__author__ = "Michael Eager, Monash Biomedical Imaging"
__version__ = str(AGILENT2DICOM_APP_VERSION)
__date__ = "$Date$"
__copyright__ = "Copyright 2014-2015 Michael Eager"


Agilent2DicomAppStamp = re.sub(
    r'\$Id$', r'\1', "$Id$")
cmd_header = 'source ~/.bashrc; if test ${MASSIVE_USERNAME+defined}; \n\
then \n\
   echo \"On Massive\" \n\
   module purge \n\
   module load massive virtualgl\n\
   module load python\n\
   module load dcmtk mrtrix \n\
   #module list \n\
   echo $PYTHONPATH \n\
else \n\
echo \"Not in MASSIVE\" \n\
fi \n\
echo \"Done\"; '

proc_header = 'source ~/.bashrc; if test ${MASSIVE_USERNAME+defined}; \n\
then \n\
   echo \"On Massive\" \n\
   module purge \n\
   module load massive virtualgl\n\
   module load matlab \n\
   module list \n\
else \n\
   echo \"Checking whether matlab in PATH.\" \n\
   if test -x `which matlab`; \n\
   then \n\
      echo \"MATLAB found.\"; \n\
   else \n\
      echo \"Cannot find MATLAB in PATH. Close program and exit PATH in local bashrc file. Exiting shell..\"; \n\
       exit 1; \n\
   fi \n\
fi \n\
echo \"Done\"; '

mrview_header = 'source ~/.bashrc;if test ${MASSIVE_USERNAME+defined}; \n\
then \n\
   echo "On Massive" \n\
   module purge \n\
   module load massive virtualgl\n\
   module load mrtrix3 \n\
   module list \n\
   GL = vglrun\n\
else GL= \n\
fi; $GL mrview '

inputdir = []


class CleanUpDialog(QtGui.QDialog, Ui_CleanUpDialog):

    def __init__(self, parent=None, inputdir=None):
        QtGui.QDialog.__init__(self, parent)
        #self.ui = Ui_CleanUpDialog
        self.setupUi(self, inputdir)


class Agilent2DicomWindow(QtGui.QMainWindow):

    """Agilent2DicomWindow GUI for FDF and FID converter
    """
    niftiflag = 0  # save to nifti flag

    def __init__(self):
        super(Agilent2DicomWindow, self).__init__()
        self.ui = Ui_MainWindow()
        # Set up the user interface from Designer.
        self.ui.setupUi(self)
        # Make some local modifications.
        # self.colorDepthCombo.addItem("2 colors (1 bit per
        self.ui.Version_tag.setText(
            _translate("MainWindow", "Agilent2Dicom v" + AGILENT2DICOM_APP_VERSION, None))

        logging.basicConfig(format='%(levelname)s:%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            filename='qtapp-agilent2dicom.log', level=logging.DEBUG)
        logging.info('Starting AgilentImageProcessing')
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
        self.ui.checkBox_stdev_cplx.setChecked(False)
        self.ui.checkBox_stdev_phase.setChecked(False)
        self.ui.checkBox_stdev_magn.setChecked(False)

        # # Connect up the buttons.
        # self.connect(self.ui.buttonBox, Qt.SIGNAL("accepted()"),self.accept)
        # self.connect(self.ui.buttonBox , Qt.SIGNAL("rejected()"),self.reject)
        self.ui.pushButton_changefdf.clicked.connect(self.ChangeFDFpath)
        self.ui.pushButton_changedicom.clicked.connect(self.ChangeFDFDicomPath)
        self.ui.pushButton_changefid.clicked.connect(self.ChangeFIDpath)
        self.ui.pushButton_changedicom2.clicked.connect(
            self.ChangeFIDDicomPath)

        self.ui.pushButton_convert.clicked.connect(self.ConvertFDF)
        self.ui.pushButton_check.clicked.connect(self.CheckFDF)
        self.ui.pushButton_view.clicked.connect(self.ViewFDF)
        self.ui.pushButton_send2daris.clicked.connect(self.Send2Daris)
        self.ui.pushButton_CleanUpDicoms.clicked.connect(self.CleanUpDicoms)
        self.ui.pushButton_convertfid.clicked.connect(self.ConvertFID)
        self.ui.pushButton_check2.clicked.connect(self.CheckFID)
        self.ui.pushButton_view2.clicked.connect(self.ViewFID)
        self.ui.pushButton_send2daris2.clicked.connect(self.Send2DarisFID)

        self.ui.pushButton_processfolder1.clicked.connect(self.ChangeProcPath1)
        self.ui.pushButton_processfolder2.clicked.connect(self.ChangeProcPath2)
        self.ui.pushButton_ChangeNifti1.clicked.connect(self.ChangeProcFile1)
        self.ui.pushButton_ChangeNifti2.clicked.connect(self.ChangeProcFile2)
        self.ui.pushButton_ClearProc1.clicked.connect(self.ClearProcLine1)
        self.ui.pushButton_ClearProc2.clicked.connect(self.ClearProcLine2)

        self.ui.pushButton_CoilSensChangeDir.clicked.connect(
            self.ChangeCoilSensDir)
        self.ui.pushButton_CoilSensChangeFile.clicked.connect(
            self.ChangeCoilSensFile)
        self.ui.pushButton_CoilSensClear.clicked.connect(self.ClearProcLine3)
        self.ui.pushButton_procout.clicked.connect(self.ChangeProcOutputPath)

        self.ui.pushButton_SWI.clicked.connect(self.ProcSWI)
        self.ui.pushButton_MEE.clicked.connect(self.ProcMEE)
        self.ui.pushButton_MCI.clicked.connect(self.ProcMCI)
        self.ui.pushButton_NLpipeline1.clicked.connect(self.ProcNLpipeline1)
        self.ui.pushButton_NLpipeline2.clicked.connect(self.ProcNLpipeline2)
        self.ui.pushButton_NLpipeline3.clicked.connect(self.ProcNLpipeline3)
        
#        self.ui.checkBox_NLinputMagPha.setChecked(False)
        self.ui.checkBox_NLinputRI.setChecked(False)
      #  QtCore.QObject.connect(self.ui.actionPRINLM,
      #                         QtCore.SIGNAL(
      #                             QtCore.QString.fromUtf8("triggered()")),
      #                         self.SetSearchArea5)
        QtCore.QObject.connect(self.ui.actionSave_Filter_Outputs_to_Nifti,
                               QtCore.SIGNAL(
                                   QtCore.QString.fromUtf8("triggered()")),
                               self.toggleNifti)
        QtCore.QObject.connect(self.ui.actionAbout,
                               QtCore.SIGNAL(
                                   QtCore.QString.fromUtf8("triggered()")),
                               self.About)
        QtCore.QObject.connect(self.ui.actionHelp,
                               QtCore.SIGNAL(
                                   QtCore.QString.fromUtf8("triggered()")),
                               self.About)
#        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def accept(self):
        '''Execute the command in response to the OK button.'''
        self.close()

    def reject(self):
        '''Cancel.'''
        self.close()

    def toggleNifti(self):
        """Toggle the nifti flag
        """
        self.niftiflag = (self.niftiflag + 1) % 2
        logging.info('Nifti toggled')

    def SetSearchArea3(self):
        self.ui.lineEdit_NLsearcharea.setText('3')
    def SetSearchArea5(self):
        self.ui.lineEdit_NLsearcharea.setText('5')
            
        
    # @QtCore.pyqtSlot()
    def ChangeFDFpath(self):
        '''ChangeFDFpath
        '''
        success = 1
        try:
            newdir = str(
                QFileDialog.getExistingDirectory(self, '''Select FDF Directory'''))
            print newdir
            self.ui.lineEdit_fdfpath.setText(newdir)
            files = os.listdir(newdir)
            if 'procpar' not in files:
                print "ChangeFDFPath Error: FDF folder doesn't contain a procpar file"
                success = 0
            fdffiles = [f for f in files if f.endswith('.fdf')]
            if len(fdffiles) == 0:
                print 'Error: FDF folder does not contain any fdf files'
                logging.info(
                    'ChangeFDFpath folder does not contain any fdf files')
                success = 0
            if success:
                if re.search('img', newdir):
                    out = re.sub('img', 'dcm', newdir)
                else:
                    out = newdir + '.dcm'

                self.ui.lineEdit_dicompath.setText(out)
                self.ui.lineEdit_darisid.setText(self.GetDarisID(newdir))
            self.ui.checkBox_gaussian3D.setChecked(False)
            self.ui.checkBox_gaussian2D.setChecked(False)
            self.ui.checkBox_median.setChecked(False)
            self.ui.checkBox_wiener.setChecked(False)
            self.UpdateGUI()
            logging.info('ChangeFDFpath ' + newdir)
        except ValueError:
            logging.info('ChangeFDFpath ValueError')

    def ChangeFDFDicomPath(self):
        """ChangeFDFDicomPath
        """
        try:
            newdir = str(QFileDialog.getExistingDirectory(self,
                                                          '''Select
                                                          DICOM
                                                          Directory'''))
            files = os.listdir(newdir)
            dcmfiles = [f for f in files if f.endswith('.dcm')]
            if len(dcmfiles) == 0:
                print 'FDF output DICOM folder does not contain any dcm files'
            else:
                quit_msg = "Are you sure you want to delete existing dicoms?"
                reply = QtGui.QMessageBox.question(self, 'Message', quit_msg,
                                                   QtGui.QMessageBox.Yes,
                                                   QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:
                    #                    event.accept()
                    self.ui.lineEdit_dicompath.setText(newdir)
                    self.ui.lineEdit_darisid.setText(
                        self.GetDarisID(str(self.ui.lineEdit_fdfpath.text())))
#                else:
#                    event.ignore()
            self.UpdateGUI()
            logging.info('ChangeFDFDicompath ' + newdir)
        except ValueError:
            logging.info('ChangeFDFDicompath error')

    def ChangeFIDpath(self):
        """ChangeFIDpath
        """
        success = 1
        try:
            newdir = str(
                QFileDialog.getExistingDirectory(self, "Select FID Directory"))
            self.ui.lineEdit_fidpath.setText(newdir)
            files = os.listdir(newdir)
            if 'procpar' not in files:
                print '''ChangeFIDPath Error: FID folder does not contain a
                         procpar file'''
                success = 0
            fidfiles = [f for f in files if f.endswith('fid')]
            if len(fidfiles) == 0:
                print 'Error: FID folder does not contain any fid files'
                success = 0
            if success:
                if re.search('fid', newdir):
                    out = re.sub('fid', 'dcm', newdir)
                else:
                    out = newdir + '.dcm'
                self.ui.lineEdit_dicompath2.setText(out)
                self.ui.lineEdit_darisid2.setText(self.GetDarisID(newdir))
            self.ui.checkBox_gaussian3D.setChecked(False)
            self.ui.checkBox_gaussian2D.setChecked(False)
            self.ui.checkBox_median.setChecked(False)
            self.ui.checkBox_wiener.setChecked(False)
            self.UpdateGUI()
            logging.info('ChangeFIDpath ' + newdir)
        except ValueError:
            logging.info('ChangeFIDpath error')

    def ChangeFIDDicomPath(self):
        """ChangeFIDDicomPath
        """
        try:
            newdir = str(QFileDialog.getExistingDirectory(self,
                                                          '''Select or
                                                          create new
                                                          DICOM
                                                          Directory
                                                          '''))
            files = os.listdir(newdir)
            dcmfiles = [f for f in files if f.endswith('.dcm')]
            if len(dcmfiles) == 0:
                print 'FDF output DICOM folder does not contain any dcm files'
            else:
                quit_msg = "Are you sure you want to delete existing dicoms?"
                reply = QtGui.QMessageBox.question(self, 'Message',
                                                   quit_msg,
                                                   QtGui.QMessageBox.Yes,
                                                   QtGui.QMessageBox.No)

                if reply == QtGui.QMessageBox.Yes:
                    self.ui.lineEdit_dicompath2.setText(newdir)
                    self.ui.lineEdit_darisid.setText(
                        self.GetDarisID(str(self.ui.lineEdit_fidpath.text())))
            self.UpdateGUI()
            logging.info('ChangeFIDDicompath ' + newdir)
        except ValueError:
            logging.info('ChangeFIDDicompath error')

    def ConvertFDF(self):
        """ConvertFDF run FDF to dicom script
        """
        try:
            input_dir = self.ui.lineEdit_fdfpath.text()
            output_dir = self.ui.lineEdit_dicompath.text()
            thispath = str(
                os.path.dirname(os.path.realpath(os.path.abspath(__file__))))
            print 'fdf2dcm path: %s' % thispath
            cmd1 = str(os.path.join(thispath, 'fdf2dcm.sh')) + \
                ' -i ' + str(input_dir) + ' -o ' + str(output_dir)
            if self.ui.checkBox_nodcmulti.isChecked():
                cmd1 += ' -d '
            if self.ui.checkBox_debugging.isChecked():
                cmd1 += ' -v '
            print cmd1
            cmd = cmd_header + cmd1
            logging.info('Convert FDF :' + cmd)
            print subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True,
                                   executable="/bin/bash").stdout.read()
            self.UpdateGUI()
            logging.info('Convert FDF : done')
        except OSError:
            logging.error('Send2Daris OSError', exc_info=True)
        except ValueError:
            logging.error(
                'ConvertFDF value error - invalid arguments to Popen.', exc_info=True)

    def GetDarisID(self, inpath):
        """GetDarisID read Daris ID from procpar
        """
        daris_id = ''
        try:
            print "GetDarisID %s " % inpath
            procpar, ptext = ReadProcpar(
                os.path.join(os.path.abspath(str(inpath)), 'procpar'))
            if 'name' in procpar.keys():
                if re.search('DaRIS', procpar['name']):
                    daris_id = re.sub(r'DaRIS\^', '', procpar['name'])
        except ValueError:
            logging.error("Value Error in GetDarisID", exc_info=True)

        return daris_id

    def CheckDicomDir(self, dpath):
        """Check the output dicom directries
        """
        if dpath != "" and os.path.isdir(dpath):
            # Check folder contains procpar and *.fdf files
            files = os.listdir(dpath)

            dcmfiles = [f for f in files if f.endswith('.dcm')]
            if len(dcmfiles) == 0:
                print 'AgilentImageProcessing Error: DICOM folder does not contain any dcm files'
                return False
            else:
                return True
        else:
            # print 'Error: DICOM folder does not exist: %s ' % dpath
            return False

    def Send2Daris(self):
        """Send DICOMs to DaRIS
        """
        try:
            daris_ID = self.ui.lineEdit_darisid.text()
            if str(daris_ID) == "":
                print "No DaRIS ID"
                QtGui.QMessageBox.warning(self, 'Warning', '''Cannot send to
                                          DaRIS without a proper ID.''')
                logging.warning('''Send2Daris Cannot send to
                                DaRIS without a proper ID.''')
                return
            dicom_dir = str(self.ui.lineEdit_dicompath.text())
            if dicom_dir == "" or not os.path.isdir(dicom_dir):
                print "No Dicom directory"
                QtGui.QMessageBox.warning(self, 'Warning', '''Cannot send to
                DaRIS without a proper Dicom path.''')
                logging.warning('''Send2Daris Cannot send to
                                DaRIS without a proper Dicom path.''')
                return

            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            cmd = os.path.join(thispath, 'dpush') + ' -c ' +\
                str(daris_ID) + ' -s mf-erc ' + str(dicom_dir)
            send_msg = '''Are you sure you want to send to DaRIS the dicom
            directory\n''' + str(dicom_dir) + "\n using the ID: " +\
                str(daris_ID) + "?"
            reply = QtGui.QMessageBox.question(self, 'Message',
                                               send_msg,
                                               QtGui.QMessageBox.Yes,
                                               QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.Yes:
                logging.info('Send2Daris ' + cmd)
                try:
                    print subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True,
                                           executable="/bin/bash").stdout.read()

                    print "dpush completed returncode" + subprocess.Popen.returncode
                    logging.info('Send2Daris ... completed.')
                except ValueError:
                    logging.error(
                        'Send2Daris invalid arguments to Popen', exc_info=True)
                except OSError:
                    logging.error('Send2Daris OSError', exc_info=True)

            else:
                logging.info('Send2Daris sending cancelled')
                print 'Send2Daris sending cancelled'
            self.UpdateGUI()
        except ValueError:
            logging.error('Send2Daris value error', exc_info=True)
            print 'Send2Daris value error'

    def CheckFDF(self):  # send_button):
        """ Run mrinfo and dcheck on dicoms
        """
        try:
            output_dir = str(self.ui.lineEdit_dicompath.text())
            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            cmd1 = 'mrinfo ' + output_dir
            print cmd1
            cmd = cmd_header + cmd1
            logging.info('CheckFDF ' + cmd)
            print subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True,
                                   executable="/bin/bash").stdout.read()
            cmd1 = os.path.join(thispath, 'dcheck.sh') + ' -o ' + output_dir
            logging.info('CheckFDF ' + cmd)
            print subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True,
                                   executable="/bin/bash").stdout.read()
            self.UpdateGUI()
        except ValueError:
            logging.error(
                'CheckFDF error. Possibly invalid arguments to Popen', exc_info=True)
        except OSError:
            logging.error('CheckFDF OSError', exc_info=True)

    def ViewFDF(self):
        """View image with mrview
        """
        try:
            output_dir = str(self.ui.lineEdit_dicompath.text())
            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            cmd1 = mrview_header + output_dir
            # print cmd1
            logging.info('ViewFDF ' + cmd1)
            print subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True,
                                   executable="/bin/bash").stdout.read()
        except ValueError:
            logging.error('ViewFDF error ', exc_info=True)
        except OSError:
            logging.error('ViewFDF OSError', exc_info=True)

    def SigmaString(self, sigmatext, combotext, source):
        """Calculate sigma
        combotext: contents of the combobox scaling
        source: must be a valid FDF or FID path with procpar file

        output: comma-separtated string without spaces
        """
        sigmafactor = np.array((1, 1, 1))
        print 'SigmaString:  %s %s' % (sigmatext, combotext)
        print 'Sigma string source', source
        try:
            if combotext != "unit voxel":
                logging.info('SigmaString ' + combotext)
                procpar, ptext = ReadProcpar(os.path.join(source, 'procpar'))
                pinfo = ProcparInfo(procpar)
                if combotext == "in mm":
                    sigmafactor = np.array(pinfo['Voxel_Res_mm'])
                elif combotext == "in um":
                    sigmafactor = np.array(pinfo['Voxel_Res_mm']) * 1000
                else:
                    print "SigmaString combobox not of 'unit voxel', 'in mm', or 'in um'"
                    logging.info(
                        "SigmaString combobox not of 'unit voxel', 'in mm', or 'in um'")
            logging.info('SigmaString ' + str(sigmafactor))
        except ValueError:
            logging.info('SigmaString combobox error')

        if len(sigmafactor) != 3:
            print "Sigma units needs to be three elements"
            logging.info('SigmaString units needs to be three elements')
            sigmafactor[1] = sigmafactor[0]
            sigmafactor[2] = sigmafactor[0]
        try:
            # if s.find(', '):
            x, y, z = map(float, sigmatext.split(','))
            argstr = ' %f,%f,%f' % (x / sigmafactor[0],
                                    y / sigmafactor[1],
                                    z / sigmafactor[2])
        except ValueError:
            argstr = ' %f,%f,%f' % (float(sigmatext) / sigmafactor[0],
                                    float(sigmatext) / sigmafactor[1],
                                    float(sigmatext) / sigmafactor[2])
        print 'SigmaString: factor %f,%f,%f : args %s' % (sigmafactor[0],
                                                          sigmafactor[1],
                                                          sigmafactor[2],
                                                          argstr)
        return argstr

    def getCommandArgs(self):
        """getCommandArgs
        """
        argstr = ''

        if self.ui.checkBox_magn.isChecked():
            argstr = ' -m'
        if self.ui.checkBox_pha.isChecked():
            argstr += ' -p'
        if self.ui.checkBox_ksp.isChecked():
            argstr += ' -k'
        if self.ui.checkBox_reimag.isChecked():
            argstr += ' -r'
        if not (self.ui.checkBox_magn.isChecked() and
                self.ui.checkBox_pha.isChecked() and
                self.ui.checkBox_ksp.isChecked() and
                self.ui.checkBox_reimag.isChecked()):
            print "Forcing magn export since no export checkboxes enabled."
            argstr = ' -m'
        print argstr
        # Complex Gaussian 3D
        if self.ui.checkBox_gaussian2D.isChecked():
            argstr += ' -2'
        if self.ui.checkBox_gaussian3D.isChecked() or \
           self.ui.checkBox_gaussian2D.isChecked():
            curUnits = str(self.ui.comboBox_gauss_sigmascale.currentText())
            sigma = self.SigmaString(str(self.ui.lineEdit_gsigma.text()),
                                     curUnits,
                                     str(self.ui.lineEdit_fidpath.text()))
            argstr += ' -g %s' % sigma
            argstr += ' -j %s' % str(self.ui.lineEdit_gorder.text())
            if self.ui.nearest.isChecked():
                argstr += ' -e nearest'
            elif self.ui.reflect.isChecked():
                argstr += ' -e reflect'
            elif self.ui.mirror.isChecked:
                argstr += ' -e mirror'
            elif self.ui.wrap.isChecked:
                argstr += ' -e wrap'
            else:
                argstr += ' -e reflect'

        # Fourier Gaussian
        if self.ui.checkBox_kspgaussian.isChecked():
            curUnits = str(self.ui.comboBox_kspgauss_sigunit.currentText())
            sigma = self.SigmaString(str(self.ui.lineEdit_gfsigma.text()),
                                     curUnits,
                                     str(self.ui.lineEdit_fidpath.text()))
            argstr += ' -G %s ' % sigma
            if self.ui.checkBox_kspgauss_super.isChecked():
                argstr += ' -D'
            if self.ui.checkBox_kspgaussshift.isChecked():
                argstr += ' -C '
        # Complex Median
        if self.ui.checkBox_median.isChecked():
            argstr += ' -n %s ' % (str(self.ui.lineEdit_median_size.text()))

        # Complex Wiener
        if self.ui.checkBox_wiener.isChecked():
            argstr += ' -w %s ' % (str(self.ui.lineEdit_wiener_size.text()))
            if self.ui.lineEdit_wiener_noise.text():
                argstr += ' -z %s ' % (
                    str(self.ui.lineEdit_wiener_noise.text()))

        # Complex Stdev
        if self.ui.checkBox_stdev_cplx.isChecked():
            argstr += ' -s 0/%s ' % (str(self.ui.lineEdit_stdev_size.text()))
        # Magn Stdev
        if self.ui.checkBox_stdev_magn.isChecked():
            argstr += ' -s 1/%s ' % (str(self.ui.lineEdit_stdev_size.text()))
        # Phase Stdev
        if self.ui.checkBox_stdev_phase.isChecked():
            argstr += ' -s 2/%s ' % (str(self.ui.lineEdit_stdev_size.text()))

        # Epanechnikov
        if self.ui.checkBox_epanechnikov2D.isChecked():
            argstr += ' -2'
        if self.ui.checkBox_epanechnikov3D.isChecked() or \
           self.ui.checkBox_epanechnikov2D.isChecked():
            curUnits = str(self.ui.comboBox_epabandwidth.currentText())
            sigma = self.SigmaString(str(self.ui.lineEdit_epaband.text()),
                                     curUnits,
                                     str(self.ui.lineEdit_fidpath.text()))
            argstr += ' -y %s' % sigma
            # argstr += ' -j %s' % str(self.ui.lineEdit_gorder.text())
            # if self.ui.nearest_epa.isChecked():
            #     argstr += ' -e nearest'
            # elif self.ui.mirror_epa.isChecked():
            #     argstr += ' -e mirror'
            # elif self.ui.reflect_epa.isChecked():
            #     argstr += ' -e reflect'
            # elif self.ui.wrap_epa.isChecked:
            #     argstr += ' -e wrap'
            # else:
            #     argstr += ' -e reflect'

        # Fourier Epanechnikov
        if self.ui.checkBox_kspepa.isChecked():
            curUnits = str(self.ui.comboBox_kspepa_scaleunit.currentText())
            sigma = self.SigmaString(str(self.ui.lineEdit_kspepa_band.text()),
                                     curUnits,
                                     str(self.ui.lineEdit_fidpath.text()))
            argstr += ' -Y %s' % sigma
            # argstr += ' -j %s' % str(self.ui.lineEdit_gorder.text())
            if self.ui.checkBox_kspepashift.isChecked():
                argstr += ' -C '
            if self.ui.checkBox_kspgauss_super.isChecked():
                argstr += ' -D '

        if self.niftiflag == 1:
            argstr += ' -N'
        logging.info('Command Args ' + argstr)
        return argstr

    def ConvertFID(self):
        """ConvertFID run FID2dicom script
        """
        try:
            input_dir = str(self.ui.lineEdit_fidpath.text())
            output_dir = str(self.ui.lineEdit_dicompath2.text())
            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            print 'fid2dcm path: %s' % thispath
            cmd1 = os.path.join(thispath, 'fid2dcm.sh')
            cmd1 = cmd1 + ' -v ' + self.getCommandArgs()
            cmd1 = cmd1 + ' -i ' + str(input_dir) + ' -o ' + str(output_dir)
            print cmd1
            cmd = cmd_header + cmd1
            # print cmd
            logging.info('ConvertFID ' + cmd)
        except ValueError:
            logging.error('ConvertFID error', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        try:
            print subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True,
                                   executable="/bin/bash").stdout.read()
            logging.info('Send2Daris ... completed.')
        except ValueError:
            logging.error(
                'ConvertFID invalid arguments to Popen', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        except OSError:
            logging.error('ConvertFID OSError', exc_info=True)
            pass
        self.UpdateGUI()
        print 'ConvertFID completed.'
        logging.info('ConvertFID done')
        # except ValueError:
        #    print 'ConvertFID error.'
        #    logging.error('ConvertFID error', exc_info=True)

    def CheckFID(self):
        """CheckFID check outputs of fid2dicom
        """
        try:
            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            print 'check path: %s' % thispath
            dicom_raw_dir = str(self.ui.lineEdit_dicompath2.text())
            # Shorten string if it ends in '/'
            if dicom_raw_dir[-1] == '/':
                dicom_raw_dir = dicom_raw_dir[:-1]
            # Get tuple of root path and raw/unfiltered dicom path
            (dicom_dir_root, dicom_raw) = os.path.split(dicom_raw_dir)
            # Do a regex and get all the dicom paths produced by Agilent2Dicom
            rgx = re.compile(r'' + re.sub('.dcm', '', dicom_raw) + ".*.dcm")
            logging.info('CheckFID ' + str(filter(rgx.match,
                                                  os.listdir(dicom_dir_root))))
            for dicom_dir in filter(rgx.match, os.listdir(dicom_dir_root)):
                # os.system("ls " + dicom_dir_root + " | grep '" +
                # re.sub('.dcm', '', dicom_raw) + ".*.dcm'"):
                dcmpath = os.path.join(dicom_dir_root, dicom_dir)
                if not os.path.isdir(dcmpath) or len(os.listdir(dcmpath)) <= 2:
                    cmd1 = os.path.join(
                        thispath, 'dcheck.sh') + ' -o ' + str(dcmpath)
                    # print cmd1
                    # cmd = cmd_header + cmd1
                    # print cmd
                    print subprocess.Popen(cmd1,
                                           stdout=subprocess.PIPE, shell=True,
                                           executable="/bin/bash").stdout.read()
            self.UpdateGUI()
            print 'CheckFID done'
            logging.info('CheckFID done')
        except OSError:
            logging.error('CheckFID OSError', exc_info=True)
        except ValueError:
            print 'CheckFID error'
            logging.error('CheckFID error', exc_info=True)

    def ViewFID(self):
        """ViewFID show dicom outputs with mrview
        """
        try:
            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            dicom_raw_dir = str(self.ui.lineEdit_dicompath2.text())
            # Shorten string if it ends in '/'
            if dicom_raw_dir[-1] == '/':
                dicom_raw_dir = dicom_raw_dir[:-1]
            # Get tuple of root path and raw/unfiltered dicom path
            (dicom_dir_root, dicom_raw) = os.path.split(dicom_raw_dir)
            # Do a regex and get all the dicom paths produced by Agilent2Dicom
            rgx = re.compile(r'' + re.sub('.dcm', '', dicom_raw) + ".*.dcm")
            for dicom_dir in filter(rgx.match, os.listdir(dicom_dir_root)):
                # os.system("ls " + dicom_dir_root + " | grep '" +
                # re.sub('.dcm', '', dicom_raw) + ".*.dcm'"):
                dcmpath = os.path.join(dicom_dir_root, dicom_dir)
                if not os.path.isdir(dcmpath) or len(os.listdir(dcmpath)) <= 2:
                    cmd1 = mrview_header + dcmpath + ' & '
                    print cmd1
                    print subprocess.Popen(cmd1,
                                           stdout=subprocess.PIPE, shell=True,
                                           executable="/bin/bash").stdout.read()

            self.UpdateGUI()
            logging.info('ViewFID done')
        except OSError:
            logging.error('ViewFID OSError', exc_info=True)
        except ValueError:
            logging.error('ViewFID error', exc_info=True)

    def Send2DarisFID(self):
        """Send2DarisFID
        """
        try:
            daris_ID = self.ui.lineEdit_darisid2.text()
            if str(daris_ID) == "":
                print "Send2DarisFID Cannot send to DaRIS without a proper ID. "
                QtGui.QMessageBox.warning(self, 'Warning', '''
                     Cannot send to DaRIS without a proper ID.''')
                logging.warning('''Send2DarisFID Cannot send to
                                DaRIS without a proper ID.''')
                return
            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))

            dicom_raw_dir = str(self.ui.lineEdit_dicompath2.text())
            # Shorten string if it ends in '/'
            if dicom_raw_dir[-1] == '/':
                dicom_raw_dir = dicom_raw_dir[:-1]
            # Get tuple of root path and raw/unfiltered dicom path
            (dicom_dir_root, dicom_raw) = os.path.split(dicom_raw_dir)
            # Do a regex and get all the dicom paths produced by Agilent2Dicom
            rgx = re.compile(r'' + re.sub('.dcm', '', dicom_raw) + ".*.dcm")
            for dicom_dir in filter(rgx.match, os.listdir(dicom_dir_root)):
                # os.system("ls " + dicom_dir_root + " | grep '" +
                # re.sub('.dcm', '', dicom_raw) + ".*.dcm'"):
                dcmpath = os.path.join(dicom_dir_root, dicom_dir)
                if not os.path.isdir(dcmpath) or len(os.listdir(dcmpath)) <= 2:
                    QtGui.QMessageBox.warning(self,
                                              'Warning', '''Cannot send
                                              to DaRIS. Directory'''
                                              + dicom_dir
                                              + " is empty")
                    logging.warning('''Send2Daris Cannot send to
                                DaRIS without a proper Dicom path.
                                    ''' + dicom_dir + " is empty")
                    print 'Send2Daris Cannot send to DaRIS without a proper Dicom path.'
                else:
                    cmd = os.path.join(thispath, 'dpush') + ' -c ' +\
                        str(daris_ID) + ' -s mf-erc ' + str(dcmpath)
                    send_msg = '''Are you sure you want to send
                    toDaRIS the dicom directory\n

                    ''' + str(dcmpath) + "\n using the ID: " +\
                        str(daris_ID) + " ?"
                    reply = QtGui.QMessageBox.question(self, 'Message',
                                                       send_msg,
                                                       QtGui.QMessageBox.Yes,
                                                       QtGui.QMessageBox.No)
                    if reply == QtGui.QMessageBox.Yes:
                        logging.info('Send2Daris ' + cmd)
                        print subprocess.Popen(cmd,
                                               stdout=subprocess.PIPE,
                                               shell=True,
                                               executable="/bin/bash").stdout.read()
                        logging.info('Send2DarisFID completed.')
                        print 'Send2DarisFID completed.'
                    else:
                        logging.info('Send2DariFIDs cancelled')
                        print 'Send2DarisFID cancelled'

        except ValueError:
            print 'Send2DarisFID ValueError'
            logging.error('Send2DarisFID error', exc_info=True)
        except OSError:
            logging.error('Send2DarisFID OSError', exc_info=True)

    def CleanUpDicoms(self):
        parentdir = os.path.dirname(str(self.ui.lineEdit_dicompath2.text()))
        dialog = CleanUpDialog(self, inputdir=parentdir)
        if dialog.exec_():
            results = dialog.getValues()
            for i in xrange(0, len(results)):
                print "Deleting ", results[i]
                import shutil
                shutil.rmtree(os.path.join(parentdir, results[i]))

    def UpdateGUI(self):
        """Update the GUI
        """
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
        """Update metadata info panel in FDF tab
        """
        procpar, procpartext = ReadProcpar(
            os.path.join(str(self.ui.lineEdit_fdfpath.text()), 'procpar'))
        SequenceText = "FDF image:"

        p = ProcparInfo(procpar)
        self.ui.FDFprocparInfo.setText(SequenceText +
                                       '\n'.join("%s:\t\t%r" %
                                                 (key, val) for (key, val) in
                                                 p.iteritems()))

    def FIDInfoUpdate(self):
        """Update metadata info panel in FID tab
        """
        procpar, procpartext = ReadProcpar(
            os.path.join(str(self.ui.lineEdit_fidpath.text()),
                         'procpar'))
        # ds, MRAcq_type=ProcparToDicomMap(procpar)
        try:
            ftext = open(
                os.path.join(str(self.ui.lineEdit_fidpath.text()),
                             'text'), 'r')
            SequenceText = "FID image:\n%s" % ftext.readline()
            ftext.close()
        except ValueError:
            SequenceText = "FID image:\n-no text-"
            pass
        # DisplayText=[ SequenceText + "\n", "StudyID: %s\n" %
        #  ds.StudyID, "Patient ID: %s\n" % ds.PatientID, "Protocol
        #  name: %s\n"% ds.ProtocolName, "Series Desc: %s\n" %
        #  ds.SeriesDescription, "Image type: %s\n" % ds.ImageType,
        #  "MR Acquisition Type: %s\n" % MRAcq_type, "Acq Matrix: %s,
        #  %s\n" % (ds.Rows, ds.Columns), "Pixel Spacing: %s, %s\n" %
        #  (ds.PixelSpacing[0], ds.PixelSpacing[1]), "Slice thickness:
        #  %s\n" % ds.SliceThickness]
        # for column, key in
        #  enumerate(self.data): # for row, item in
        #  enumerate(self.data[key]):
        # newitem =
        #  QTableWidgetItem(item)
        # self.ui.listViewsetItem(row,
        #  column, newitem) print DisplayText
        #  QtGui.QMessageBox.question(self, 'Info',
        #  '\n'.join(DisplayText), QtGui.QMessageBox.Ok)

        p = ProcparInfo(procpar)
        self.ui.FIDprocparInfo.setText(
            SequenceText + '\n'.join("%s:\t\t%r" % (key,
                                                    val) for (key, val) in
                                     p.iteritems()))

    def About(self):
        """ Help dialog
        """
        msg = '''                      Agilent2Dicom
        
Agilent2Dicom converts FDF and FID images from the Agilent 9.4T MR scanner at
Monash Biomedical Imaging (MBI) into enhanced MR DICOM images.
        
Homepage:
https://confluence-vre.its.monash.edu.au/display/MBI/Agilent+FDF+to+Dicom+converter
App Version: %s
%s
Source Version: %s
Author: %s
Copyright: %s''' % (__version__, Agilent2DicomAppStamp, AGILENT2DICOM_VERSION, __author__, __copyright__)

        reply = QtGui.QMessageBox.question(self, 'Message', msg,
                                           QtGui.QMessageBox.Ok)

    def close(self):
        """close GUI
        """
        print "Closing Agilent2Dicom application."
        logging.info('Closing Agilent2DicomAppQt')

    def ChangeProcPath1(self):
        """ChangeProc1path
        """
        success = 1
        try:
            newdir = str(
                QFileDialog.getExistingDirectory(self, "Select FID or FDF Directory"))
            self.ui.lineEdit_ProcInfolder1.setText(newdir)
            files = os.listdir(newdir)
            if 'procpar' not in files:
                print '''ChangeFIDPath Error: FID folder does not contain a
                         procpar file'''
                success = 0
            fidfiles = [f for f in files if f.endswith('fid')]
            if len(fidfiles) == 0:
                print 'Input folder does not contain any fid files'

                fdffiles = [f for f in files if f.endswith('.fdf')]
                if len(fdffiles) == 0:
                    print 'Error: Folder does not contain any fid of fdf files'
                    success = 0

            if success:
                if re.search('fid', newdir):
                    out = re.sub('fid', 'nii', newdir)
                elif re.search('img', newdir):
                    out = re.sub('img', 'nii', newdir)
                else:
                    out = newdir + '.nii'
                self.ui.lineEdit_ProcOutputfolder.setText(out)

            # Reset checkboxes
            # self.ui.checkBox_gaussian3D.setChecked(False)
            # self.ui.checkBox_gaussian2D.setChecked(False)
            # self.ui.checkBox_median.setChecked(False)
            # self.ui.checkBox_wiener.setChecked(False)
            # self.UpdateGUI()
            logging.info('ChangeProc1path ' + newdir)
        except ValueError:
            logging.info('ChangeProc1path error')

    def ChangeProcPath2(self):
        """ChangeProc2path
        """
        success = 1
        try:
            newdir = str(
                QFileDialog.getExistingDirectory(self, "Select FID or FDF Directory"))
            self.ui.lineEdit_ProcInfolder2.setText(newdir)
            files = os.listdir(newdir)
            if 'procpar' not in files:
                print '''ChangeProc2Path Error: FID folder does not contain a
                         procpar file'''
                success = 0
            fidfiles = [f for f in files if f.endswith('fid')]
            if len(fidfiles) == 0:
                print 'Input folder does not contain any fid files'
                fdffiles = [f for f in files if f.endswith('.fdf')]
                if len(fdffiles) == 0:
                    print 'Error: Folder does not contain any fid of fdf files'
                    success = 0
                else:
                    print 'FDF folder found ', newdir
            else:
                print 'FID folder found ', newdir
            logging.info('ChangeProc2path ' + newdir)
        except ValueError:
            logging.info('ChangeProc2path error')

    def ChangeCoilSensDir(self):
        """ChangeCoilSensDir
        """
        success = 1
        try:
            newdir = str(
                QFileDialog.getExistingDirectory(self, "Select FID or FDF Directory"))
            self.ui.lineEdit_CoilSensDir.setText(newdir)
            files = os.listdir(newdir)
            if 'procpar' not in files:
                print '''ChangeCoilSensDir Error: FID folder does not contain a
                         procpar file'''
                success = 0
            fidfiles = [f for f in files if f.endswith('fid')]
            if len(fidfiles) == 0:
                print 'Input folder does not contain any fid files'
                fdffiles = [f for f in files if f.endswith('.fdf')]
                if len(fdffiles) == 0:
                    print 'Error: Folder does not contain any fid of fdf files'
                    success = 0
                else:
                    print 'FDF folder found ', newdir
            else:
                print 'FID folder found ', newdir
            logging.info('ChangeCoilSensDir ' + newdir)
        except ValueError:
            self.ui.lineEdit_CoilSensDir.setText("")
            logging.info('ChangeCoilSensDir error')

    def ChangeProcFile1(self):
        """ChangeProc1path
        """
        success = 1
        try:
            newnifti = str(
                QFileDialog.getOpenFileName(self, "Select Nifti file"))
            self.ui.lineEdit_ProcInfolder1.setText(newnifti)
            if not os.path.isfile(newnifti):
                print "File not found"
            if success:
                out = os.path.dirname(newnifti)
                self.ui.lineEdit_ProcOutputfolder.setText(out)

            # Reset checkboxes
            # self.ui.checkBox_gaussian3D.setChecked(False)
            # self.ui.checkBox_gaussian2D.setChecked(False)
            # self.ui.checkBox_median.setChecked(False)
            # self.ui.checkBox_wiener.setChecked(False)
            # self.UpdateGUI()
            logging.info('ChangeProc1file ' + newnifti)
        except ValueError:
            logging.info('ChangeProc1file error')

    def ChangeProcFile2(self):
        """ChangeProc2path
        """
        success = 1
        try:
            newnifti = str(
                QFileDialog.getOpenFileName(self, "Select Nifti file"))
            self.ui.lineEdit_ProcInfolder2.setText(newnifti)
            if not os.path.isfile(newnifti):
                print "File not found"
            logging.info('ChangeProc2file ' + newdir)
        except ValueError:
            logging.info('ChangeProc2file error')

    def ChangeCoilSensFile(self):
        """ChangeProc1path
        """
        success = 1
        try:
            newnifti = str(
                QFileDialog.getOpenFileName(self, "Select Nifti file"))
            self.ui.lineEdit_CoilSens.setText(newnifti)
            if not os.path.isfile(newnifti):
                print "File not found"
            if success:
                out = os.path.dirname(newnifti)
                self.ui.lineEdit_CoilSens.setText(out)

            logging.info('ChangeCoilSensFile ' + newnifti)
        except ValueError:
            logging.info('ChangeCoilSensFile error')

    def ClearProcLine1(self):
        """ChangeProc1path
        """
        self.ui.lineEdit_ProcInfolder1.setText("")

    def ClearProcLine2(self):
        """ChangeProc2path
        """
        self.ui.lineEdit_ProcInfolder2.setText("")

    def ClearProcLine3(self):
        """ChangeProc2path
        """
        self.ui.lineEdit_CoilSens.setText("")

    def ChangeProcOutputPath(self):
        """ChangeProcOutputPath
        """
        try:
            newdir = str(QFileDialog.getExistingDirectory(self,
                                                          '''Select or
                                                          create new
                                                          Nifti
                                                          Directory
                                                          '''))
            files = os.listdir(newdir)
            niifiles = [f for f in files if f.endswith('.nii.gz')]
            if not len(niifiles) == 0:
                quit_msg = "Do you want to delete existing niftis?"
                reply = QtGui.QMessageBox.question(self, 'Message',
                                                   quit_msg,
                                                   QtGui.QMessageBox.Yes,
                                                   QtGui.QMessageBox.No)

                if reply == QtGui.QMessageBox.Yes:
                    print subprocess.Popen('rm -f ' + newdir + '/*.nii.gz',
                                           stdout=subprocess.PIPE,
                                           shell=True,
                                           executable="/bin/bash").stdout.read()
                    logging.info('ProcOutputpath ... deleted Niftis.')

            self.ui.lineEdit_ProcOutputfolder.setText(newdir)
            self.UpdateGUI()
            logging.info('ChangeProcOutputpath ' + newdir)
        except ValueError:
            logging.info('ChangeProcOutputpath error')

    def ProcSWI(self):
        """Process SWI
        """
        try:
            input1_dir = str(self.ui.lineEdit_ProcInfolder1.text())
            input2_dir = str(self.ui.lineEdit_ProcInfolder2.text())
            output_dir = str(self.ui.lineEdit_ProcOutputfolder.text())
            preprocess = int(self.ui.checkBox_swipreprocess.isChecked())
            swiorder = str(self.ui.lineEdit_swiorder.text())
            saveRI = int(self.ui.checkBox_swisaveRI.isChecked())
            swineg = int(self.ui.checkBox_swineg.isChecked())
            swipos = int(self.ui.checkBox_swipos.isChecked())
            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            # print 'swifid path: %s' % thispath
            # cmd1 = os.path.join(thispath, 'swifid.sh')
            # cmd1 = cmd1 + ' -v '
            # cmd1 = cmd1 + ' -i ' + str(input1_dir) + ' -o ' + str(output_dir)
            # print cmd1
            # cmd = cmd_header + cmd1
            cmd = """%s matlab -nodesktop -nosplash -r "addpath %s/matlab;call_swi('"'"'%s'"'"','"'"'%s'"'"','"'"'%s'"'"',%s,%d,%d,%d,%d);quit"  """ % ( proc_header, str(thispath),str(input1_dir), str(input2_dir), str(output_dir), swiorder, preprocess, saveRI, swineg, swipos)

            print cmd
            logging.info('Processing SWI ' + cmd)
        except ValueError:
            logging.error('ProcSWI error', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        try:
            print subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True,
                                   executable="/bin/bash").stdout.read()
            logging.info('ProcSWI ... completed.')
        except ValueError:
            logging.error(
                'ProcSWI invalid arguments to Popen', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        except OSError:
            logging.error('ProcSWI OSError', exc_info=True)
            pass
        self.UpdateGUI()
        print 'ProcSWI completed.'
        logging.info('ProcSWI done')

    def ProcMCI(self):
        """Process MCI
        """
        try:
            input1_dir = str(self.ui.lineEdit_ProcInfolder1.text())
            input2_dir = str(self.ui.lineEdit_ProcInfolder2.text())
            output_dir = str(self.ui.lineEdit_ProcOutputfolder.text())
            saveRI = int(self.ui.radioButton_MCI_RI.isChecked())

            cmd = """%s matlab -nodesktop -nosplash -r "addpath ./matlab;call_mci('"'"'%s'"'"','"'"'%s'"'"','"'"'%s'"'"',%d);quit"  """ % (
                proc_header, str(input1_dir), str(input2_dir), str(output_dir), saveRI)
            # print cmd1
            # cmd = cmd_header + cmd1
            print cmd
            logging.info('Processing MCI ' + cmd)
        except ValueError:
            logging.error('ProcMCI error', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        try:
            print subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True,
                                   executable="/bin/bash").stdout.read()
            logging.info('ProcMCI ... completed.')
        except ValueError:
            logging.error(
                'ProcMCI invalid arguments to Popen', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        except OSError:
            logging.error('ProcMCI OSError', exc_info=True)
            pass
        self.UpdateGUI()
        print 'ProcMCI completed.'
        logging.info('ProcMCI done')

    def ProcMEE(self):
        """Process MEE
        """
        try:
            input1_dir = str(self.ui.lineEdit_ProcInfolder1.text())
            input2_dir = str(self.ui.lineEdit_ProcInfolder2.text())
            output_dir = str(self.ui.lineEdit_ProcOutputfolder.text())
            meeorder = str(self.ui.lineEdit_meeorder.text())
            preprocess = int(self.ui.checkBox_meepreproc.isChecked())
            saveRI = int(self.ui.checkBox_meesaveRI.isChecked())
            useswi = int(self.ui.checkBox_meeswi.isChecked())
            if useswi and saveRI:
                msg = "SWI is not complex, so real/imag output not possible. Do you want the SWI disabled and save the real/imag from the homodyne filtered image?"
                reply = QtGui.QMessageBox.question(
                    self, 'Message', msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
                if reply == QtGui.QMessageBox.Yes:
                    #                    event.accept()
                    self.ui.checkBox_meeswi.setChecked(False)
                else:
                    self.checkBox_meesaveRI.setChecked(False)
                    #                    event.ignore()

            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            cmd = """%s matlab -nodesktop -nosplash -r "addpath ./matlab;call_swi('"'"'%s'"'"','"'"'%s'"'"','"'"'%s'"'"',%s,%d,%d,%d);quit"  """ % (
                proc_header, str(input1_dir), str(input2_dir), os.path.join(str(output_dir), 'swi.nii.gz'), meeorder, preprocess, saveRI, useswi)

            print cmd
            logging.info('Processing MEE ' + cmd)
        except ValueError:
            logging.error('ProcMEE error', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        try:
            print subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True,
                                   executable="/bin/bash").stdout.read()
            logging.info('ProcMEE ... completed.')
        except ValueError:
            logging.error(
                'ProcMEE invalid arguments to Popen', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        except OSError:
            logging.error('ProcMEE OSError', exc_info=True)
            pass
        self.UpdateGUI()
        print 'ProcMEE completed.'
        logging.info('ProcMEE done')

    def GetNLfilter(self):
        """GetNLfilter select which NL filter is to be used
        """
        if self.ui.MRONLM.isChecked():
            return 0
        elif self.ui.PRINLM.isChecked():
            return 1
        elif self.ui.AONLM.isChecked():
            return 2
        elif self.ui.ONLM.isChecked():
            return 3
        elif self.ui.ODCT.isChecked():
            return 4
        elif self.ui.MRONLM2.isChecked():
            return -5
        elif self.ui.PRINLM2.isChecked():
            return -1
        elif self.ui.AONLM2.isChecked():
            return -2
        elif self.ui.ONLM2.isChecked():
            return -3
        elif self.ui.ODCT2.isChecked():
            return -4
        else:
            return 0

    def GetNLparams(self):
        sigmatext = str(self.ui.lineEdit_NLsigma.text())
        print sigmatext
        if sigmatext and re.search('[0-9]', sigmatext[0]):
            sigma = float(sigmatext)
        else:
            sigma = '[]'

        sigmaratiotext = str(self.ui.linedit_NLsigmaratio.text())
        print sigmaratiotext
        if sigmaratiotext and re.search('[0-9]', sigmaratiotext[0]):
            sigmaratio = float(sigmaratiotext)
        else:
            sigmaratio = '[]'
        searchtext = str(self.ui.lineEdit_NLsearcharea.text())
        if searchtext and re.search('[0-9]', searchtext[0]):
            searcharea = float(searchtext)
        else:
            searcharea = '[]'
        patchtext = str(self.ui.lineEdit_NLpatcharea.text())
        if patchtext and re.search('[0-9]', patchtext[0]):
            patcharea = float(patchtext)
        else:
            patcharea = '[]'
        rician = 1
        if not self.ui.checkBox_NLrician.isChecked():
            rician = 0
        return (sigma, sigmaratio, searcharea, patcharea, rician)

    def ProcNLpipeline1(self):
        """Process NLpipeline1
        """
        try:
            inputRI=0
            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            # print 'this path: %s' % thispath
            input1_dir = str(self.ui.lineEdit_ProcInfolder1.text())
            input2_dir = str(self.ui.lineEdit_ProcInfolder2.text())
            inputRI = int(self.ui.checkBox_NLinputRI.isChecked())
          #  if self.ui.checkBox_NLinputMagPha.isChecked():
          #      inputRI=inputRI+4
            output_dir = str(self.ui.lineEdit_ProcOutputfolder.text())
            NLfilter = self.GetNLfilter()

            if not input2_dir or inputRI:
                input2_dir = '[]'
                inputRI = 0

            if self.ui.checkBox_NLenableparams.isChecked():
                (sigma, sigmaratio, searcharea,
                 patcharea, rician) = self.GetNLparams()
                cmd = """ %s matlab -nodesktop -nosplash -r \"addpath %s/matlab/NLmeans; call_pipeline1('\"'\"'%s'\"'\"','\"'\"'%s'\"'\"','\"'\"'%s'\"'\"',%d,%s,%s,%s,%s,%s,%s);quit\" """ % (
                    proc_header, thispath, str(input1_dir), output_dir, str(input2_dir), NLfilter, str(inputRI),str(sigma), str(sigmaratio), str(searcharea), str(patcharea), str(rician))
            else:
                cmd = """ %s matlab -nodesktop -nosplash -r \"addpath %s/matlab/NLmeans; call_pipeline1('\"'\"'%s'\"'\"','\"'\"'%s'\"'\"','\"'\"'%s'\"'\"',%d,%s,[],[],[],[],1);quit\" """ % (
                    proc_header, thispath, str(input1_dir), output_dir, str(input2_dir), NLfilter, str(inputRI))
            print cmd

            # print cmd
            logging.info('Processing NLpipeline1 ' + cmd)
        except ValueError:
            logging.error('ProcNLpipeline1 setup error', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        try:
            print subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True,
                                   executable="/bin/bash").stdout.read()
            logging.info('ProcNLpipeline1 ... completed.')
        except ValueError:
            logging.error(
                'ProcNLpipeline1 invalid arguments to Popen', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        except OSError:
            logging.error('ProcNLpipeline1 OSError', exc_info=True)
            pass
        self.UpdateGUI()
        print 'ProcNLpipeline1 completed.'
        logging.info('ProcNLpipeline1 done')

    def ProcNLpipeline2(self):
        """Process NLpipeline2
        """
        try:
            input1_dir = str(self.ui.lineEdit_ProcInfolder1.text())
            input2_dir = str(self.ui.lineEdit_ProcInfolder2.text())
            if input2_dir=='':
                pp2_msg = "Pipeline 2 requires two input images. Click Ok to continue. "
                reply = QtGui.QMessageBox.question(self, 'Message', pp2_msg,
                                                   QtGui.QMessageBox.Ok,
                                                   QtGui.QMessageBox.Cancel)
                return
            output_dir = str(self.ui.lineEdit_ProcOutputfolder.text())
            NLfilter = self.GetNLfilter()

            thispath =os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            print 'this path: %s' % thispath
            if self.ui.checkBox_NLenableparams.isChecked():
                (sigma, sigmaratio, searcharea,
                 patcharea, rician) = self.GetNLparams()

                cmd = """ %s matlab -nodesktop -nosplash -r \"addpath %s/matlab/NLmeans;
                call_pipeline2('\"'\"'%s'\"'\"','\"'\"'%s'\"'\"','\"'\"'%s'\"'\"',%d,%s,%s,%s,%s,%s);quit\" """ % (
                    proc_header, thispath, str(input1_dir), str(input2_dir), str(output_dir),
                    NLfilter, str(sigma), str(sigmaratio), str(searcharea), str(patcharea), str(rician))
            else:
                cmd = """ %s matlab -nodesktop -nosplash -r \"addpath %s/matlab/NLmeans;
                call_pipeline2('\"'\"'%s'\"'\"','\"'\"'%s'\"'\"','\"'\"'%s'\"'\"',%d);quit\" """ % (
                    proc_header, thispath, str(input1_dir), str(input2_dir), str(output_dir), NLfilter)

            print cmd
            logging.info('Processing NLpipeline2 ' + cmd)
        except ValueError:
            logging.error('ProcNLpipeline2 error', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        try:
            print subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True,
                                   executable="/bin/bash").stdout.read()
            logging.info('ProcNLpipeline2 ... completed.')
        except ValueError:
            logging.error(
                'ProcNLpipeline2 invalid arguments to Popen', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        except OSError:
            logging.error('ProcNLpipeline2 OSError', exc_info=True)
            pass
        self.UpdateGUI()
        print 'ProcNLpipeline2 completed.'
        logging.info('ProcNLpipeline2 done')

    def ProcNLpipeline3(self):
        """Process NLpipeline3
        """
        try:
            input1_dir = str(self.ui.lineEdit_ProcInfolder1.text())
            input2_dir = str(self.ui.lineEdit_ProcInfolder2.text())
            output_dir = str(self.ui.lineEdit_ProcOutputfolder.text())
            saveRI = int(self.ui.checkBox_NLsaveRI.isChecked())
            savePhase = int(self.ui.checkBox_NLsavePhase.isChecked())
            saveSWI = int(self.ui.checkBox_SaveSWIProc3.isChecked())
            NLfilter = self.GetNLfilter()
            thispath = os.path.dirname(
                os.path.realpath(os.path.abspath(__file__)))
            print 'this path: %s' % thispath
            if self.ui.checkBox_NLenableparams.isChecked():
                (sigma, sigmaratio, searcharea,
                 patcharea, rician) = self.GetNLparams()
                if rician == 1:
                    pp3_msg = """Pipeline 3 uses complex images and the rician variable should ideally be set to zero.
                    Click Yes to continue with Rician=1.
                    Click No to set Rician to zero then continue pipeline.
                    Click Cancel to stop pipeline. """
                    reply = QtGui.QMessageBox.question(self, 'Message', pp3_msg,
                                                       QtGui.QMessageBox.Yes,
                                                       QtGui.QMessageBox.No, QtGui.QMessageBox.Cancel)
                    if reply == QtGui.QMessageBox.No:
                    #                    event.accept()
                        rician=0
                        self.ui.checkBox_NLmeansrician.setChecked(False)
                    elif reply == QtGui.QMessageBox.Cancel:
                        return
                cmd = """ %s matlab -nodesktop -nosplash -r \"addpath %s/matlab/NLmeans;
                call_pipeline3('\"'\"'%s'\"'\"','\"'\"'%s'\"'\"','\"'\"'%s'\"'\"',%d,%d,%s,%s,%s,%s,%s);quit\" """ % (
                    proc_header, thispath, str(input1_dir), str(input2_dir),
                    str(output_dir), saveRI + 4 * savePhase +
                    8 * saveSWI, NLfilter, str(sigma),
                    str(sigmaratio), str(searcharea), str(patcharea), str(rician))
            else:
                cmd = """ %s matlab -nodesktop -nosplash -r \"addpath %s/matlab/NLmeans;
                call_pipeline3('\"'\"'%s'\"'\"','\"'\"'%s'\"'\"','\"'\"'%s'\"'\"',%s,%d);quit\" """ % (
                    proc_header, thispath, str(input1_dir), str(input2_dir),
                    str(output_dir), str(saveRI + 4 * savePhase + 8 * saveSWI), NLfilter)

            print cmd
            logging.info('Processing NLpipeline3 ' + cmd)
        except ValueError:
            logging.error('ProcNLpipeline3 error', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        try:
            print subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   shell=True,
                                   executable="/bin/bash").stdout.read()
            logging.info('ProcNLpipeline3 ... completed.')
        except ValueError:
            logging.error(
                'ProcNLpipeline3 invalid arguments to Popen', exc_info=True)
            # this is necessary to avoid the parent try from catching this
            # exception
            pass
        except OSError:
            logging.error('ProcNLpipeline3 OSError', exc_info=True)
            pass
        self.UpdateGUI()
        print 'ProcNLpipeline3 completed.'
        logging.info('ProcNLpipeline3 done')


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main = Agilent2DicomWindow()
    main.show()
    sys.exit(app.exec_())

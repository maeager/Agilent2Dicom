import os,sys,re
from PyQt4 import Qt, QtGui, QtCore
from PyQt4.QtGui import QDialog,QFileDialog,QApplication
from Agilent2DicomWidget import Ui_Form
from fdf2dcm_global import *
DEBUGGING=1

cmd_header='(if test ${MASSIVEUSERNAME+defined} \n\
then \n\
echo ''On Massive'' \n\
module unload python \n\
module load python/2.7.1-gcc \n\
module load python/2.7.3-gcc dcmtk mrtrix \n\
module list \n\
fi \n\
echo ''Done''\n '

class Agilent2DicomWindow(QtGui.QWidget):
    def __init__(self):
        super(Agilent2DicomWindow,self).__init__()
        self.ui=Ui_Form()
    # Set up the user interface from Designer.
        self.ui.setupUi(self)
        
        # Make some local modifications.
        # self.colorDepthCombo.addItem("2 colors (1 bit per pixel)")
        
        # Disable some features
        if DEBUGGING == 1:
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
        try:
            newdir = str(QFileDialog.getExistingDirectory(self, "Select FDF Directory"))
            self.ui.lineEdit_fdfpath.setText(newdir)
	    if re.search('img',newdir):
              out = re.sub('img','dcm',newdir)
	    else:
	      out = dir_+'.dcm'
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
            newdir = file = str(QFileDialog.getOpenFile(self, "Select DICOM Directory"))
            self.ui.lineEdit_dicompath.setText(newdir)
            self.ui.lineEdit_darisid.setText(self.GetDarisID(out))
            self.UpdateGUI()
        except ValueError:
            pass
		  
    def ChangeFIDpath(self):
        try:
            newdir = str(QFileDialog.getExistingDirectory(self, "Select FID Directory"))
            self.ui.lineEdit_fidpath.setText(newdir)
            if re.search('img',newdir):
                out = re.sub('img','dcm',newdir)
            else:
                out = dir_+'.dcm'
            self.ui.lineEdit_dicompath2.setText(out)
            self.ui.lineEdit_darisid2.setText(self.GetDarisID(newdir))
            self.UpdateGUI()
        except ValueError:
            pass
			  
    def ChangeFIDDicomPath(self):
        try:
            newdir = str(QFileDialog.getOpenFile(self, "Select DICOM Directory"))
            self.ui.lineEdit_dicompath2.setText(newdir)
            self.ui.lineEdit_darisid.setText(self.GetDarisID(out))
            self.UpdateGUI()
        except ValueError:
            pass
			      
    def ConvertFDF(self):
        try:
            import os
            input_dir = self.ui.lineEdit_fdfpath.getText()
            output_dir = self.ui.lineEdit_dicompath.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'fdf2dcm path: %s' % thispath
            cmd1 = os.path.join(thispath,'fdf2dcm.sh') + ' -v  -i ' + input_dir + ' -o ' + output_dir
            print(cmd1)
            cmd= cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
				  
    def GetDarisID(inpath):
        daris_id=''
        try:
            procpar,pep = ReadProcpar.ReadProcpar(os.path.join(inpath,'procpar'))
            if 'name' in procpar.keys():
                if re.search('DaRIS',procpar['name']):
                    daris_id = re.sub('DaRIS\^','',procpar['name'])
        except ValueError:
            pass
        return daris_id

    def CheckDicmoDir(dpath):
        if dpath and os.path.isdir(dpath) and os.path.exists(os.path.join(dpath,'0001.dcm')):
            return True
        else:
            return False
					  					  
    def Send2Daris(*args):
        try:
            daris_ID = self.ui.lineEdit_darisid.getText()
            dicom_dir = self.ui.lineEdit_dicompath.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd = os.path.join(thispath,'dpush') + ' -c ' + daris_ID + ' -s mf-erc ' + dicom_dir
            os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
					      
    def CheckFDF(*args):  #send_button):
        try:
            
            output_dir = self.ui.lineEdit_dicompath.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'check path: %s' % thispath
            cmd1 = os.path.join(thispath,'dcheck.sh') + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
						  
    def ViewFDF(self):
        try:
            output_dir = self.ui.lineEdit_dicompath.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd1 = os.path.join(thispath,'mrview') + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
        except ValueError:
            pass
						      
    def getCommandArgs(self):
        argstr=''
        if self.ui.checkBox_magn.checked:
            argstr=' -m'
        if self.ui.checkBox_pha.checked:
            argstr+=' -p'
        if self.ui.checkBox_ksp.checked:
            argstr+=' -k'
        if self.ui.checkBox_reimag.checked:
            argstr+=' -r'
        
        if self.ui.checkBox_gaussian:
            argstr+=' -r -s %s -go %s' % (self.ui.lineEdit_gsigma,self.ui.lineEdit_gorder)
        if self.ui.nearest.checked:
            argstr+=' -gm nearest'
        elif self.ui.reflect.checked:
            argstr+=' -gm reflect'
        elif self.ui.wrap.checked:
            argstr+=' -gm wrap'
        if self.ui.checkBox_median:
            argstr+=' -d -n %s ' % (self.ui.lineEdit_median_size)
        if self.ui.checkBox_wiener:
            argstr+=' -w -ws %s -wn %s' % (self.ui.lineEdit_wiener_size, self.ui.lineEdit_wiener_noise)
									    
									    
									    
									    
    def ConvertFID(self):
        try:
            input_dir = self.ui.lineEdit_fidpath.getText()
            output_dir = self.ui.lineEdit_dicompath2.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'fid2dcm path: %s' % thispath
            cmd1 = os.path.join(thispath,'fid2dcm.sh')
            cmd1 = cmd1+' -v '  
            cmd1 = cmd+'-i ' + input_dir + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
										
    def CheckFID(*args):  #send_button):
        try:
            
            output_dir = self.ui.lineEdit_dicompath2.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'check path: %s' % thispath
            cmd1 = os.path.join(thispath,'dcheck.sh') + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            self.UpdateGUI()
        except ValueError:
            pass
										    
    def ViewFID(self):
        try:
            output_dir = self.ui.lineEdit_dicompath2.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd1 = os.path.join(thispath,'mrview') + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            UpdateGUI()
        except ValueError:
            pass
											
    def Send2Daris2(*args):
        try:
            daris_ID = self.ui.lineEdit_darisid2.getText()
            dicom_dir = self.ui.lineEdit_dicompath2.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd = os.path.join(thispath,'dpush') + ' -c ' + daris_ID + ' -s mf-erc ' + dicom_dir
            os.system(cmd)	
        except ValueError:
            pass
											    
    def UpdateGUI(self):
        self.ui.pushButton_check.setEnabled(False)
        self.ui.pushButton_view.setEnabled(False)
        self.ui.pushButton_send2daris.setEnabled(False)
        self.ui.pushButton_check2.setEnabled(False)
        self.ui.pushButton_view2.setEnabled(False)
        self.ui.pushButton_send2daris2.setEnabled(False)
        if CheckDicmoDir(self.ui.lineEdit_dicompath.getText()):
            self.ui.pushButton_check.setEnabled(True)
            self.ui.pushButton_view.setEnabled(True)
            self.ui.pushButton_send2daris.setEnabled(True)
        if CheckDicmoDir(self.ui.lineEdit_dicompath2.getText()):
            self.ui.pushButton_check2.setEnabled(True)
            self.ui.pushButton_view2.setEnabled(True)
            self.ui.pushButton_send2daris2.setEnabled(True)
												  


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main = Agilent2DicomWindow()
    main.show()
    sys.exit(app.exec_())
												  

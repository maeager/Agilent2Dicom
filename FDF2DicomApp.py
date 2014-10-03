import os,sys
from PyQt4 import Qt
from PyQt4.QtGui import QDialog,QFileDialog,QApplication
from FDF2DicomQt import Ui_Dialog
from fdf2dcm_global import *

cmd_header='(if test ${MASSIVEUSERNAME+defined} \n\
then \n\
  echo ''On Massive'' \n\
  module unload python \n\
  module load python/2.7.1-gcc \n\
  module load python/2.7.3-gcc dcmtk mrtrix \n\
  module list \n\
fi \n\
echo ''Done''\n '

class ImageConverterDialog(QDialog, Ui_Dialog):
    def __init__(self):
        QDialog.__init__(self)

        # Set up the user interface from Designer.
        self.setupUi(self)

        # Make some local modifications.
        # self.colorDepthCombo.addItem("2 colors (1 bit per pixel)")

	# Disable some features
	#self.tab_diffusion.setEnabled(False)
	#self.tab_multiecho.setEnabled(False)
	#self.pushButton_check.setEnabled(False)
	#self.pushButton_view.setEnabled(False)
	#self.pushButton_send2daris.setEnabled(False)
	#self.pushButton_check2.setEnabled(False)
	#self.pushButton_view2.setEnabled(False)
	#self.pushButton_send2daris2.setEnabled(False)
	#self.checkBox_median.setChecked(False)
	#self.checkBox_median.setEnabled(False)
	#self.lineEdit_median_size.setEnabled(False)
	#self.checkBox_wiener.setEnabled(False)
	#self.lineEdit_wiener_size.setEnabled(False)
	#self.lineEdit_wienernoise.setEnabled(False)
	#self.checkBox_magn.setEnabled(True)
	#self.checkBox_magn.setChecked(True)
	#self.checkBox_ksp.setChecked(False)
	#self.checkBox_reimag.setChecked(False)
	#self.checkBox_pha.setChecked(False)


        # Connect up the buttons.
        self.connect(self.buttonBox, Qt.SIGNAL("accepted()"), self.accept)
        self.connect(self.buttonBox ,Qt.SIGNAL("rejected()"), self.reject)
        #self.buttonBox.Ok.clicked.connect(self.accept)
        #self.buttonBox.cancelButton.clicked.connect(self.reject)
        self.pushButton_changefdf.clicked.connect(self.ChangeFDFpath)
        self.pushButton_changedicom.clicked.connect(self.ChangeFDFDicomPath)
        self.pushButton_changefid.clicked.connect(self.ChangeFIDpath)
        self.pushButton_changedicom2.clicked.connect(self.ChangeFIDDicomPath)
        self.pushButton_convert.clicked.connect(self.ConvertFDF)
        self.pushButton_check.clicked.connect(self.CheckFDF)
        self.pushButton_view.clicked.connect(self.ViewFDF)
        self.pushButton_send2daris.clicked.connect(self.Send2Daris)
	self.pushButton_convertfid.clicked.connect(self.ConvertFID)
        self.pushButton_check2.clicked.connect(self.CheckFID)
        self.pushButton_view2.clicked.connect(self.ViewFID)
        self.pushButton_send2daris2.clicked.connect(self.Send2Daris2)

    def ChangeFDFpath(self):
      try:
        newdir = file = str(QFileDialog.getExistingDirectory(self, "Select FDF Directory"))
        self.lineEdit_fdfpath.setText(newdir)
        if re.search('img',newdir):
            out = re.sub('img','dcm',newdir)
        else:
            out = dir_+'.dcm'
        self.lineEdit_dicompath.setText(out)
        self.lineEdit_darisid.setText(self.GetDarisID(newdir))
        UpdateGUI()
     except ValueError:
        pass

    def ChangeFDFDicomPath(self):
      try:
        newdir = file = str(QFileDialog.getOpenFile(self, "Select DICOM Directory"))
        self.lineEdit_dicompath.setText(newdir)
        self.lineEdit_darisid.setText(self.GetDarisID(out))
        UpdateGUI()
     except ValueError:
        pass

    def ChangeFIDpath(self):
      try:
	newdir = str(QFileDialog.getExistingDirectory(self, "Select FID Directory"))
        self.lineEdit_fidpath.setText(newdir)
        if re.search('img',newdir):
            out = re.sub('img','dcm',newdir)
        else:
            out = dir_+'.dcm'
        self.lineEdit_dicompath2.setText(out)
        self.lineEdit_darisid2.setText(self.GetDarisID(newdir))
        UpdateGUI()
      except ValueError:
        pass

    def ChangeFIDDicomPath(self):
      try:
	newdir = file = str(QFileDialog.getOpenFile(self, "Select DICOM Directory"))
        self.lineEdit_dicompath2.setText(newdir)
        self.lineEdit_darisid.setText(self.GetDarisID(out))
	UpdateGUI()
      except ValueError:
        pass

    def ConvertFDF(self):
      try:
        import os
        input_dir = self.lineEdit_fdfpath.getText()
        output_dir = self.lineEdit_dicompath.getText()
        thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
        print 'fdf2dcm path: %s' % thispath
        cmd1 = os.path.join(thispath,'fdf2dcm.sh') + ' -v  -i ' + input_dir + ' -o ' + output_dir
        print(cmd1)
        cmd= cmd_header + cmd1 +')'
        print(cmd)
        os.system(cmd)
        UpdateGUI()
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
            daris_ID = self.lineEdit_darisid.getText()
            dicom_dir = self.lineEdit_dicompath.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd = os.path.join(thispath,'dpush') + ' -c ' + daris_ID + ' -s mf-erc ' + dicom_dir
            os.system(cmd)
            UpdateGUI()
        except ValueError:
            pass

    def CheckFDF(*args):  #send_button):
        try:
            import os
            output_dir = self.lineEdit_dicompath.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'check path: %s' % thispath
            cmd1 = os.path.join(thispath,'dcheck.sh') + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            UpdateGUI()
        except ValueError:
            pass

    def ViewFDF(self):
	try:
	    output_dir = self.lineEdit_dicompath.getText()
	    thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd1 = os.path.join(thispath,'mrview') + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
        except ValueError:
	    pass
      
      
    def ConvertFID(self):
        try:
            import os
            input_dir = self.lineEdit_fidpath.getText()
            output_dir = self.lineEdit_dicompath2.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'fid2dcm path: %s' % thispath
            cmd1 = os.path.join(thispath,'fid2dcm.sh') + ' -v  -i ' + input_dir + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
            UpdateGUI()
        except ValueError:
            pass

    def CheckFID(*args):  #send_button):
        try:
            import os
            output_dir = self.lineEdit_dicompath2.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'check path: %s' % thispath
            cmd1 = os.path.join(thispath,'dcheck.sh') + ' -o ' + output_dir
            print(cmd1)
            cmd=cmd_header + cmd1 +')'
            print(cmd)
            os.system(cmd)
	    UpdateGUI()
	except ValueError:
            pass

    def ViewFID(self):
	try:
	    output_dir = self.lineEdit_dicompath2.getText()
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
            daris_ID = self.lineEdit_darisid2.getText()
            dicom_dir = self.lineEdit_dicompath2.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            cmd = os.path.join(thispath,'dpush') + ' -c ' + daris_ID + ' -s mf-erc ' + dicom_dir
            os.system(cmd)	
        except ValueError:
            pass

    def UpdateGUI(self):
      	self.pushButton_check.setEnabled(False)
	self.pushButton_view.setEnabled(False)
	self.pushButton_send2daris.setEnabled(False)
	self.pushButton_check2.setEnabled(False)
	self.pushButton_view2.setEnabled(False)
	self.pushButton_send2daris2.setEnabled(False)
	if CheckDicmoDir(self.lineEdit_dicompath.getText()):
	  self.pushButton_check.setEnabled(True)
	  self.pushButton_view.setEnabled(True)
	  self.pushButton_send2daris.setEnabled(True)
	if CheckDicmoDir(self.lineEdit_dicompath2.getText()):
	  self.pushButton_check2.setEnabled(True)
	  self.pushButton_view2.setEnabled(True)
	  self.pushButton_send2daris2.setEnabled(True)
	  

app = QApplication(sys.argv)
window = QDialog()
ui = ImageConverterDialog()
ui.setupUi(window)

window.show()
sys.exit(app.exec_())

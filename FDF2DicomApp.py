import os,sys
from PyQt4 import Qt
from PyQt4.QtGui import QDialog,QFileDialog,QApplication
from FDF2DicomQt import Ui_Dialog

class ImageConverterDialog(QDialog, Ui_Dialog):
    def __init__(self):
        QDialog.__init__(self)

        # Set up the user interface from Designer.
        self.setupUi(self)

        # Make some local modifications.
        # self.colorDepthCombo.addItem("2 colors (1 bit per pixel)")

        # Connect up the buttons.
        self.connect(self.buttonBox, Qt.SIGNAL("accepted()"), self.accept)
        self.connect(self.buttonBox ,Qt.SIGNAL("rejected()"), self.reject)
        #self.buttonBox.Ok.clicked.connect(self.accept)
        #self.buttonBox.cancelButton.clicked.connect(self.reject)
        self.pushButton_changefdf.clicked.connect(self.ChangeFDFpath)
        self.pushButton_changedicom.clicked.connect(self.ChangeFDFDicomPath)
        self.pushButton_changefid.clicked.connect(self.ChangeFIDpath)
        self.pushButton_changedicom2.clicked.connect(self.ChangeFIDDicomPath)
        self.pushButton_send2daris.clicked.connect(self.Send2Daris)

    def ChangeFDFpath(self):
        newdir = file = str(QFileDialog.getExistingDirectory(self, "Select FDF Directory"))
        self.lineEdit_fdfpath.setText(newdir)
        if re.search('img',newdir):
            out = re.sub('img','dcm',newdir)
        else:
            out = dir_+'.dcm'
        self.lineEdit_dicompath.setText(out)
        self.lineEdit_darisid.setText(self.GetDarisID(newdir))

    def ChangeFDFDicomPath(self):
        newdir = file = str(QFileDialog.getOpenFile(self, "Select DICOM Directory"))
        self.lineEdit_dicompath.setText(newdir)
        self.lineEdit_darisid.setText(self.GetDarisID(out))

    def ChangeFIDpath(self):
        newdir = str(QFileDialog.getExistingDirectory(self, "Select FID Directory"))
        self.lineEdit_fidpath.setText(newdir)
        if re.search('img',newdir):
            out = re.sub('img','dcm',newdir)
        else:
            out = dir_+'.dcm'
        self.lineEdit_dicompath2.setText(out)
        self.lineEdit_darisid2.setText(self.GetDarisID(newdir))

    def ChangeFIDDicomPath(self):
        newdir = file = str(QFileDialog.getOpenFile(self, "Select DICOM Directory"))
        self.lineEdit_dicompath2.setText(newdir)
        self.lineEdit_darisid.setText(self.GetDarisID(out))

    def ConvertFDF(self):
      try:
        import os
        input_dir = self.lineEdit_fdfpath.getText()
        output_dir = self.lineEdit_dicompath.get()
        thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
        print 'dconvert path: %s' % thispath
        cmd1 = os.path.join(thispath,'fdf2dcm.sh') + ' -v  -i ' + input_dir + ' -o ' + output_dir
        print(cmd1)
        cmd='(if test ${MASSIVEUSERNAME+defined} \n\
then \n\
  echo ''On Massive'' \n\
  module unload python/3.3.5-gcc \n\
  module load python \n\
fi \n\
echo ''Done''\n ' + cmd1 +')'
        print(cmd)
        os.system(cmd)
        if check_dir(output_dir):
	      print 'Ready to send dicoms to DaRIS'
 #           send_button.foreground="dark green"
 #           send_button.state='active'
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
        if os.path.isdir(dpath) and os.path.exists(os.path.join(dpath,'0001.dcm')):
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
        except ValueError:
            pass

    def Check(*args):  #send_button):
        try:
            import os
            output_dir = outputdir.get()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'check path: %s' % thispath
            cmd1 = os.path.join(thispath,'dcheck.sh') + ' -o ' + output_dir
            print(cmd1)
            cmd='(if test ${MASSIVEUSERNAME+defined} \n\
then \n\
  echo ''On Massive'' \n\
  module unload python/3.3.5-gcc \n\
  module load python \n\
fi \n\
echo ''Done''\n ' + cmd1 +')'
            print(cmd)
            os.system(cmd)
            if check_dir(output_dir):
                print 'Ready to send dicoms to DaRIS'
 #           send_button.foreground="dark green"
 #           send_button.state='active'
        except ValueError:
            pass

    def onButton_ConvertFID(self):
        try:
            import os
            input_dir = self.lineEdit_fidpath.getText()
            output_dir = self.lineEdit_dicompath2.getText()
            thispath = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
            print 'fid2dcm path: %s' % thispath
            cmd1 = os.path.join(thispath,'fid2dcm.sh') + ' -v  -i ' + input_dir + ' -o ' + output_dir
            print(cmd1)
            cmd='(if test ${MASSIVEUSERNAME+defined} \n\
then \n\
  echo ''On Massive'' \n\
  module unload python/3.3.5-gcc \n\
  module load python \n\
fi \n\
echo ''Done''\n ' + cmd1 +')'
            print(cmd)
            os.system(cmd)
            if check_dir(output_dir):
                print 'Ready to send dicoms to DaRIS'
 #           send_button.foreground="dark green"
 #           send_button.state='active'
        except ValueError:
            pass


app = QApplication(sys.argv)
window = QDialog()
ui = ImageConverterDialog()
ui.setupUi(window)

window.show()
sys.exit(app.exec_())

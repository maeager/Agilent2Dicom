# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dcmcleanup.ui'
#
# Created: Tue Feb 24 10:05:31 2015
#      by: PyQt4 UI code generator 4.11
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
import os

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8

    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)


class Ui_CleanUpDialog(object):
    dcmdirs = []
    numcheckboxes = 0

    def setupUi(self, Dialog, inputdir):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(400, 300)
        self.verticalLayout_2 = QtGui.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(Dialog)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        if not os.path.exists(inputdir):
            inputdir = str(
                os.path.dirname(os.path.realpath(os.path.abspath(__file__))))
        files = os.listdir(inputdir)
        print files
        self.dcmdirs = [
            f for f in files if f.endswith('.dcm') or f.endswith('.nii')]
        print self.dcmdirs
        indx = 0
        if self.dcmdirs:

            print "Creating check boxes"
            for outdir in self.dcmdirs:
                print "Creating check box ", outdir
                # if os.path.exists(outdir):
                #    print "Creating check box ", outdir
                try:
                    cmdline = "self.checkBox_" + \
                        str(indx) + " = QtGui.QCheckBox(Dialog);"
                    exec(cmdline)  # , "<string>", "exec")
                    cmdline = "self.checkBox_" + \
                        str(indx) + ".setObjectName(_fromUtf8(\"checkbox\"));"
                    exec(cmdline)  # , "<string>", "exec")
                    cmdline = "self.verticalLayout.addWidget(self.checkBox_" + str(
                        indx) + ")"
                    exec(cmdline)  # , "<string>", "exec")
                except SyntaxError:
                    print cmdline
                    print "Syntax Error in retranslateUi"

                print cmdline
                #compile(cmdline, "<string>", "exec")
                indx = indx + 1
        self.numcheckboxes = indx
        # self.checkBox_2 = QtGui.QCheckBox(Dialog)
        # self.checkBox_2.setObjectName(_fromUtf8("checkBox_2"))
        # self.verticalLayout.addWidget(self.checkBox_2)

        spacerItem = QtGui.QSpacerItem(
            20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(
            QtGui.QDialogButtonBox.Cancel | QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Help)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout_2.addWidget(self.buttonBox)

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(
            self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Dialog.accept)
        QtCore.QObject.connect(
            self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Dialog.reject)
        QtCore.QObject.connect(
            self.buttonBox, QtCore.SIGNAL(_fromUtf8("helpRequested()")), Dialog.help)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(
            _translate("Dialog", "Agilent2Dicom Clean up", None))
        self.label.setText(
            _translate("Dialog", "Delete existing DICOM folders", None))
        indx = 0
        for outdir in self.dcmdirs:

            try:
                if eval('self.checkBox_' + str(indx)):
                    cmdline = 'self.checkBox_' + \
                        str(indx) + \
                        '.setText(_translate("Dialog", str(outdir), None))'
                    print cmdline
                    exec(cmdline)  # , "<string>", "exec")
                else:
                    print "check box ", indx, "not valid"
            except SyntaxError:
                print "Syntax Error in retranslateUi"

            indx = indx + 1

    def getValues(self):
        """
        """

        cindx = 0
        checked_dcmdirs = []
        for indx in xrange(0, self.numcheckboxes):
            try:
                if eval("self.checkBox_" + str(indx) + ".isChecked()"):
                    print "checkBox ", indx, self.dcmdirs[indx], " checked"
                    checked_dcmdirs.append(self.dcmdirs[indx])
                    cindx = cindx + 1
#                exec(cmdline) #, "<string>", "exec")
            except SyntaxError:
                print "Syntax Error in getValues"

#            if os.path.exists(outdir):
#                exec("if self.checkBox_"+str(indx)+".isChecked(): checked_dcmdirs[cindx] = dcmdirs[indx]")

        return checked_dcmdirs

    def help(self):
        '''Execute the command in response to the OK button.'''
        for indx in xrange(0, self.numcheckboxes):
            try:
                if eval("self.checkBox_" + str(indx) + ".isChecked()"):
                    print "checkBox ", indx, self.dcmdirs[indx], " checked ", eval("self.checkBox_" + str(indx) + ".isChecked()")

#                exec(cmdline) #, "<string>", "exec")
            except SyntaxError:
                print "Syntax Error in getValues"
        self.close()

    def accept(self):
        '''Execute the command in response to the OK button.'''
        self.close()

    def reject(self):
        '''Cancel.'''
        self.close()

    def close(self):
        self.dcmdirs = []

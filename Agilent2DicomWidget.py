# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Agilent2DicomWidget.ui'
#
# Created: Tue Oct 28 13:58:01 2014
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

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


class Ui_Form(object):

    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(565, 471)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(Form)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.tabWidget = QtGui.QTabWidget(Form)
        self.tabWidget.setMouseTracking(True)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab_fdf = QtGui.QWidget()
        self.tab_fdf.setObjectName(_fromUtf8("tab_fdf"))
        self.verticalLayout_12 = QtGui.QVBoxLayout(self.tab_fdf)
        self.verticalLayout_12.setObjectName(_fromUtf8("verticalLayout_12"))
        self.gridLayout_4 = QtGui.QGridLayout()
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.label_12 = QtGui.QLabel(self.tab_fdf)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.gridLayout_4.addWidget(self.label_12, 3, 1, 1, 1)
        self.pushButton_changefdf = QtGui.QPushButton(self.tab_fdf)
        self.pushButton_changefdf.setObjectName(
            _fromUtf8("pushButton_changefdf"))
        self.gridLayout_4.addWidget(self.pushButton_changefdf, 0, 4, 1, 1)
        self.label_13 = QtGui.QLabel(self.tab_fdf)
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.gridLayout_4.addWidget(self.label_13, 0, 1, 1, 1)
        self.label_16 = QtGui.QLabel(self.tab_fdf)
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.gridLayout_4.addWidget(self.label_16, 1, 1, 1, 1)
        self.label_17 = QtGui.QLabel(self.tab_fdf)
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.gridLayout_4.addWidget(self.label_17, 2, 1, 1, 1)
        self.pushButton_changedicom = QtGui.QPushButton(self.tab_fdf)
        self.pushButton_changedicom.setObjectName(
            _fromUtf8("pushButton_changedicom"))
        self.gridLayout_4.addWidget(self.pushButton_changedicom, 1, 4, 1, 1)
        self.lineEdit_darisid = QtGui.QLineEdit(self.tab_fdf)
        self.lineEdit_darisid.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_darisid.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_darisid.setObjectName(_fromUtf8("lineEdit_darisid"))
        self.gridLayout_4.addWidget(self.lineEdit_darisid, 2, 3, 1, 1)
        self.lineEdit_fdfpath = QtGui.QLineEdit(self.tab_fdf)
        self.lineEdit_fdfpath.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_fdfpath.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_fdfpath.setObjectName(_fromUtf8("lineEdit_fdfpath"))
        self.gridLayout_4.addWidget(self.lineEdit_fdfpath, 0, 3, 1, 1)
        self.lineEdit_dicompath = QtGui.QLineEdit(self.tab_fdf)
        self.lineEdit_dicompath.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_dicompath.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_dicompath.setObjectName(_fromUtf8("lineEdit_dicompath"))
        self.gridLayout_4.addWidget(self.lineEdit_dicompath, 1, 3, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.gridLayout_4.addLayout(self.horizontalLayout, 5, 3, 1, 1)
        self.tab_fdfoptions = QtGui.QTabWidget(self.tab_fdf)
        self.tab_fdfoptions.setEnabled(True)
        self.tab_fdfoptions.setObjectName(_fromUtf8("tab_fdfoptions"))
        self.tab_generic = QtGui.QWidget()
        self.tab_generic.setObjectName(_fromUtf8("tab_generic"))
        self.verticalLayout_11 = QtGui.QVBoxLayout(self.tab_generic)
        self.verticalLayout_11.setObjectName(_fromUtf8("verticalLayout_11"))
        self.gridLayout_7 = QtGui.QGridLayout()
        self.gridLayout_7.setObjectName(_fromUtf8("gridLayout_7"))
        self.checkBox_debugging = QtGui.QCheckBox(self.tab_generic)
        self.checkBox_debugging.setObjectName(_fromUtf8("checkBox_debugging"))
        self.gridLayout_7.addWidget(self.checkBox_debugging, 0, 0, 1, 1)
        self.checkBox_nodcmulti = QtGui.QCheckBox(self.tab_generic)
        self.checkBox_nodcmulti.setObjectName(_fromUtf8("checkBox_nodcmulti"))
        self.gridLayout_7.addWidget(self.checkBox_nodcmulti, 1, 0, 1, 1)
        self.verticalLayout_11.addLayout(self.gridLayout_7)
        self.tab_fdfoptions.addTab(self.tab_generic, _fromUtf8(""))
        self.tab_diffusion = QtGui.QWidget()
        self.tab_diffusion.setEnabled(False)
        self.tab_diffusion.setObjectName(_fromUtf8("tab_diffusion"))
        self.verticalLayout_13 = QtGui.QVBoxLayout(self.tab_diffusion)
        self.verticalLayout_13.setObjectName(_fromUtf8("verticalLayout_13"))
        self.gridLayout_8 = QtGui.QGridLayout()
        self.gridLayout_8.setObjectName(_fromUtf8("gridLayout_8"))
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setSpacing(5)
        self.horizontalLayout_6.setContentsMargins(-1, 5, -1, 5)
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.radioButton_bmatrixdiff = QtGui.QRadioButton(self.tab_diffusion)
        self.radioButton_bmatrixdiff.setObjectName(
            _fromUtf8("radioButton_bmatrixdiff"))
        self.horizontalLayout_6.addWidget(self.radioButton_bmatrixdiff)
        self.radioButton_directionaldiff = QtGui.QRadioButton(
            self.tab_diffusion)
        self.radioButton_directionaldiff.setObjectName(
            _fromUtf8("radioButton_directionaldiff"))
        self.horizontalLayout_6.addWidget(self.radioButton_directionaldiff)
        self.gridLayout_8.addLayout(self.horizontalLayout_6, 7, 1, 1, 1)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setSpacing(5)
        self.horizontalLayout_7.setContentsMargins(-1, 5, -1, 5)
        self.horizontalLayout_7.setObjectName(_fromUtf8("horizontalLayout_7"))
        self.label_18 = QtGui.QLabel(self.tab_diffusion)
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.horizontalLayout_7.addWidget(self.label_18)
        self.lineEdit_8 = QtGui.QLineEdit(self.tab_diffusion)
        self.lineEdit_8.setObjectName(_fromUtf8("lineEdit_8"))
        self.horizontalLayout_7.addWidget(self.lineEdit_8)
        self.gridLayout_8.addLayout(self.horizontalLayout_7, 7, 0, 1, 1)
        self.checkBox_diffusion1 = QtGui.QCheckBox(self.tab_diffusion)
        self.checkBox_diffusion1.setObjectName(
            _fromUtf8("checkBox_diffusion1"))
        self.gridLayout_8.addWidget(self.checkBox_diffusion1, 0, 0, 1, 1)
        self.verticalLayout_13.addLayout(self.gridLayout_8)
        self.tab_fdfoptions.addTab(self.tab_diffusion, _fromUtf8(""))
        self.tab_multiecho = QtGui.QWidget()
        self.tab_multiecho.setEnabled(False)
        self.tab_multiecho.setObjectName(_fromUtf8("tab_multiecho"))
        self.verticalLayout_14 = QtGui.QVBoxLayout(self.tab_multiecho)
        self.verticalLayout_14.setObjectName(_fromUtf8("verticalLayout_14"))
        self.gridLayout_10 = QtGui.QGridLayout()
        self.gridLayout_10.setObjectName(_fromUtf8("gridLayout_10"))
        self.horizontalLayout_12 = QtGui.QHBoxLayout()
        self.horizontalLayout_12.setSpacing(5)
        self.horizontalLayout_12.setContentsMargins(-1, 5, -1, 5)
        self.horizontalLayout_12.setObjectName(
            _fromUtf8("horizontalLayout_12"))
        self.gridLayout_10.addLayout(self.horizontalLayout_12, 6, 0, 1, 1)
        self.horizontalLayout_13 = QtGui.QHBoxLayout()
        self.horizontalLayout_13.setSpacing(5)
        self.horizontalLayout_13.setContentsMargins(-1, 5, -1, 5)
        self.horizontalLayout_13.setObjectName(
            _fromUtf8("horizontalLayout_13"))
        self.label_21 = QtGui.QLabel(self.tab_multiecho)
        self.label_21.setObjectName(_fromUtf8("label_21"))
        self.horizontalLayout_13.addWidget(self.label_21)
        self.lineEdit_13 = QtGui.QLineEdit(self.tab_multiecho)
        self.lineEdit_13.setObjectName(_fromUtf8("lineEdit_13"))
        self.horizontalLayout_13.addWidget(self.lineEdit_13)
        self.gridLayout_10.addLayout(self.horizontalLayout_13, 6, 1, 1, 1)
        self.checkBox_8 = QtGui.QCheckBox(self.tab_multiecho)
        self.checkBox_8.setObjectName(_fromUtf8("checkBox_8"))
        self.gridLayout_10.addWidget(self.checkBox_8, 0, 1, 1, 1)
        self.verticalLayout_14.addLayout(self.gridLayout_10)
        self.tab_fdfoptions.addTab(self.tab_multiecho, _fromUtf8(""))
        self.gridLayout_4.addWidget(self.tab_fdfoptions, 3, 3, 1, 1)
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.pushButton_convert = QtGui.QPushButton(self.tab_fdf)
        self.pushButton_convert.setObjectName(_fromUtf8("pushButton_convert"))
        self.verticalLayout_4.addWidget(self.pushButton_convert)
        self.pushButton_check = QtGui.QPushButton(self.tab_fdf)
        self.pushButton_check.setEnabled(False)
        self.pushButton_check.setObjectName(_fromUtf8("pushButton_check"))
        self.verticalLayout_4.addWidget(self.pushButton_check)
        self.pushButton_view = QtGui.QPushButton(self.tab_fdf)
        self.pushButton_view.setEnabled(False)
        self.pushButton_view.setObjectName(_fromUtf8("pushButton_view"))
        self.verticalLayout_4.addWidget(self.pushButton_view)
        self.pushButton_send2daris = QtGui.QPushButton(self.tab_fdf)
        self.pushButton_send2daris.setEnabled(False)
        self.pushButton_send2daris.setObjectName(
            _fromUtf8("pushButton_send2daris"))
        self.verticalLayout_4.addWidget(self.pushButton_send2daris)
        self.gridLayout_4.addLayout(self.verticalLayout_4, 3, 4, 1, 1)
        self.verticalLayout_12.addLayout(self.gridLayout_4)
        self.tabWidget.addTab(self.tab_fdf, _fromUtf8(""))
        self.tab_fid_2 = QtGui.QWidget()
        self.tab_fid_2.setObjectName(_fromUtf8("tab_fid_2"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.tab_fid_2)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.gridLayout_11 = QtGui.QGridLayout()
        self.gridLayout_11.setObjectName(_fromUtf8("gridLayout_11"))
        self.label_22 = QtGui.QLabel(self.tab_fid_2)
        self.label_22.setObjectName(_fromUtf8("label_22"))
        self.gridLayout_11.addWidget(self.label_22, 2, 1, 1, 1)
        self.label_23 = QtGui.QLabel(self.tab_fid_2)
        self.label_23.setObjectName(_fromUtf8("label_23"))
        self.gridLayout_11.addWidget(self.label_23, 3, 1, 1, 1)
        self.pushButton_changefid = QtGui.QPushButton(self.tab_fid_2)
        self.pushButton_changefid.setObjectName(
            _fromUtf8("pushButton_changefid"))
        self.gridLayout_11.addWidget(self.pushButton_changefid, 0, 3, 1, 1)
        self.label_24 = QtGui.QLabel(self.tab_fid_2)
        self.label_24.setObjectName(_fromUtf8("label_24"))
        self.gridLayout_11.addWidget(self.label_24, 0, 1, 1, 1)
        self.label_27 = QtGui.QLabel(self.tab_fid_2)
        self.label_27.setObjectName(_fromUtf8("label_27"))
        self.gridLayout_11.addWidget(self.label_27, 1, 1, 1, 1)
        self.pushButton_changedicom2 = QtGui.QPushButton(self.tab_fid_2)
        self.pushButton_changedicom2.setObjectName(
            _fromUtf8("pushButton_changedicom2"))
        self.gridLayout_11.addWidget(self.pushButton_changedicom2, 1, 3, 1, 1)
        self.lineEdit_darisid2 = QtGui.QLineEdit(self.tab_fid_2)
        self.lineEdit_darisid2.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_darisid2.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_darisid2.setObjectName(_fromUtf8("lineEdit_darisid2"))
        self.gridLayout_11.addWidget(self.lineEdit_darisid2, 2, 2, 1, 1)
        self.lineEdit_fidpath = QtGui.QLineEdit(self.tab_fid_2)
        self.lineEdit_fidpath.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_fidpath.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_fidpath.setObjectName(_fromUtf8("lineEdit_fidpath"))
        self.gridLayout_11.addWidget(self.lineEdit_fidpath, 0, 2, 1, 1)
        self.lineEdit_dicompath2 = QtGui.QLineEdit(self.tab_fid_2)
        self.lineEdit_dicompath2.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_dicompath2.setText(_fromUtf8(""))
        self.lineEdit_dicompath2.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_dicompath2.setObjectName(
            _fromUtf8("lineEdit_dicompath2"))
        self.gridLayout_11.addWidget(self.lineEdit_dicompath2, 1, 2, 1, 1)
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setObjectName(_fromUtf8("verticalLayout_5"))
        self.pushButton_convertfid = QtGui.QPushButton(self.tab_fid_2)
        self.pushButton_convertfid.setObjectName(
            _fromUtf8("pushButton_convertfid"))
        self.verticalLayout_5.addWidget(self.pushButton_convertfid)
        self.pushButton_check2 = QtGui.QPushButton(self.tab_fid_2)
        self.pushButton_check2.setEnabled(False)
        self.pushButton_check2.setObjectName(_fromUtf8("pushButton_check2"))
        self.verticalLayout_5.addWidget(self.pushButton_check2)
        self.pushButton_view2 = QtGui.QPushButton(self.tab_fid_2)
        self.pushButton_view2.setEnabled(False)
        self.pushButton_view2.setObjectName(_fromUtf8("pushButton_view2"))
        self.verticalLayout_5.addWidget(self.pushButton_view2)
        self.pushButton_send2daris2 = QtGui.QPushButton(self.tab_fid_2)
        self.pushButton_send2daris2.setEnabled(False)
        self.pushButton_send2daris2.setObjectName(
            _fromUtf8("pushButton_send2daris2"))
        self.verticalLayout_5.addWidget(self.pushButton_send2daris2)
        self.gridLayout_11.addLayout(self.verticalLayout_5, 3, 3, 1, 1)
        self.verticalLayout_8 = QtGui.QVBoxLayout()
        self.verticalLayout_8.setObjectName(_fromUtf8("verticalLayout_8"))
        self.tabWidget_3 = QtGui.QTabWidget(self.tab_fid_2)
        self.tabWidget_3.setEnabled(True)
        self.tabWidget_3.setObjectName(_fromUtf8("tabWidget_3"))
        self.tab_3 = QtGui.QWidget()
        self.tab_3.setObjectName(_fromUtf8("tab_3"))
        self.verticalLayout_6 = QtGui.QVBoxLayout(self.tab_3)
        self.verticalLayout_6.setObjectName(_fromUtf8("verticalLayout_6"))
        self.gridLayout_12 = QtGui.QGridLayout()
        self.gridLayout_12.setObjectName(_fromUtf8("gridLayout_12"))
        self.horizontalLayout_28 = QtGui.QHBoxLayout()
        self.horizontalLayout_28.setSpacing(5)
        self.horizontalLayout_28.setContentsMargins(-1, 5, -1, 5)
        self.horizontalLayout_28.setObjectName(
            _fromUtf8("horizontalLayout_28"))
        self.gridLayout_12.addLayout(self.horizontalLayout_28, 1, 1, 1, 1)
        self.horizontalLayout_32 = QtGui.QHBoxLayout()
        self.horizontalLayout_32.setSpacing(5)
        self.horizontalLayout_32.setContentsMargins(-1, 5, -1, 5)
        self.horizontalLayout_32.setObjectName(
            _fromUtf8("horizontalLayout_32"))
        self.label_30 = QtGui.QLabel(self.tab_3)
        self.label_30.setObjectName(_fromUtf8("label_30"))
        self.horizontalLayout_32.addWidget(self.label_30)
        spacerItem = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_32.addItem(spacerItem)
        self.lineEdit_gorder = QtGui.QLineEdit(self.tab_3)
        self.lineEdit_gorder.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_gorder.setObjectName(_fromUtf8("lineEdit_gorder"))
        self.horizontalLayout_32.addWidget(self.lineEdit_gorder)
        self.gridLayout_12.addLayout(self.horizontalLayout_32, 2, 0, 1, 1)
        self.horizontalLayout_33 = QtGui.QHBoxLayout()
        self.horizontalLayout_33.setSpacing(5)
        self.horizontalLayout_33.setContentsMargins(-1, 5, -1, 5)
        self.horizontalLayout_33.setObjectName(
            _fromUtf8("horizontalLayout_33"))
        self.label_31 = QtGui.QLabel(self.tab_3)
        self.label_31.setObjectName(_fromUtf8("label_31"))
        self.horizontalLayout_33.addWidget(self.label_31)
        spacerItem1 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_33.addItem(spacerItem1)
        self.lineEdit_gsigma = QtGui.QLineEdit(self.tab_3)
        self.lineEdit_gsigma.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_gsigma.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_gsigma.setObjectName(_fromUtf8("lineEdit_gsigma"))
        self.horizontalLayout_33.addWidget(self.lineEdit_gsigma)
        self.gridLayout_12.addLayout(self.horizontalLayout_33, 1, 0, 1, 1)
        self.checkBox_gaussian = QtGui.QCheckBox(self.tab_3)
        self.checkBox_gaussian.setObjectName(_fromUtf8("checkBox_gaussian"))
        self.gridLayout_12.addWidget(self.checkBox_gaussian, 0, 0, 1, 1)
        self.horizontalLayout_34 = QtGui.QHBoxLayout()
        self.horizontalLayout_34.setSpacing(5)
        self.horizontalLayout_34.setContentsMargins(-1, 5, -1, 5)
        self.horizontalLayout_34.setObjectName(
            _fromUtf8("horizontalLayout_34"))
        self.label_32 = QtGui.QLabel(self.tab_3)
        self.label_32.setObjectName(_fromUtf8("label_32"))
        self.horizontalLayout_34.addWidget(self.label_32)
        self.nearest = QtGui.QRadioButton(self.tab_3)
        self.nearest.setChecked(True)
        self.nearest.setObjectName(_fromUtf8("nearest"))
        self.horizontalLayout_34.addWidget(self.nearest)
        self.reflect = QtGui.QRadioButton(self.tab_3)
        self.reflect.setObjectName(_fromUtf8("reflect"))
        self.horizontalLayout_34.addWidget(self.reflect)
        self.wrap = QtGui.QRadioButton(self.tab_3)
        self.wrap.setObjectName(_fromUtf8("wrap"))
        self.horizontalLayout_34.addWidget(self.wrap)
        spacerItem2 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_34.addItem(spacerItem2)
        self.gridLayout_12.addLayout(self.horizontalLayout_34, 3, 0, 1, 1)
        self.verticalLayout_6.addLayout(self.gridLayout_12)
        self.tabWidget_3.addTab(self.tab_3, _fromUtf8(""))
        self.tab_4 = QtGui.QWidget()
        self.tab_4.setObjectName(_fromUtf8("tab_4"))
        self.verticalLayout_7 = QtGui.QVBoxLayout(self.tab_4)
        self.verticalLayout_7.setObjectName(_fromUtf8("verticalLayout_7"))
        self.verticalLayout_9 = QtGui.QVBoxLayout()
        self.verticalLayout_9.setObjectName(_fromUtf8("verticalLayout_9"))
        self.checkBox_median = QtGui.QCheckBox(self.tab_4)
        self.checkBox_median.setObjectName(_fromUtf8("checkBox_median"))
        self.verticalLayout_9.addWidget(self.checkBox_median)
        self.horizontalLayout_41 = QtGui.QHBoxLayout()
        self.horizontalLayout_41.setObjectName(
            _fromUtf8("horizontalLayout_41"))
        self.label_33 = QtGui.QLabel(self.tab_4)
        self.label_33.setObjectName(_fromUtf8("label_33"))
        self.horizontalLayout_41.addWidget(self.label_33)
        spacerItem3 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_41.addItem(spacerItem3)
        self.lineEdit_median_size = QtGui.QLineEdit(self.tab_4)
        self.lineEdit_median_size.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_median_size.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_median_size.setObjectName(
            _fromUtf8("lineEdit_median_size"))
        self.horizontalLayout_41.addWidget(self.lineEdit_median_size)
        self.verticalLayout_9.addLayout(self.horizontalLayout_41)
        self.verticalLayout_7.addLayout(self.verticalLayout_9)
        self.tabWidget_3.addTab(self.tab_4, _fromUtf8(""))
        self.tab_8 = QtGui.QWidget()
        self.tab_8.setEnabled(True)
        self.tab_8.setObjectName(_fromUtf8("tab_8"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.tab_8)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.verticalLayout_10 = QtGui.QVBoxLayout()
        self.verticalLayout_10.setObjectName(_fromUtf8("verticalLayout_10"))
        self.checkBox_wiener = QtGui.QCheckBox(self.tab_8)
        self.checkBox_wiener.setEnabled(True)
        self.checkBox_wiener.setObjectName(_fromUtf8("checkBox_wiener"))
        self.verticalLayout_10.addWidget(self.checkBox_wiener)
        self.horizontalLayout_42 = QtGui.QHBoxLayout()
        self.horizontalLayout_42.setObjectName(
            _fromUtf8("horizontalLayout_42"))
        self.label_34 = QtGui.QLabel(self.tab_8)
        self.label_34.setEnabled(True)
        self.label_34.setObjectName(_fromUtf8("label_34"))
        self.horizontalLayout_42.addWidget(self.label_34)
        spacerItem4 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_42.addItem(spacerItem4)
        self.lineEdit_wiener_size = QtGui.QLineEdit(self.tab_8)
        self.lineEdit_wiener_size.setEnabled(True)
        self.lineEdit_wiener_size.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_wiener_size.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_wiener_size.setObjectName(
            _fromUtf8("lineEdit_wiener_size"))
        self.horizontalLayout_42.addWidget(self.lineEdit_wiener_size)
        self.verticalLayout_10.addLayout(self.horizontalLayout_42)
        self.horizontalLayout_43 = QtGui.QHBoxLayout()
        self.horizontalLayout_43.setObjectName(
            _fromUtf8("horizontalLayout_43"))
        self.label_35 = QtGui.QLabel(self.tab_8)
        self.label_35.setEnabled(True)
        self.label_35.setObjectName(_fromUtf8("label_35"))
        self.horizontalLayout_43.addWidget(self.label_35)
        spacerItem5 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_43.addItem(spacerItem5)
        self.lineEdit_wiener_noise = QtGui.QLineEdit(self.tab_8)
        self.lineEdit_wiener_noise.setEnabled(True)
        self.lineEdit_wiener_noise.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.lineEdit_wiener_noise.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.lineEdit_wiener_noise.setObjectName(
            _fromUtf8("lineEdit_wiener_noise"))
        self.horizontalLayout_43.addWidget(self.lineEdit_wiener_noise)
        self.verticalLayout_10.addLayout(self.horizontalLayout_43)
        self.verticalLayout_3.addLayout(self.verticalLayout_10)
        self.tabWidget_3.addTab(self.tab_8, _fromUtf8(""))
        self.verticalLayout_8.addWidget(self.tabWidget_3)
        self.label_36 = QtGui.QLabel(self.tab_fid_2)
        self.label_36.setObjectName(_fromUtf8("label_36"))
        self.verticalLayout_8.addWidget(self.label_36)
        self.horizontalLayout_35 = QtGui.QHBoxLayout()
        self.horizontalLayout_35.setObjectName(
            _fromUtf8("horizontalLayout_35"))
        self.checkBox_magn = QtGui.QCheckBox(self.tab_fid_2)
        self.checkBox_magn.setEnabled(True)
        self.checkBox_magn.setCheckable(True)
        self.checkBox_magn.setChecked(True)
        self.checkBox_magn.setObjectName(_fromUtf8("checkBox_magn"))
        self.horizontalLayout_35.addWidget(self.checkBox_magn)
        self.checkBox_pha = QtGui.QCheckBox(self.tab_fid_2)
        self.checkBox_pha.setObjectName(_fromUtf8("checkBox_pha"))
        self.horizontalLayout_35.addWidget(self.checkBox_pha)
        self.verticalLayout_8.addLayout(self.horizontalLayout_35)
        self.horizontalLayout_36 = QtGui.QHBoxLayout()
        self.horizontalLayout_36.setObjectName(
            _fromUtf8("horizontalLayout_36"))
        self.checkBox_reimag = QtGui.QCheckBox(self.tab_fid_2)
        self.checkBox_reimag.setObjectName(_fromUtf8("checkBox_reimag"))
        self.horizontalLayout_36.addWidget(self.checkBox_reimag)
        self.checkBox_ksp = QtGui.QCheckBox(self.tab_fid_2)
        self.checkBox_ksp.setObjectName(_fromUtf8("checkBox_ksp"))
        self.horizontalLayout_36.addWidget(self.checkBox_ksp)
        self.verticalLayout_8.addLayout(self.horizontalLayout_36)
        self.gridLayout_11.addLayout(self.verticalLayout_8, 3, 2, 1, 1)
        self.verticalLayout_2.addLayout(self.gridLayout_11)
        self.tabWidget.addTab(self.tab_fid_2, _fromUtf8(""))
        self.verticalLayout.addWidget(self.tabWidget)
        self.buttonBox = QtGui.QDialogButtonBox(Form)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(
            QtGui.QDialogButtonBox.Close | QtGui.QDialogButtonBox.Help)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Form)
        self.tabWidget.setCurrentIndex(0)
        self.tab_fdfoptions.setCurrentIndex(0)
        self.tabWidget_3.setCurrentIndex(2)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate(
            "Form", "MBI\'s Agilent to Dicom converter application (1.2.1)", None))
        self.tabWidget.setToolTip(
            _translate("Form", "MBI\'s Agilent to Dicom application", None))
        self.label_12.setText(_translate("Form", "Options", None))
        self.pushButton_changefdf.setToolTip(
            _translate("Form", "Change the FDF input directory", None))
        self.pushButton_changefdf.setText(
            _translate("Form", "Change Dir", None))
        self.label_13.setText(_translate("Form", "FDF folder", None))
        self.label_16.setText(_translate("Form", "Dicom folder", None))
        self.label_17.setText(_translate("Form", "DaRIS ID", None))
        self.pushButton_changedicom.setToolTip(_translate(
            "Form", "Change the output DICOM directory.  Only do this if the automatic folder name is already taken.", None))
        self.pushButton_changedicom.setText(
            _translate("Form", "Change Dir", None))
        self.lineEdit_darisid.setToolTip(_translate(
            "Form", "DaRIS ID string should automatically update with new FDF folder", None))
        self.checkBox_debugging.setToolTip(
            _translate("Form", "Display more verbose infomration for debugging", None))
        self.checkBox_debugging.setText(
            _translate("Form", "Show debugging", None))
        self.checkBox_nodcmulti.setToolTip(_translate(
            "Form", "Do not convert 2D slices to enhance DICOM format.  Some sequences may prefer this including ASL.", None))
        self.checkBox_nodcmulti.setText(
            _translate("Form", "Keep 2D DICOM slices", None))
        self.tab_fdfoptions.setTabText(self.tab_fdfoptions.indexOf(
            self.tab_generic), _translate("Form", "Generic", None))
        self.radioButton_bmatrixdiff.setText(
            _translate("Form", "Bmatrix", None))
        self.radioButton_directionaldiff.setText(
            _translate("Form", "Directional", None))
        self.label_18.setText(_translate("Form", "TextLabel", None))
        self.checkBox_diffusion1.setText(_translate("Form", "Diffusion", None))
        self.tab_fdfoptions.setTabText(self.tab_fdfoptions.indexOf(
            self.tab_diffusion), _translate("Form", "Diffusion", None))
        self.label_21.setText(_translate("Form", "Echo Train List", None))
        self.checkBox_8.setText(_translate("Form", "Multi-echo", None))
        self.tab_fdfoptions.setTabText(self.tab_fdfoptions.indexOf(
            self.tab_multiecho), _translate("Form", "Multi echo", None))
        self.pushButton_convert.setToolTip(
            _translate("Form", "Convert FDF folder to DICOM", None))
        self.pushButton_convert.setText(_translate("Form", "Convert", None))
        self.pushButton_check.setText(_translate("Form", "Check", None))
        self.pushButton_view.setText(_translate("Form", "View", None))
        self.pushButton_send2daris.setToolTip(
            _translate("Form", "Send Dicom folder to DaRIS ", None))
        self.pushButton_send2daris.setText(
            _translate("Form", "Send to DaRIS", None))
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab_fdf), _translate("Form", "FDF converter", None))
        self.label_22.setText(_translate("Form", "DaRIS ID", None))
        self.label_23.setToolTip(
            _translate("Form", "Smoothing filter options ", None))
        self.label_23.setText(_translate("Form", "Filter \n"
                                         " Options", None))
        self.pushButton_changefid.setToolTip(
            _translate("Form", "Change the input FID directory", None))
        self.pushButton_changefid.setText(
            _translate("Form", "Change Dir", None))
        self.label_24.setText(_translate("Form", "FID folder", None))
        self.label_27.setText(_translate("Form", "Dicom folder", None))
        self.pushButton_changedicom2.setToolTip(
            _translate("Form", "Change the output DICOM directory.", None))
        self.pushButton_changedicom2.setText(
            _translate("Form", "Change Dir", None))
        self.lineEdit_darisid2.setToolTip(_translate(
            "Form", "DaRIS ID generated from header information in FID procpar.  Do not edit unless you know what you are doing!", None))
        self.lineEdit_fidpath.setToolTip(
            _translate("Form", "FID directory. To change, use the change dir button.", None))
        self.lineEdit_dicompath2.setToolTip(_translate(
            "Form", "Output dicom directory. Allow input FID directory to determine this path.  Otherwise use the change dir button.", None))
        self.pushButton_convertfid.setToolTip(
            _translate("Form", "Run the fid2dcm script", None))
        self.pushButton_convertfid.setText(_translate("Form", "Convert", None))
        self.pushButton_check2.setToolTip(
            _translate("Form", "Check the unfiltered reconstruction DICM", None))
        self.pushButton_check2.setText(_translate("Form", "Check", None))
        self.pushButton_view2.setToolTip(
            _translate("Form", "View the filtered and unfiltered images", None))
        self.pushButton_view2.setText(_translate("Form", "View", None))
        self.pushButton_send2daris2.setToolTip(
            _translate("Form", "Send DICOMs to Daris", None))
        self.pushButton_send2daris2.setText(
            _translate("Form", "Send to DaRIS", None))
        self.tabWidget_3.setToolTip(_translate(
            "Form", "Multidimensional median, Gaussian or Wiener  filters.", None))
        self.label_30.setToolTip(_translate("Form", "order : {0, 1, 2, 3} or sequence from same set, optional\n"
                                            "\n"
                                            "    The order of the filter along each axis is given as a sequence of integers, or as a single number. An order of 0 corresponds to convolution with a Gaussian kernel. An order of 1, 2, or 3 corresponds to convolution with the first, second or third derivatives of a Gaussian. Higher order derivatives are not implemented\n"
                                            "", None))
        self.label_30.setText(_translate("Form", "Order", None))
        self.lineEdit_gorder.setToolTip(_translate("Form", "order : {0, 1, 2, 3} or sequence from same set, optional\n"
                                                   "\n"
                                                   "    The order of the filter along each axis is given as a sequence of integers, or as a single number. An order of 0 corresponds to convolution with a Gaussian kernel. An order of 1, 2, or 3 corresponds to convolution with the first, second or third derivatives of a Gaussian. Higher order derivatives are not implemented\n"
                                                   "", None))
        self.lineEdit_gorder.setText(_translate("Form", "0", None))
        self.label_31.setToolTip(_translate("Form", "scalar or sequence of scalars\n"
                                            "\n"
                                            "    Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.", None))
        self.label_31.setText(_translate("Form", "Sigma", None))
        self.lineEdit_gsigma.setToolTip(_translate("Form", "scalar or sequence of scalars\n"
                                                   "\n"
                                                   "    Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.", None))
        self.lineEdit_gsigma.setText(_translate("Form", "0.707", None))
        self.checkBox_gaussian.setToolTip(_translate("Form", " scipy.ndimage.filters.gaussian_filter(input, sigma, order=0, output=None, mode=\'reflect\', cval=0.0, truncate=4.0)[source]\n"
                                                     "\n"
                                                     "    Multidimensional Gaussian filter.", None))
        self.checkBox_gaussian.setText(
            _translate("Form", "Use Gaussian Filter", None))
        self.label_32.setToolTip(_translate("Form", "scalar or sequence of scalars\n"
                                            "\n"
                                            "    Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.", None))
        self.label_32.setText(_translate("Form", "Mode:", None))
        self.nearest.setText(_translate("Form", "Nearest", None))
        self.reflect.setText(_translate("Form", "Reflect", None))
        self.wrap.setText(_translate("Form", "Wrap", None))
        self.tabWidget_3.setTabText(
            self.tabWidget_3.indexOf(self.tab_3), _translate("Form", "Gaussian", None))
        self.checkBox_median.setToolTip(_translate("Form", "NOT IN USE:\n"
                                                   "scipy.ndimage.filters.median_filter¶\n"
                                                   "\n"
                                                   "scipy.ndimage.filters.median_filter(input, size=None, footprint=None, output=None, mode=\'reflect\', cval=0.0, origin=0)[source]\n"
                                                   "\n"
                                                   "    Calculates a multidimensional median filter.", None))
        self.checkBox_median.setText(
            _translate("Form", "Use Median Filter", None))
        self.label_33.setToolTip(_translate("Form", "     \n"
                                            "\n"
                                            "size : scalar or tuple, optional\n"
                                            "\n"
                                            "    See footprint, below\n"
                                            "\n"
                                            "footprint : array, optional\n"
                                            "\n"
                                            "    Either size or footprint must be defined. size gives the shape that is taken from the input array, at every element position,\n"
                                            " to define the input to the filter function. footprint is a boolean array that specifies (implicitly) a shape,\n"
                                            " but also which of the elements within this shape will get passed to the filter function. Thus size=(n, m) is equivalent to footprint = np.ones((n, m)).\n"
                                            " We adjust size to the number of dimensions of the input array, so that, if the input array is shape (10, 10, 10),\n"
                                            " and size is 2, then the actual size used is (2, 2, 2).\n"
                                            "", None))
        self.label_33.setText(_translate("Form", "Window Size(s)", None))
        self.lineEdit_median_size.setToolTip(_translate("Form", "     \n"
                                                        "\n"
                                                        "size : scalar or tuple, optional\n"
                                                        "\n"
                                                        "    See footprint, below\n"
                                                        "\n"
                                                        "footprint : array, optional\n"
                                                        "\n"
                                                        "    Either size or footprint must be defined. size gives the shape that is taken from the input array, at every element position,\n"
                                                        " to define the input to the filter function. footprint is a boolean array that specifies (implicitly) a shape,\n"
                                                        " but also which of the elements within this shape will get passed to the filter function. Thus size=(n, m) is equivalent to footprint = np.ones((n, m)).\n"
                                                        " We adjust size to the number of dimensions of the input array, so that, if the input array is shape (10, 10, 10),\n"
                                                        " and size is 2, then the actual size used is (2, 2, 2).\n"
                                                        "", None))
        self.lineEdit_median_size.setText(_translate("Form", "5", None))
        self.tabWidget_3.setTabText(
            self.tabWidget_3.indexOf(self.tab_4), _translate("Form", "Median", None))
        self.checkBox_wiener.setToolTip(_translate("Form", "NOT ENABLED - Speak to MBI Imaging Team\n"
                                                   "\n"
                                                   "scipy.signal.wiener(im, mysize = None, noise = None)\n"
                                                   "\n"
                                                   "    Perform a Wiener filter on an N-dimensional array.\n"
                                                   "\n"
                                                   "    Apply a Wiener filter to the N-dimensional array im.", None))
        self.checkBox_wiener.setText(
            _translate("Form", "Use Wiener Filter", None))
        self.label_34.setToolTip(_translate("Form", "mysize : int or arraylike, optional\n"
                                            "\n"
                                            "    A scalar or an N-length list giving the size of the Wiener filter window in each dimension. Elements of mysize should be odd. If mysize is a scalar, then this scalar is used as the size in each dimension.\n"
                                            "", None))
        self.label_34.setText(_translate("Form", " Window Size(s)", None))
        self.lineEdit_wiener_size.setToolTip(_translate("Form", "mysize : int or arraylike, optional\n"
                                                        "\n"
                                                        "    A scalar or an N-length list giving the size of the Wiener filter window in each dimension. Elements of mysize should be odd. If mysize is a scalar, then this scalar is used as the size in each dimension.\n"
                                                        "", None))
        self.lineEdit_wiener_size.setText(_translate("Form", "5", None))
        self.label_35.setToolTip(_translate(
            "Form", " < html> < head/> < body> < p > noise: Estimation of noise, Set to 0 for local variance to be used. < /p> < /body> < /html > ", None))
        self.label_35.setText(
            _translate("Form", "Noise (Est. variance)", None))
        self.lineEdit_wiener_noise.setToolTip(_translate("Form", "mysize : int or arraylike, optional\n"
                                                         "\n"
                                                         "    A scalar or an N-length list giving the size of the Wiener filter window in each dimension. Elements of mysize should be odd. If mysize is a scalar, then this scalar is used as the size in each dimension.\n"
                                                         "", None))
        self.lineEdit_wiener_noise.setText(_translate("Form", "0", None))
        self.tabWidget_3.setTabText(
            self.tabWidget_3.indexOf(self.tab_8), _translate("Form", "Wiener", None))
        self.label_36.setText(_translate("Form", "Outputs of Filter", None))
        self.checkBox_magn.setToolTip(
            _translate("Form", "Save the magnitude of the complex filtered image", None))
        self.checkBox_magn.setText(_translate("Form", "Save Magnitude", None))
        self.checkBox_pha.setToolTip(_translate(
            "Form", "Save the phase of the complex filtered image to \" < outputpath > - < filtertype > -pha.dcm\"", None))
        self.checkBox_pha.setText(_translate("Form", "Save Phase", None))
        self.checkBox_reimag.setToolTip(_translate(
            "Form", "Save the real and imaginary components of the complex filtered image to \" < outputpath > - < filtertype > - < real or imag > .dcm\"", None))
        self.checkBox_reimag.setText(
            _translate("Form", "Save REAL and IMAG", None))
        self.checkBox_ksp.setToolTip(
            _translate("Form", "Save k space data to MATLAB mat file.", None))
        self.checkBox_ksp.setText(
            _translate("Form", "Save K space data", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(
            self.tab_fid_2), _translate("Form", "FID converter", None))

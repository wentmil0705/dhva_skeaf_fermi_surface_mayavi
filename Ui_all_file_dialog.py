# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/Users/wentworth/Desktop/skeaf_demo_sw/skeaf_demo_2/all_file_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_all_file_dialog(object):
    def setupUi(self, all_file_dialog):
        all_file_dialog.setObjectName("all_file_dialog")
        all_file_dialog.resize(396, 395)
        self.layoutWidget = QtWidgets.QWidget(all_file_dialog)
        self.layoutWidget.setGeometry(QtCore.QRect(21, 22, 361, 361))
        self.layoutWidget.setObjectName("layoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.layoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.bxsf_layout = QtWidgets.QVBoxLayout()
        self.bxsf_layout.setObjectName("bxsf_layout")
        self.bxsf_button = QtWidgets.QPushButton(self.layoutWidget)
        self.bxsf_button.setObjectName("bxsf_button")
        self.bxsf_layout.addWidget(self.bxsf_button)
        self.bxsf_label = QtWidgets.QLabel(self.layoutWidget)
        self.bxsf_label.setMinimumSize(QtCore.QSize(100, 0))
        self.bxsf_label.setText("")
        self.bxsf_label.setObjectName("bxsf_label")
        self.bxsf_layout.addWidget(self.bxsf_label)
        self.gridLayout.addLayout(self.bxsf_layout, 0, 0, 1, 1)
        self.long_layout = QtWidgets.QVBoxLayout()
        self.long_layout.setObjectName("long_layout")
        self.long_button = QtWidgets.QPushButton(self.layoutWidget)
        self.long_button.setObjectName("long_button")
        self.long_layout.addWidget(self.long_button)
        self.long_label = QtWidgets.QLabel(self.layoutWidget)
        self.long_label.setMinimumSize(QtCore.QSize(100, 0))
        self.long_label.setText("")
        self.long_label.setObjectName("long_label")
        self.long_layout.addWidget(self.long_label)
        self.gridLayout.addLayout(self.long_layout, 1, 0, 1, 1)
        self.ang_layout = QtWidgets.QVBoxLayout()
        self.ang_layout.setSpacing(20)
        self.ang_layout.setObjectName("ang_layout")
        self.ang_button = QtWidgets.QPushButton(self.layoutWidget)
        self.ang_button.setObjectName("ang_button")
        self.ang_layout.addWidget(self.ang_button)
        self.ang_label = QtWidgets.QLabel(self.layoutWidget)
        self.ang_label.setMinimumSize(QtCore.QSize(100, 0))
        self.ang_label.setText("")
        self.ang_label.setObjectName("ang_label")
        self.ang_layout.addWidget(self.ang_label)
        self.gridLayout.addLayout(self.ang_layout, 2, 0, 1, 1)
        self.au_layout = QtWidgets.QVBoxLayout()
        self.au_layout.setSpacing(20)
        self.au_layout.setObjectName("au_layout")
        self.au_button = QtWidgets.QPushButton(self.layoutWidget)
        self.au_button.setObjectName("au_button")
        self.au_layout.addWidget(self.au_button)
        self.au_label = QtWidgets.QLabel(self.layoutWidget)
        self.au_label.setMinimumSize(QtCore.QSize(100, 0))
        self.au_label.setText("")
        self.au_label.setObjectName("au_label")
        self.au_layout.addWidget(self.au_label)
        self.gridLayout.addLayout(self.au_layout, 3, 0, 1, 1)
        self.import_button = QtWidgets.QPushButton(self.layoutWidget)
        self.import_button.setObjectName("import_button")
        self.gridLayout.addWidget(self.import_button, 4, 0, 1, 1)

        self.retranslateUi(all_file_dialog)
        QtCore.QMetaObject.connectSlotsByName(all_file_dialog)

    def retranslateUi(self, all_file_dialog):
        _translate = QtCore.QCoreApplication.translate
        all_file_dialog.setWindowTitle(_translate("all_file_dialog", "Dialog"))
        self.bxsf_button.setText(_translate("all_file_dialog", "BXSF file"))
        self.long_button.setText(_translate("all_file_dialog", "results_long.out"))
        self.ang_button.setText(_translate("all_file_dialog", "results_orbitoutlines_invAng.out"))
        self.au_button.setText(_translate("all_file_dialog", "results_orbitoutlines_invau.out"))
        self.import_button.setText(_translate("all_file_dialog", "Import"))

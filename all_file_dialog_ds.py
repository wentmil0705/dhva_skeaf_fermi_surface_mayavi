import imp
from PyQt5 import QtCore, QtGui, QtWidgets
from Ui_all_file_dialog import Ui_all_file_dialog
    
class all_file_dialog(QtWidgets.QDialog, Ui_all_file_dialog):
    _signal = QtCore.pyqtSignal(str,str,str,str)
    def __init__(self):
        super(all_file_dialog, self).__init__()
        self.setupUi(self)
        self.retranslateUi(self)
        self.import_button.clicked.connect(self.returnFilePath)
        self.bxsf_button.clicked.connect(lambda: self.showNormalFile(self.bxsf_label))
        self.long_button.clicked.connect(lambda: self.showNormalFile(self.long_label))
        self.ang_button.clicked.connect(lambda: self.showNormalFile(self.ang_label))
        self.au_button.clicked.connect(lambda: self.showNormalFile(self.au_label))

    def showNormalFile(self, lbl):
        fname = QtWidgets.QFileDialog.getOpenFileName(None, 'Open file', '/home')


        if fname[0]:
            lbl.setText(fname[0])
            return fname[0]
        else:
            return ''


    def returnFilePath(self):
        self._signal.emit(self.bxsf_label.text(), self.long_label.text(), self.ang_label.text(), self.au_label.text())
        self.close()
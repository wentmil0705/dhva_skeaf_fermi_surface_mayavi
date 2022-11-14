__author__ = "Guanzhang Liu"
__maintainer__ = "Guanzhang Liu"
__email__ = "MG21340029@smail.nju.edu.cn"
__date__ = "October 23, 2022"

'''
all_file_dialog_ds.py为Ui_all_file_dialog.py界面对应功能的设计文件
'''

from PyQt5 import QtCore, QtGui, QtWidgets
from ui_py_files.Ui_all_file_dialog import Ui_all_file_dialog   # 导入Ui_all_file_dialog界面设计类

'''all_file_dialog类继承自QDialog属性，并承接Ui_all_file_dialog的界面'''    
class all_file_dialog(QtWidgets.QDialog, Ui_all_file_dialog):

    '''利用_signal对主窗口与子窗口进行信息传递'''
    _signal = QtCore.pyqtSignal(str,str,str,str)


    '''初始化，继承父类Ui_all_file_dialog并进行细节设计'''
    def __init__(self):
        super(all_file_dialog, self).__init__()
        self.setupUi(self)
        self.retranslateUi(self)
        ## import_button用于向主窗口发射_signal并关闭子窗口
        self.import_button.clicked.connect(self.returnFilePath)
        ## 每一个文件输入按钮都会调用文件选取函数
        self.bxsf_button.clicked.connect(lambda: self.selectFile(self.bxsf_label))
        self.long_button.clicked.connect(lambda: self.selectFile(self.long_label))
        self.ang_button.clicked.connect(lambda: self.selectFile(self.ang_label))
        self.au_button.clicked.connect(lambda: self.selectFile(self.au_label))

    '''
    文件选择函数：
    lbl：选择文件后所用来显示绝对路径的区域
    '''
    def selectFile(self, lbl):
        fname = QtWidgets.QFileDialog.getOpenFileName(None, 'Open file', '/Users/wentworth/Desktop')

        if fname[0]:
            lbl.setText(fname[0])
            return fname[0]
        else:
            return ''

    '''
    向主窗口发送所选文件路径信息函数
    '''
    def returnFilePath(self):
        self._signal.emit(self.bxsf_label.text(), self.long_label.text(), self.ang_label.text(), self.au_label.text())  # 发送信息
        self.close()    # 关闭子窗口
__author__ = "Guanzhang Liu"
__maintainer__ = "Guanzhang Liu"
__email__ = "MG21340029@smail.nju.edu.cn"
__date__ = "October 23, 2022"

'''
B_direction_ds.py为ui_B_direction.py界面对应功能的设计文件
'''

from ui_py_files.ui_B_direction import *    # 导入ui_B_direction界面设计类
from PyQt5 import QtWidgets, QtCore
import numpy as np
from share import *

'''B_direction类继承自QWidget属性，并承接ui_B_direction的界面'''    
class B_direction(QtWidgets.QWidget, Ui_B_direction):
    _angle_signal = QtCore.pyqtSignal(float, float) # 磁场方向角度信号
    _close_signal = QtCore.pyqtSignal(bool)         # 取消显示磁场方向信号

    def __init__(self, bcell, parent=None):
        super(B_direction, self).__init__(parent)
        self.setupUi(self)

        self.bcell = bcell
        self.xb1 = saveVal(0.)
        self.xb2 = saveVal(0.)
        self.xb3 = saveVal(1.)
        self.theta = saveVal(0.)
        self.phi = saveVal(0.)
        self.Theta_line.editingFinished.connect(lambda: self.showFinished(self.Theta_line.text(), self.theta))
        self.Phi_line.editingFinished.connect(lambda: self.showFinished(self.Phi_line.text(), self.phi))
        self.b1_line.editingFinished.connect(lambda: self.showLineFinished(self.b1_line.text(), self.xb1))
        self.b2_line.editingFinished.connect(lambda: self.showLineFinished(self.b2_line.text(), self.xb2))
        self.b3_line.editingFinished.connect(lambda: self.showLineFinished(self.b3_line.text(), self.xb3))
        self.updateB_button.clicked.connect(self.emitSignal)

        self.Theta_line.setText(str(0.))
        self.Phi_line.setText(str(0.))
        self.b1_line.setText(str(0.))
        self.b2_line.setText(str(0.))
        self.b3_line.setText(str(1.))
    
    # 发射信号更新所设定磁场方向
    def emitSignal(self):
        self._angle_signal.emit(self.theta.getVal(), self.phi.getVal())
    
    # 输入完成提示，并保存值
    def showFinished(self, text, angle):
        if isFloat(text):
            angle.savVal(float(text))
        else:
            QtWidgets.QMessageBox.critical(self, "Error", "Input data is not correct!")

    # 输入完成提示，并保存值
    def showLineFinished(self, text, xb):
        if isFloat(text):
            xb.savVal(float(text))
            h = np.array([self.xb1.getVal(), self.xb2.getVal(), self.xb3.getVal()])
            h = np.dot(h, self.bcell)
            theta, phi = calc_angle(self.bcell, h)
            self.theta.savVal(theta)
            self.phi.savVal(phi)
            self.Theta_line.setText('{:.5f}'.format(self.theta.getVal()))
            self.Phi_line.setText('{:.5f}'.format(self.phi.getVal()))
        else:
            QtWidgets.QMessageBox.critical(self, "Error", "Input data is not correct!")
    
    # 重置关闭函数，不再显示磁场方向
    def closeEvent(self, event):
        self._close_signal.emit(True)
        event.accept()
 

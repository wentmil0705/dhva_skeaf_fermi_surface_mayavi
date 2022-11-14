__author__ = "Guanzhang Liu"
__maintainer__ = "Guanzhang Liu"
__email__ = "MG21340029@smail.nju.edu.cn"
__date__ = "October 23, 2022"

'''
more_ds.py为Ui_more.py界面对应功能的设计文件
'''

from PyQt5 import QtCore, QtGui, QtWidgets
from ui_py_files.Ui_more import Ui_More # 导入Ui_more界面设计类
from share import *

'''more类继承自QDialog属性，并承接Ui_more的界面'''    
class more(QtWidgets.QDialog, Ui_More):

    '''利用_signal对主窗口与子窗口进行信息传递'''
    _magnetic_field_signal = QtCore.pyqtSignal(bool)
    _BG_COLOR_signal = QtCore.pyqtSignal(tuple)
    _SLICE_COLOR_signal = QtCore.pyqtSignal(tuple)
    _fs_opacity_signal = QtCore.pyqtSignal(float)
    _line_opacity_signal = QtCore.pyqtSignal(float)
    _slice_opacity_signal = QtCore.pyqtSignal(float)
    _trait_opacity_signal = QtCore.pyqtSignal(float)


    '''初始化，继承父类Ui_more并进行细节设计'''
    def __init__(self):
        super(more, self).__init__()
        self.setupUi(self)
        self.retranslateUi(self)

        ## mayavi界面显示的default背景颜色(1,1,1)
        self.BG_COLOR = saveVal((1, 1, 1))
        ## 切片显示的default颜色(0,0,0)
        self.SLICE_COLOR = saveVal((0, 0, 0))
        ## 磁场显示default不开启
        self.magnetic_field_checkbox.setChecked(False)
        self.magnetic_field_checkbox.stateChanged.connect(self.showMag)
        ## 费米面显示透明度
        self.fs_opacity_slider.valueChanged.connect(self.showFsOpacity)
        ## 边界线显示透明度
        self.line_opacity_slider.valueChanged.connect(self.showLineOpacity)
        ## 轨道显示透明度
        self.trait_opacity_slider.valueChanged.connect(self.showTraitOpacity)
        ## 切片显示透明度
        self.slice_opacity_slider.valueChanged.connect(self.showSliceOpacity)

        ## 对应滑块更改时调用函数
        self.fs_opacity_line.editingFinished.connect(self.showFsOpacity2)
        self.line_opacity_line.editingFinished.connect(self.showLineOpacity2)
        self.trait_opacity_line.editingFinished.connect(self.showTraitOpacity2)
        self.slice_opacity_line.editingFinished.connect(self.showSliceOpacity2)
        
        ## 对应滑块的数字显示区域
        self.fs_opacity_line.setText(str(1.0))
        self.line_opacity_line.setText(str(1.0))
        self.trait_opacity_line.setText(str(1.0))
        self.slice_opacity_line.setText(str(0.4))
        self.fs_opacity_slider.setValue(10)
        self.line_opacity_slider.setValue(10)
        self.trait_opacity_slider.setValue(10)
        self.slice_opacity_slider.setValue(4)

        ## 对应颜色更改调用函数
        self.background_button.clicked.connect(self.showBgColor)
        self.slice_color_button.clicked.connect(self.showSliceColor)


    def showMag(self):
        self._magnetic_field_signal.emit(self.magnetic_field_checkbox.isChecked())
    '''
    滑块滑动时更改滑块数字显示区域函数
    text：滑块值
    line：数字显示区域
    '''
    def showFsOpacity(self):
        self.showValue(self.fs_opacity_slider.value(), self.fs_opacity_line)
        self._fs_opacity_signal.emit(self.fs_opacity_slider.value()/10)
    
    def showLineOpacity(self):
        self.showValue(self.line_opacity_slider.value(), self.line_opacity_line)
        self._line_opacity_signal.emit(self.line_opacity_slider.value()/10)
    
    def showTraitOpacity(self):
        self.showValue(self.trait_opacity_slider.value(), self.trait_opacity_line)
        self._trait_opacity_signal.emit(self.trait_opacity_slider.value()/10)
    
    def showSliceOpacity(self):
        self.showValue(self.slice_opacity_slider.value(), self.slice_opacity_line)
        self._slice_opacity_signal.emit(self.slice_opacity_slider.value()/10)
    
    def showValue(self, text, line):
        line.setText(str(text / 10))
    
    '''
    数字显示区域更改滑块位置函数
    text：数字显示区域值
    slider：所需改变的滑块
    line：数字显示区域
    '''
    def showFsOpacity2(self):
        self.showText(self.fs_opacity_line.text(), self.fs_opacity_slider, self.fs_opacity_line)
        self._fs_opacity_signal.emit(float(self.fs_opacity_line.text()))
    
    def showLineOpacity2(self):
        self.showText(self.line_opacity_line.text(), self.line_opacity_slider, self.line_opacity_line)
        self._line_opacity_signal.emit(float(self.line_opacity_line.text()))
    
    def showTraitOpacity2(self):
        self.showText(self.trait_opacity_line.text(), self.trait_opacity_slider, self.trait_opacity_line)
        self._trait_opacity_signal.emit(float(self.trait_opacity_line.text()))
    
    def showSliceOpacity2(self):
        self.showText(self.slice_opacity_line.text(), self.slice_opacity_slider, self.slice_opacity_line)
        self._slice_opacity_signal.emit(float(self.slice_opacity_line.text()))
    
    def showText(self, text, slider, line):
        ## 判断输入值是否合理
        if notFloat(text):
            QtWidgets.QMessageBox.warning(self, "Error", "Input is WRONG!\nPlease Input number between 0 and 1!")
            line.setText(str(slider.value()/10))
        elif (float(text) > 1) | (float(text) < 0):
            QtWidgets.QMessageBox.warning(self, "Error", "Input is WRONG!\nPlease Input number between 0 and 1!")
            line.setText(str(slider.value()/10))

        else:
            slider.setValue(10*float(text))
            line.setText(str(slider.value()/10))
    
    '''
    颜色选择并显示函数
    frm：颜色显示区域
    colorVal：颜色值记录
    '''
    def showBgColor(self):
        self.selectColor(self.background_color_frame, self.BG_COLOR)
        self._BG_COLOR_signal.emit(self.BG_COLOR.getVal()[0:3])
    
    def showSliceColor(self):
        self.selectColor(self.slice_color_frame, self.SLICE_COLOR)
        self._SLICE_COLOR_signal.emit(self.SLICE_COLOR.getVal()[0:3])

    def selectColor(self, frm, colorVal):
        col = QtWidgets.QColorDialog.getColor()

        if col.isValid():
            frm.setStyleSheet("QWidget { background-color: %s}" % col.name())
            colorVal.savVal(col.getRgbF())
    
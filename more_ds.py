from PyQt5 import QtCore, QtGui, QtWidgets
from Ui_more import Ui_More

class saveVal(object):
    
    def __init__(self, v=0.):
        self.val = v
    
    def savVal(self, v):
        self.val = v
    
    def getVal(self):
        return self.val

def notFloat(x):
    try:
        float(x)
        return False
    except:
        return True

class more(QtWidgets.QDialog, Ui_More):
    _signal = QtCore.pyqtSignal(dict)
    def __init__(self):
        super(more, self).__init__()
        self.setupUi(self)
        self.retranslateUi(self)
        self.BG_COLOR = saveVal((1, 1, 1))
        self.SLICE_COLOR = saveVal((0, 0, 0))

        self.magnetic_field_checkbox.setChecked(False)

        self.fs_opacity_slider.valueChanged.connect(lambda: self.showValue(self.fs_opacity_slider.value(), self.fs_opacity_line))
        self.line_opacity_slider.valueChanged.connect(lambda: self.showValue(self.line_opacity_slider.value(), self.line_opacity_line))
        self.trait_opacity_slider.valueChanged.connect(lambda: self.showValue(self.trait_opacity_slider.value(), self.trait_opacity_line))
        self.slice_opacity_slider.valueChanged.connect(lambda: self.showValue(self.slice_opacity_slider.value(), self.slice_opacity_line))
        
        self.fs_opacity_line.editingFinished.connect(lambda: self.showText(self.fs_opacity_line.text(), self.fs_opacity_slider, self.fs_opacity_line))
        self.line_opacity_line.editingFinished.connect(lambda: self.showText(self.line_opacity_line.text(), self.line_opacity_slider, self.line_opacity_line))
        self.trait_opacity_line.editingFinished.connect(lambda: self.showText(self.trait_opacity_line.text(), self.trait_opacity_slider, self.trait_opacity_line))
        self.slice_opacity_line.editingFinished.connect(lambda: self.showText(self.slice_opacity_line.text(), self.slice_opacity_slider, self.slice_opacity_line))

        self.fs_opacity_line.setText(str(1.0))
        self.line_opacity_line.setText(str(1.0))
        self.trait_opacity_line.setText(str(1.0))
        self.slice_opacity_line.setText(str(0.4))

        self.background_button.clicked.connect(lambda: self.showColor(self.background_color_frame, self.BG_COLOR))
        self.slice_color_button.clicked.connect(lambda: self.showColor(self.slice_color_frame, self.SLICE_COLOR))


    def showValue(self, text, line):
        line.setText(str(text / 10))
    
    def showText(self, text, slider, line):
        if notFloat(text):
            QtWidgets.QMessageBox.warning(self, "Error", "Input is WRONG!\nPlease Input number between 0 and 1!")
            line.setText(str(slider.value()/10))
        elif (float(text) > 1) | (float(text) < 0):
            QtWidgets.QMessageBox.warning(self, "Error", "Input is WRONG!\nPlease Input number between 0 and 1!")
            line.setText(str(slider.value()/10))

        else:
            slider.setValue(10*float(text))
            line.setText(str(slider.value()/10))
    
    def showColor(self, frm, colorVal):
        col = QtWidgets.QColorDialog.getColor()

        if col.isValid():
            frm.setStyleSheet("QWidget { background-color: %s}" % col.name())
            colorVal.savVal(col.getRgbF())
    
    def closeEvent(self, event):
        para_dict = {}
        para_dict['fs_opacity'] = float(self.fs_opacity_line.text())
        para_dict['line_opacity'] = float(self.line_opacity_line.text())
        para_dict['trait_opacity'] = float(self.trait_opacity_line.text())
        para_dict['slice_opacity'] = float(self.slice_opacity_line.text())
        para_dict['background_color'] = self.BG_COLOR.getVal()[0:3]
        para_dict['slice_color'] = self.SLICE_COLOR.getVal()[0:3]
        para_dict['magnetic_field'] = self.magnetic_field_checkbox.isChecked()
        self._signal.emit(para_dict)
        event.accept()

    def keyPressEvent(self, event) -> None:
        if event.key() == QtCore.Qt.Key_Enter:
            self.close()
        else:
            super().keyPressEvent(event)
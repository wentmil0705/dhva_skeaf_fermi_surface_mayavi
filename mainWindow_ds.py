from Ui_mainWindow import Ui_MainWindow
from pickle import TRUE
from posixpath import dirname
from re import T
from PyQt5 import QtCore, QtGui, QtWidgets
from matplotlib import container
from matplotlib.image import imsave
from mayavi_show import *
import numpy as np
from scipy.interpolate import interpn
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import os
from scipy.interpolate import splprep, splev
from matplot_show import *
import pandas as pd
from comboCheckbox import *
from all_file_dialog_ds import *
from more_ds import *

class saveVal(object):
    
    def __init__(self, v=0.):
        self.val = v
    
    def savVal(self, v):
        self.val = v
    
    def getVal(self):
        return self.val

def isFloat(x):
    try:
        float(x)
        return True
    except:
        return False

class myWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    
    def __init__(self, parent=None):
        super(myWindow, self).__init__(parent)
        self.setupUi(self)
        
        ## basic parameters
        self.BXSF_FILE_PATH = ''
        self.SKEAF_DIRECTORY_PATH = ''
        self.LONG_FILE_PATH = ''
        self.ANG_FILE_PATH = ''
        self.AU_FILE_PATH = ''
        self.CONFIG_FILE_PATH = ''
        self.BCELL = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.EBANDS = []
        self.EFERMI = 0
        self.PREDICTED_FREQS = []
        self.PREDICTED_M = []
        self.PREDICTED_CURV = []
        self.PREDICTED_ORBITTYPE = []
        self.PREDICTED_RUC_COORDS = []
        self.PREDICTED_SLICE = []
        self.PREDICTED_ORBIT = []
        self.PREDICTED_POINT = []
        self.ORBIT_DATA = []
        self.CONVAU2ANG = 0.529177209
        self.ANG_FLAG = False
        self.LONG_DF = pd.DataFrame()
        self.LONG_F = ''
        self.LONG_FLAG = False
        self.LINE_COLOR = saveVal((0., 0., 0.))
        self.INNER_COLOR = saveVal((0.97647059, 0.87843137, 0.5372549))
        self.OUTER_COLOR = saveVal((0.1254902, 0.14509804, 0.35294118))
        self.TRAIT_COLOR = saveVal((0, 1, 0))
        self.MODE = 0    # 0 is simple, 1 is seperate
        self.BZ = 2 # 2 is first bz, 3 is primitive bz
        self.INTERPOL_RATIO = saveVal(1.)
        self.FERMI_ENERGY = saveVal(0.)
        self.EPS = saveVal(1.)
        self.SECTION_V1 = saveVal(0.)
        self.SECTION_V2 = saveVal(0.)
        self.SECTION_V3 = saveVal(0.)
        self.BZ_NUMBER1 = saveVal(1.)
        self.BZ_NUMBER2 = saveVal(1.)
        self.BZ_NUMBER3 = saveVal(1.)
        self.ROTATE1 = saveVal(0.)
        self.ROTATE2 = saveVal(0.)
        self.ROTATE3 = saveVal(0.)
        self.SHOW_SLICE = saveVal(False)
        self.INNER_OUTER = saveVal(False)
        self.LINE_WIDTH = saveVal(1.)
        self.TRAIT_WIDTH = saveVal(1.)
        self.AXES = saveVal(False)
        self.TRAIT = saveVal(False)
        self.SELECT_ROW = saveVal([])
        self.PARA_DICT = {}
        self.THETA = 0.
        self.PHI = 0.
        self.BXSF_FERMI_ENERGY = 0.
        self.LONG_FERMI_ENERGY = 0.
        self.RATIO_CHANGED = False
        self.FULL_TRAIT = saveVal(False)
        self.FULL_FS_FULL_TRAIT = saveVal(False)
        self.LINE = saveVal(True)
        self.FS_OP = 1
        self.SLICE_OP = 0.4
        self.TRAIT_OP = 1
        self.LINE_OP = 1
        self.BG_COLOR = saveVal((1, 1, 1))
        self.SLICE_COLOR = saveVal((0, 0, 0))
        self.MAGNETIC_FIELD = False

        self.predicted_frequencies_tabel_view.itemClicked.connect(self.showPic)

        self.extreme_value_set_parameters_button.clicked.connect(self.showSaveExtremeData)

        self.orbit_outline_set_parameters_button.clicked.connect(self.showSaveOrbitData)

        self.mode_button_group.setId(self.simple_button, 0)
        self.mode_button_group.setId(self.seperate_button, 1)
        self.mode_button_group.buttonClicked[int].connect(self.showButton)
        self.simple_button.setChecked(True)
        self.interpol_ratio_line.setPlaceholderText('1.0')
        self.interpol_ratio_line.textChanged.connect(lambda: self.showRatio(self.interpol_ratio_line.text(), self.INTERPOL_RATIO))
        self.fermi_energy_line.textChanged.connect(lambda: self.showChanged(self.fermi_energy_line.text(), self.FERMI_ENERGY))
        self.first_bz_button.setChecked(True)
        self.bz_button_group.setId(self.first_bz_button, 2)
        self.bz_button_group.setId(self.primitive_bz_button, 3)
        self.bz_button_group.buttonClicked[int].connect(self.showButton)
        self.show_slice_checkbox.setDisabled(False)
        self.show_slice_checkbox.stateChanged.connect(lambda: self.showState(self.show_slice_checkbox.isChecked(), self.SHOW_SLICE))
        self.eps_line.setReadOnly(True)
        self.eps_show_slice_layout.addWidget(self.show_slice_checkbox)
        self.eps_line.setPlaceholderText(str(self.EPS.getVal()))
        self.eps_line.textChanged.connect(lambda: self.showChanged(self.eps_line.text(), self.EPS))
        self.section_v_line1.setReadOnly(True)
        self.section_v_line2.setReadOnly(True)
        self.section_v_line3.setReadOnly(True)
        self.section_v_line1.textChanged.connect(lambda: self.showChanged(self.section_v_line1.text(), self.SECTION_V1))
        self.section_v_line2.textChanged.connect(lambda: self.showChanged(self.section_v_line2.text(), self.SECTION_V2))
        self.section_v_line3.textChanged.connect(lambda: self.showChanged(self.section_v_line3.text(), self.SECTION_V3))
        self.bz_number_line1.textChanged.connect(lambda: self.showChanged(self.bz_number_line1.text(), self.BZ_NUMBER1))
        self.bz_number_line2.textChanged.connect(lambda: self.showChanged(self.bz_number_line2.text(), self.BZ_NUMBER2))
        self.bz_number_line3.textChanged.connect(lambda: self.showChanged(self.bz_number_line3.text(), self.BZ_NUMBER3))
        self.rotate_line1.textChanged.connect(lambda: self.showChanged(self.rotate_line1.text(), self.ROTATE1))
        self.rotate_line2.textChanged.connect(lambda: self.showChanged(self.rotate_line2.text(), self.ROTATE2))
        self.rotate_line3.textChanged.connect(lambda: self.showChanged(self.rotate_line3.text(), self.ROTATE3))
        self.line_checkbox.stateChanged.connect(lambda: self.showState(self.line_checkbox.isChecked(), self.LINE))
        self.line_checkbox.setChecked(True)
        self.line_color_button.clicked.connect(lambda: self.showColor(self.line_color_frame, self.LINE_COLOR))
        self.inner_outer_checkbox.setChecked(False)
        self.inner_color_button.setDisabled(True)
        self.outer_color_button.setDisabled(True)
        self.inner_color_button.clicked.connect(lambda: self.showColor(self.inner_color_frame, self.INNER_COLOR))
        self.outer_color_button.clicked.connect(lambda: self.showColor(self.outer_color_frame, self.OUTER_COLOR))
        self.inner_outer_checkbox.stateChanged.connect(lambda: self.showState(self.inner_outer_checkbox.isChecked(), self.INNER_OUTER))
        self.line_width_line.setPlaceholderText(str(self.LINE_WIDTH.getVal()))
        self.line_width_line.textChanged.connect(lambda: self.showChanged(self.line_width_line.text(), self.LINE_WIDTH))
        self.trait_select_combobox.setDisabled(True)
        self.trait_check_box.setDisabled(True)
        self.trait_check_box.stateChanged.connect(lambda: self.showState(self.trait_check_box.isChecked(), self.TRAIT))
        self.trait_width_line.setPlaceholderText(str(self.TRAIT_WIDTH.getVal()))
        self.trait_width_line.setReadOnly(True)
        self.trait_width_line.textChanged.connect(lambda: self.showChanged(self.trait_width_line.text(), self.TRAIT_WIDTH))
        self.trait_color_button.setDisabled(True)
        self.trait_color_button.clicked.connect(lambda: self.showColor(self.trait_color_frame, self.TRAIT_COLOR))
        self.axes_checkbox.stateChanged.connect(lambda: self.showState(self.axes_checkbox.isChecked(), self.AXES))
        self.full_trait_checkbox.stateChanged.connect(lambda: self.showState(self.full_trait_checkbox.isChecked(), self.FULL_TRAIT))
        self.full_trait_checkbox.setDisabled(True)
        self.full_fs_full_trait_checkbox.setDisabled(True)
        self.full_fs_full_trait_checkbox.stateChanged.connect(lambda: self.showState(self.full_fs_full_trait_checkbox.isChecked(), self.FULL_FS_FULL_TRAIT))
        self.more_dialog = more()
        self.more_button.clicked.connect(self.showMore)
        self.update_button.clicked.connect(self.myUpdate)

        self.actionImport_bxsf_file.triggered.connect(self.showBxsf)
        self.actionImport_results_long_out.triggered.connect(self.showLong)
        self.actionImport_results_orbitoutlines_invAng_out.triggered.connect(self.showAng)
        self.actionImport_results_orbitoutlines_invau_out.triggered.connect(self.showAu)
        self.actionImport_all_results_file.triggered.connect(self.showSkeaf)
        self.actionImport_config_in.triggered.connect(self.showNormalFile)
        



    
    def showBxsf(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '/home')

        if '.bxsf' in fname[0]:
            self.BCELL, self.EBANDS, efermi = self.interpol(str(fname[0]))
            self.BXSF_FERMI_ENERGY = efermi
            self.FERMI_ENERGY.savVal(efermi)
            self.fermi_energy_line.setText(str(self.FERMI_ENERGY.getVal()))
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.bxsfMayavi(self.BCELL, self.EBANDS, self.FERMI_ENERGY.getVal(), self.PARA_DICT)
            self.BXSF_FILE_PATH = fname[0]
            return fname[0]
        
        elif '.bxsf' not in fname[0] and fname[0] != '':
            reply = QtWidgets.QMessageBox.warning(self, 'Import bxsf file', 'Wrong type of bxsf file!\nStill import?', QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            
            if reply == QtWidgets.QMessageBox.Yes:
                # self.mayavi_widget.bxsfMayavi(str(fname[0]))
                self.BXSF_FILE_PATH = fname[0]
                return fname[0]
            
            else:
                return ''
    

    def showSkeaf(self):
        # dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Open directory', '/home')
        
        # if len(dirname) == 0:
        #     return
        
        # # print(dirname)
        # fileList = os.listdir(dirname)
        # # print(fileList)

        # if 'results_long.out' not in fileList:
        #     QtWidgets.QMessageBox.information(self, 'Import all file', 'results_long.out did not exist!')

        # elif 'results_orbitoutlines_invAng.out' not in fileList and 'results_orbitoutlines_invau.out' not in fileList:
        #     QtWidgets.QMessageBox.information(self, 'Import all file', 'None of results_orbitoutlines_invau.out or results_orbitoutlines.invAng.out exists in the directory!')
        
        # else:
        #     self.SKEAF_DIRECTORY_PATH = dirname
        #     return dirname
        
        # return ''
        self.import_all_results_file = all_file_dialog()
        self.import_all_results_file.show()
        self.import_all_results_file._signal.connect(self.getFilePath)
        # self.BXSF_FILE_PATH, self.LONG_FILE_PATH, self.ANG_FILE_PATH, self.AU_FILE_PATH = self.import_all_results_file_ui.closeEvent()

    def showMore(self):
        if self.BXSF_FILE_PATH != '':
            self.more_dialog.fs_opacity_line.setReadOnly(False)
            self.more_dialog.fs_opacity_slider.setDisabled(False)
            if self.LINE.getVal():
                self.more_dialog.line_opacity_line.setReadOnly(False)
                self.more_dialog.line_opacity_slider.setDisabled(False)
            else:
                self.more_dialog.line_opacity_line.setReadOnly(True)
                self.more_dialog.line_opacity_slider.setDisabled(True)
            if self.TRAIT.getVal():
                self.more_dialog.trait_opacity_line.setReadOnly(False)
                self.more_dialog.trait_opacity_slider.setDisabled(False)
            else:
                self.more_dialog.trait_opacity_line.setReadOnly(True)
                self.more_dialog.trait_opacity_slider.setDisabled(True)
            if self.SHOW_SLICE.getVal():
                self.more_dialog.slice_opacity_line.setReadOnly(False)
                self.more_dialog.slice_opacity_slider.setDisabled(False)
                self.more_dialog.slice_color_button.setDisabled(False)
                self.more_dialog.magnetic_field_checkbox.setDisabled(False)
            else:
                self.more_dialog.slice_opacity_line.setReadOnly(True)
                self.more_dialog.slice_opacity_slider.setDisabled(True)
                self.more_dialog.slice_color_button.setDisabled(True)
                self.more_dialog.magnetic_field_checkbox.setDisabled(True)
        else:
            self.more_dialog.fs_opacity_line.setReadOnly(True)
            self.more_dialog.fs_opacity_slider.setDisabled(True)
            self.more_dialog.line_opacity_line.setReadOnly(True)
            self.more_dialog.line_opacity_slider.setDisabled(True) 
            self.more_dialog.trait_opacity_line.setReadOnly(True)
            self.more_dialog.trait_opacity_slider.setDisabled(True)
            self.more_dialog.slice_opacity_line.setReadOnly(True)
            self.more_dialog.slice_opacity_slider.setDisabled(True)
            self.more_dialog.slice_color_button.setDisabled(True)
            self.more_dialog.magnetic_field_checkbox.setDisabled(True)
        
        self.more_dialog.show()
        self.more_dialog._signal.connect(self.getMore)
    
    def getMore(self, tmp_dict):
        self.FS_OP = tmp_dict['fs_opacity']
        self.LINE_OP = tmp_dict['line_opacity']
        self.TRAIT_OP = tmp_dict['trait_opacity']
        self.SLICE_OP = tmp_dict['slice_opacity']
        self.BG_COLOR.savVal(tmp_dict['background_color'])
        self.SLICE_COLOR.savVal(tmp_dict['slice_color'])
        self.MAGNETIC_FIELD = tmp_dict['magnetic_field']

    def getFilePath(self, bxsf, long, ang, au):
        self.BXSF_FILE_PATH = bxsf
        self.LONG_FILE_PATH = long
        self.ANG_FILE_PATH = ang
        self.AU_FILE_PATH = au
        if self.BXSF_FILE_PATH != '':
            self.BCELL, self.EBANDS, efermi = self.interpol(self.BXSF_FILE_PATH)
            self.BXSF_FERMI_ENERGY = efermi
            self.FERMI_ENERGY.savVal(efermi)
            self.fermi_energy_line.setText(str(self.FERMI_ENERGY.getVal()))
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.bxsfMayavi(self.BCELL, self.EBANDS, self.FERMI_ENERGY.getVal(), self.PARA_DICT)
        if self.LONG_FILE_PATH != '':
            self.LONG_DF, self.LONG_F = self.read_long_file(self.LONG_FILE_PATH)
        if self.AU_FILE_PATH != '':
            self.ANG_FLAG = False
            self.read_au_Ang_file(self.AU_FILE_PATH)
            self.trait_check_box.setChecked(False)
            self.trait_check_box.setDisabled(False) 
        if self.ANG_FILE_PATH != '':
            self.ANG_FLAG = True
            self.read_au_Ang_file(self.ANG_FILE_PATH)
            self.trait_check_box.setChecked(False)
            self.trait_check_box.setDisabled(False) 

    

    def showNormalFile(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '/home')

        if fname[0]:
            return fname[0]
        else:
            return ''
    

    def showLong(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '/home')

        if fname[0]:
            self.LONG_DF, self.LONG_F = self.read_long_file(fname[0])
            self.LONG_FILE_PATH = fname[0]
            return fname[0]
        else:
            return ''
    

    def showAu(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '/home')

        if fname[0]:
            self.ANG_FLAG = False
            self.read_au_Ang_file(fname[0])
            self.AU_FILE_PATH = fname[0]
            self.trait_check_box.setChecked(False)
            self.trait_check_box.setDisabled(False)
            return fname[0]
        else:
            return ''
    
    def showAng(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '/home')

        if fname[0]:
            self.ANG_FLAG = True
            self.read_au_Ang_file(fname[0])
            self.ANG_FILE_PATH = fname[0]
            self.trait_check_box.setChecked(False)
            self.trait_check_box.setDisabled(False)
            return fname[0]
        else:
            return ''

    def showPic(self, Item=None):
        
        if Item == None:
            return
        
        else:
            row = Item.row()
            if (self.AU_FILE_PATH != '') | (self.ANG_FILE_PATH != ''):
                self.orbit_outline_pic_widget.drawOrbit(row, self.ORBIT_DATA)
            if (self.LONG_FILE_PATH != ''):
                self.extreme_value_pic_widget.drawExtreme(self.PREDICTED_SLICE[row], self.PREDICTED_ORBIT[row], self.LONG_DF)

    def showColor(self, frm, colorVal):

        col = QtWidgets.QColorDialog.getColor()

        if col.isValid():
            frm.setStyleSheet("QWidget { background-color: %s}" % col.name())
            colorVal.savVal(col.getRgbF())

    
    def showChanged(self, text, textVal):

        # print(text)
        if isFloat(text):
            textVal.savVal(float(text))
        # print(self.INTERPOL_RATIO.getVal())
        # print(self.FERMI_ENERGY.getVal())
    
    def showRatio(self, text, textVal):

        # print(text)
        if isFloat(text):
            if float(text) != self.INTERPOL_RATIO:
                textVal.savVal(float(text))
        # print(self.INTERPOL_RATIO.getVal())
        # print(self.FERMI_ENERGY.getVal())
                self.RATIO_CHANGED = True

        
    def showButton(self, id):

        if id < 2:
            self.MODE = id
            self.fermi_energy_line.setText(str(self.FERMI_ENERGY.getVal()))
            self.eps_line.setText(str(self.EPS.getVal()))
            if id == 0:
                self.trait_check_box.setChecked(False)
                self.bz_number_line1.setReadOnly(False)
                self.bz_number_line2.setReadOnly(False)
                self.bz_number_line3.setReadOnly(False)
                if (self.ANG_FILE_PATH != '') | (self.AU_FILE_PATH != ''):
                    self.trait_check_box.setDisabled(False)
                else:
                    self.trait_check_box.setDisabled(True)
                if self.trait_check_box.isChecked():
                    self.fermi_energy_line.setReadOnly(True)
                else:
                    self.fermi_energy_line.setReadOnly(False)
                self.eps_line.setReadOnly(True)
                self.inner_outer_checkbox.setDisabled(False)
                if self.inner_outer_checkbox.isChecked():
                    self.inner_color_button.setDisabled(False)
                    self.inner_color_button.setDisabled(False)
                else:
                    self.inner_color_button.setDisabled(True)
                    self.outer_color_button.setDisabled(True)
            else:
                self.trait_check_box.setDisabled(False)
                self.trait_check_box.setChecked(True)
                self.fermi_energy_line.setReadOnly(True)
                self.eps_line.setReadOnly(False)
                self.inner_outer_checkbox.setChecked(False)
                self.inner_outer_checkbox.setDisabled(True)
                self.inner_color_button.setDisabled(True)
                self.outer_color_button.setDisabled(True)
                self.bz_number_line1.setText(str(1))
                # self.bz_number_line1.setReadOnly(True)
                self.bz_number_line2.setText(str(1))
                # self.bz_number_line2.setReadOnly(True)                
                self.bz_number_line3.setText(str(1))
                # self.bz_number_line3.setReadOnly(True)       

        else:
            self.BZ = id
        # print(self.MODE)
        # print(self.BZ)
    

    def showState(self, isChecked, stateVal):

        stateVal.savVal(isChecked)
        # print(self.INNER_OUTER.getVal())

        if self.trait_check_box.isChecked():
            self.trait_width_line.setReadOnly(False)
        else:
            self.trait_width_line.setReadOnly(True)
        
        if self.inner_outer_checkbox.isChecked():
            self.inner_color_button.setDisabled(False)
            self.outer_color_button.setDisabled(False)
        else:
            self.inner_color_button.setDisabled(True)
            self.outer_color_button.setDisabled(True)
        
        if self.show_slice_checkbox.isChecked():
            self.section_v_line1.setReadOnly(False)
            self.section_v_line2.setReadOnly(False)
            self.section_v_line3.setReadOnly(False)
        else:
            self.section_v_line1.setReadOnly(True)
            self.section_v_line2.setReadOnly(True)
            self.section_v_line3.setReadOnly(True)
        
        if self.line_checkbox.isChecked():
            self.line_width_line.setReadOnly(False)
            self.line_color_button.setDisabled(False)
        else:
            self.line_width_line.setReadOnly(True)
            self.line_color_button.setDisabled(True)

        if self.trait_check_box.isChecked():
            self.full_trait_checkbox.setDisabled(False)
            self.trait_color_button.setDisabled(False)
            self.trait_select_combobox.setDisabled(False)
            self.full_fs_full_trait_checkbox.setDisabled(False)
            if self.LONG_FILE_PATH != '':
                self.fermi_energy_line.setText(str(self.LONG_FERMI_ENERGY))
            self.fermi_energy_line.setReadOnly(True)
            if self.full_fs_full_trait_checkbox.isChecked():
                self.full_trait_checkbox.setChecked(True)
                self.full_trait_checkbox.setDisabled(True)
                self.bz_number_line1.setText('')
                self.bz_number_line2.setText('')
                self.bz_number_line3.setText('')
                self.BZ_NUMBER1.savVal(1)
                self.BZ_NUMBER2.savVal(1)
                self.BZ_NUMBER3.savVal(1)
                self.bz_number_line1.setReadOnly(True)
                self.bz_number_line2.setReadOnly(True)
                self.bz_number_line3.setReadOnly(True)
            else:
                self.full_trait_checkbox.setDisabled(False)
                self.bz_number_line1.setReadOnly(False)
                self.bz_number_line2.setReadOnly(False)
                self.bz_number_line3.setReadOnly(False)
        else:
            self.full_trait_checkbox.setDisabled(True)
            self.full_fs_full_trait_checkbox.setDisabled(True)
            self.trait_color_button.setDisabled(True)
            self.trait_select_combobox.setDisabled(True)
            self.fermi_energy_line.setReadOnly(False)
        
    def showSaveExtremeData(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Open directory', '/home')

        if len(dirname) != 0:
            self.extreme_value_pic_widget.saveExtremeData(dirname)
    
    def showSaveOrbitData(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Open directory', '/home')

        if len(dirname) != 0:
            self.orbit_outline_pic_widget.saveOrbitData(dirname, self.PREDICTED_FREQS)



    def read_para_dict(self):
        b1 = np.linalg.norm(self.BCELL, axis=1)[0]
        para_dict = {}
        para_dict['mode'] = self.MODE
        para_dict['bz_mode'] = self.BZ
        if self.RATIO_CHANGED & (self.BXSF_FILE_PATH != ''):
            self.BCELL, self.EBANDS, efermi = self.interpol(self.BXSF_FILE_PATH, int(self.INTERPOL_RATIO.getVal() * 80))
            self.RATIO_CHANGED = False
        para_dict['interpol_ratio'] = self.INTERPOL_RATIO.getVal() * 80
        if self.TRAIT.getVal() & (self.LONG_FILE_PATH != ''):
            self.FERMI_ENERGY.savVal(self.LONG_FERMI_ENERGY)
            self.fermi_energy_line.setText(str(self.FERMI_ENERGY.getVal()))
        para_dict['efermi'] = self.FERMI_ENERGY.getVal()
        para_dict['line_color'] = self.LINE_COLOR.getVal()[0:3]
        para_dict['inner_color'] = self.INNER_COLOR.getVal()[0:3]
        para_dict['outer_color'] = self.OUTER_COLOR.getVal()[0:3]
        para_dict['eps'] = self.EPS.getVal() * (b1 / 10 / 2)
        para_dict['section_v'] = np.dot(np.array([self.SECTION_V1.getVal(), self.SECTION_V2.getVal(), self.SECTION_V3.getVal()]), self.BCELL) * 2 * np.pi
        if self.ANG_FLAG == True:
            para_dict['section_v'] = para_dict['section_v'] / self.CONVAU2ANG
        para_dict['bz_number'] = np.array([self.BZ_NUMBER1.getVal(), self.BZ_NUMBER2.getVal(), self.BZ_NUMBER3.getVal()])
        para_dict['rotate'] = np.array([self.ROTATE1.getVal(), self.ROTATE2.getVal(), self.ROTATE3.getVal()])
        para_dict['show_slice'] = self.SHOW_SLICE.getVal()
        para_dict['inner_outer'] = self.INNER_OUTER.getVal()
        para_dict['line_width'] = self.LINE_WIDTH.getVal() * (b1 / 200 / 1.5)
        para_dict['trait_width'] = self.TRAIT_WIDTH.getVal() * (b1 / 200 / 1.5)
        para_dict['trait_color'] = self.TRAIT_COLOR.getVal()[0:3]
        para_dict['axes'] = self.AXES.getVal()
        para_dict['trait'] = self.TRAIT.getVal()
        if (self.ANG_FILE_PATH != '') | (self.AU_FILE_PATH != '') | (self.LONG_FILE_PATH != ''):
            para_dict['mag_h'] = self.ruc2sc()
            para_dict['is_mag_h'] = True
        else:
            para_dict['mag_h'] = np.array([0., 0., 0.])
            para_dict['is_mag_h'] = False
        para_dict['bxsf_file_path'] = self.BXSF_FILE_PATH
        para_dict['orbit_data'] = self.ORBIT_DATA
        para_dict['select_row'] = np.array(self.trait_select_combobox.getCheckItem()) - 1
        para_dict['ang_flag'] = self.ANG_FLAG
        para_dict['freqs'] = self.PREDICTED_FREQS
        para_dict['full_trait'] = self.FULL_TRAIT.getVal()
        para_dict['full_fs_full_trait'] = self.FULL_FS_FULL_TRAIT.getVal()
        para_dict['line'] = self.LINE.getVal()
        para_dict['fs_opacity'] = self.FS_OP
        para_dict['slice_opacity'] = self.SLICE_OP
        para_dict['trait_opacity'] = self.TRAIT_OP
        para_dict['line_opacity'] = self.LINE_OP
        para_dict['slice_color'] = self.SLICE_COLOR.getVal()
        para_dict['background_color'] = self.BG_COLOR.getVal()
        para_dict['magnetic_field'] = self.MAGNETIC_FIELD
        return para_dict

    def myUpdate(self):
        self.PARA_DICT = self.read_para_dict()
        if self.BXSF_FILE_PATH != '':
            self.mayavi_widget.bxsfMayavi(self.BCELL, self.EBANDS, self.FERMI_ENERGY.getVal(), self.PARA_DICT)


    def read_au_Ang_file(self, file_path):
        with open(file_path, 'r') as filereader:
                f = filereader.read()                     
        orbit_data = re.findall(r'kx[\sE\-\.\dkyz\+]+', f)
        orbit_data  = list(map(lambda data: list(map(lambda x: list(map(float, x.split())), data.split('\n')[1:-1])), orbit_data))
        self.ORBIT_DATA = orbit_data
        self.PREDICTED_POINT = list(map(lambda x: int(x), re.findall(r'Points\s+\=\s+(\d+)', f)))
        
        if self.LONG_FILE_PATH == '':
            self.THETA = float(re.findall(r'Theta\(deg\)\s+\=\s+(\d+\.\d+)', f)[0])
            self.PHI = float(re.findall(r'Phi\(deg\)\s+\=\s+(\d+\.\d+)', f)[0])
            self.PREDICTED_FREQS = list(map(lambda x: float(x), re.findall(r'Freq\(kT, average of all copies\)\s+\=\s+(\d+\.\d+)', f)))
            self.PREDICTED_SLICE = list(map(lambda x: int(x), re.findall(r'Slice\s+\=\s+(\d+)', f)))
            self.PREDICTED_ORBIT = list(map(lambda x: int(x), re.findall(r'Orbit\s+\#\s+\=\s+(\d+)', f)))
            
            self.predicted_frequencies_tabel_view.clear()
            self.predicted_frequencies_tabel_view.setRowCount(0)
            self.predicted_frequencies_tabel_view.setColumnCount(4)
            # self.predicted_frequencies_tabel_view.setRowCount(len(self.PREDICTED_FREQS))
            self.predicted_frequencies_tabel_view.setHorizontalHeaderLabels(['Freq(kT)', 'Slice', 'Orbit', 'Points'])
            self.predicted_frequencies_tabel_view.setVerticalHeaderLabels(list(map(lambda x: str(x), np.arange(1, len(self.PREDICTED_FREQS)+1))))
            items = []
            for i in range(len(self.PREDICTED_FREQS)):
                item = []
                item.append(self.PREDICTED_FREQS[i])
                item.append(self.PREDICTED_SLICE[i])
                item.append(self.PREDICTED_ORBIT[i])
                item.append(self.PREDICTED_POINT[i])
                items.append(item)
            for i in range(len(items)):
                item = items[i]
                row = self.predicted_frequencies_tabel_view.rowCount()
                self.predicted_frequencies_tabel_view.insertRow(row)
                for j in range(len(item)):
                    item = QtWidgets.QTableWidgetItem(str(items[i][j]))
                    self.predicted_frequencies_tabel_view.setItem(row,j,item)
            self.predicted_frequencies_tabel_view.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)

            self.trait_select_combobox.clear()
            items = ['All'] + self.PREDICTED_FREQS
            for i in items:
                self.trait_select_combobox.addItem(str(i))

    def interpol(self, file_path, mesh_num = 80):
        self.BZ_NUMBER1.savVal(1)
        self.BZ_NUMBER2.savVal(1)
        self.BZ_NUMBER3.savVal(1)
        self.bz_number_line1.setText(str(1))
        self.bz_number_line2.setText(str(1))
        self.bz_number_line3.setText(str(1))
        with open(file_path, 'r') as filereader:
            bxsf_orig = filereader.read()
        bxsf = bxsf_orig.split('\n')
        bcell = np.array(list(map(lambda x: float(x), ' '.join(bxsf[9:12]).split()))).reshape(3,3)  ## 这个读取不准确，实际上与kVectors_xyz一致

        ## 读取所有采样点的能量值
        fermi_ebands3d_uc = ' '.join(bxsf[13:-3]).split()
        fermi_ebands3d_uc = np.array(list(map(lambda x: float(x), fermi_ebands3d_uc)))

        ebands = fermi_ebands3d_uc

        ## 插值处理过程
        nx = ny = nz = np.linspace(0, 1, mesh_num)
        xi = []
        for x  in nx:
            for y in ny:
                for z in nz:
                    # coord_nowpos = tuple(np.dot(np.array([x, y, z]), bcell))
                    coord_nowpos = [x, y, z]
                    xi.append(coord_nowpos)

        ebands = np.array(list(map(lambda x: float(x), ebands)))
        mesh = np.asarray(list(map(lambda x: int(x), bxsf[7].split())))
        ebands = ebands.reshape(mesh[0], mesh[1], mesh[2])

        points = (np.linspace(0, 1, mesh[0]), np.linspace(0, 1, mesh[1]), np.linspace(0, 1, mesh[2]))
        values = ebands
        xi = np.asarray(xi)
        res = interpn(points, values, xi, method='linear', fill_value=True, bounds_error=False)

        ebands = res

        ## 获取画画费米面所需参数
        ebands = np.array(list(map(lambda x: float(x), ebands)))
        efermi = float(re.findall(r'[\d\.]+', bxsf[1])[0])
        mesh = np.array([mesh_num, mesh_num, mesh_num])
        ebands = ebands.reshape(mesh[0], mesh[1], mesh[2])

        return bcell, ebands, efermi
    
    def read_long_file(self, file_path):
        with open(file_path, 'r') as filereader:
            f = filereader.read()

        self.LONG_FERMI_ENERGY = float(re.findall(r'Fermi energy:\s+(\d+.\d+)', f)[0])

        raw_data = re.findall(r'No orbits|Orbit[\s\d:,a-zA-KM-Z^\(\)\-\_\.\=\*/+]+', f)
        raw_data = raw_data[0:-1]

        df = pd.DataFrame([], columns=['slice', 'orbit', 'freq', 'm', 'avgx', 'stdx', 'avgy', 'stdy', 'maxx', 'minx', 'maxy', 'miny'])
        for i, text in enumerate(raw_data):
            slice = i+1
            if 'No orbits' in text:
                raw_slice = {'slice': [slice]}
            else:
                orbits = list(map(int, re.findall(r'Orbit\s+(\d+)*', text)))
                areas = list(map(float, re.findall(r'Area\s+=\s+(\d+.\d+)*', text)))
                freqs = list(map(float, re.findall(r'Freq.\s+=\s+(\d+.\d+)*', text)))
                ms = list(map(float, re.findall(r'm\*\s+=\s+(\d+.\d+)*', text)))
                avgxs = list(map(float, re.findall(r'Average\s+\(x,y\)\s+=\s+\(\s+(\d+.\d+)*', text)))
                stdxs = list(map(float, re.findall(r'Average\s+\(x,y\)\s+=\s+\(\s+[\d\.]+[\s\+\-/]+(\d+.\d+)*', text)))
                avgys = list(map(float, re.findall(r'Average\s+\(x,y\)\s+=\s+\(\s+[\d\.\s\+\-/]+,\s+(\d+.\d+)*', text)))
                stdys = list(map(float, re.findall(r'Average\s+\(x,y\)\s+=\s+\(\s+[\d\.\s\+\-/]+,\s+[\d.]+\s+\+/\-\s+(\d+.\d+)*', text)))
                maxxs = list(map(float, re.findall(r'maxx\s+=\s+(\d+.\d+)*', text))) 
                minxs = list(map(float, re.findall(r'minx\s+=\s+(\d+.\d+)*', text))) 
                maxys = list(map(float, re.findall(r'maxy\s+=\s+(\d+.\d+)*', text))) 
                minys = list(map(float, re.findall(r'miny\s+=\s+(\d+.\d+)*', text)))
                raw_slice = {'slice': [slice] * len(orbits),
                            'orbit': orbits,
                            'area': areas,
                            'freq': freqs,
                            'm': ms,
                            'avgx': avgxs,
                            'avgy': avgys,
                            'stdx': stdxs,
                            'stdy': stdys,
                            'maxx': maxxs,
                            'maxy': maxys,
                            'minx': minxs,
                            'miny': minys}
            df = df.append(pd.DataFrame(raw_slice))

        self.THETA = float(re.findall(r'Theta\s+\=\s+(\d+.\d+)', f)[0])
        self.PHI = float(re.findall(r'Phi\s+\=\s+(\d+.\d+)', f)[0])
        predicted_area = re.findall(r'Predicted dHvA frequencies:[\sa-zA-Z\d\.\+\-\=\(\)\:\*\,\_\/\^\%]+', f)[0]
        without_predicted_area = f.replace(predicted_area, '')
        self.PREDICTED_FREQS = list(map(lambda x: float(x), re.findall(r'Freq.\s+\=\s+(\-?\d+\.\d+)', predicted_area)))
        self.PREDICTED_M = list(map(lambda x: float(x), re.findall(r'm\*\s+\=\s+(\-?\d+\.\d+)', predicted_area)))
        self.PREDICTED_CURV = list(map(lambda x: float(x), re.findall(r'Curv\.\s+\=\s+(\-?\d+\.\d+)', predicted_area)))
        self.PREDICTED_ORBITTYPE = list(map(lambda x: float(x), re.findall(r'orbit \(1=e,-1=h\):\s+(\-?\d\.\d)', predicted_area)))
        ruc_coords_raw = re.findall(r'RUC avg coords: \(\s+(\-?\d.\d+)\+\/\-\s+[\d\.Na]+,\s+(\-?\d\.\d+)\+\/\-\s+[\d\.Na]+,\s+(\-?\d\.\d+)', predicted_area)
        self.PREDICTED_RUC_COORDS = []
        for i in ruc_coords_raw:
            self.PREDICTED_RUC_COORDS.append(list(map(lambda x: float(x), list(i))))
        self.PREDICTED_SLICE = list(map(lambda x: int(x), re.findall(r'slice\s+(\d+)', predicted_area)))
        self.PREDICTED_ORBIT = list(map(lambda x: int(x), re.findall(r'orbit\s+(\d+)', predicted_area)))
        
        # print(self.PREDICTED_FREQS)
        # print(self.PREDICTED_RUC_COORDS)
        self.predicted_frequencies_tabel_view.clear()
        self.predicted_frequencies_tabel_view.setRowCount(0)
        self.predicted_frequencies_tabel_view.setColumnCount(7)
        # self.predicted_frequencies_tabel_view.setRowCount(len(self.PREDICTED_FREQS))
        self.predicted_frequencies_tabel_view.setHorizontalHeaderLabels(['Freq(kT)', 'm*(m_e)', 'Curv(kT A^2)', 'orbit(1=e,-1=h)', 'RUC avg coords', 'Slice', 'Orbit'])
        self.predicted_frequencies_tabel_view.setVerticalHeaderLabels(list(map(lambda x: str(x), np.arange(1, len(self.PREDICTED_FREQS)+1))))
        items = []
        for i in range(len(self.PREDICTED_FREQS)):
            item = []
            item.append(self.PREDICTED_FREQS[i])
            item.append(self.PREDICTED_M[i])
            item.append(self.PREDICTED_CURV[i])
            item.append(self.PREDICTED_ORBITTYPE[i])
            item.append(self.PREDICTED_RUC_COORDS[i])
            item.append(self.PREDICTED_SLICE[i])
            item.append(self.PREDICTED_ORBIT[i])
            items.append(item)
        for i in range(len(items)):
            item = items[i]
            row = self.predicted_frequencies_tabel_view.rowCount()
            self.predicted_frequencies_tabel_view.insertRow(row)
            for j in range(len(item)):
                item = QtWidgets.QTableWidgetItem(str(items[i][j]))
                self.predicted_frequencies_tabel_view.setItem(row,j,item)
        self.predicted_frequencies_tabel_view.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        
        self.trait_select_combobox.clear()
        items = ['All'] + self.PREDICTED_FREQS
        for i in items:
            self.trait_select_combobox.addItem(str(i))       
        return df, f
    
    def ruc2sc(self):
        phi = np.radians(self.PHI)
        theta = np.radians(self.THETA)
        s = np.sin(phi)
        t = np.cos(phi)
        u = 1 - np.cos(phi)
        v = np.sin(theta)
        w = np.cos(theta)

        sc_to_ruc = np.array([[v*v*u+t, -v*w*u, -w*s],
                            [-v*w*u, w*w*u+t, -v*s],
                            [w*s, v*s, t]
                            ])
        
        kVectors_xyz = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        sc = np.dot(sc_to_ruc, kVectors_xyz)
        return sc[2]

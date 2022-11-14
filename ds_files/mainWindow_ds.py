__author__ = "Guanzhang Liu"
__maintainer__ = "Guanzhang Liu"
__email__ = "MG21340029@smail.nju.edu.cn"
__date__ = "October 23, 2022"

'''
mainWindow_ds.为Ui_mainWindow.py界面对应功能的设计文件
'''

from turtle import bgcolor
from ui_py_files.Ui_mainWindow import Ui_MainWindow
from PyQt5 import QtWidgets
from mayavi_show import *
import numpy as np
from scipy.interpolate import interpn
from matplot_show import *
import pandas as pd
from comboCheckbox import *
from ds_files.all_file_dialog_ds import *
from ds_files.more_ds import *
from share import *
from ds_files.normal_calculate_ds import *
from ds_files.B_direction_ds import *
from ds_files.slice_ds import *



'''myWindow类继承自QMainWindow属性，并承接Ui_MainWindow的界面'''    
class myWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    
    '''初始化，继承父类Ui_mainWindow并进行细节设计'''
    def __init__(self, parent=None):
        super(myWindow, self).__init__(parent)
        self.setupUi(self)

        '''basic parameters'''
        self.CONVAU2ANG = 0.529177209   # 单位转换
        self.BXSF_FILE_PATH = ''    # 费米面文件路径
        self.SKEAF_DIRECTORY_PATH = ''  # skeaf目录路径
        self.LONG_FILE_PATH = ''    # results_long.out文件路径
        self.ANG_FILE_PATH = ''     # results_orbitoutlines_invAng.out文件路径
        self.AU_FILE_PATH = ''      # results_orbitoutlines_invau.out文件路径
        self.CONFIG_FILE_PATH = ''  # config.in文件路径
        self.BCELL = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])    # 扩胞值保存，初始为(1,1,1)
        self.EBANDS = []            # 插值后的能带保存
        self.EFERMI = 0             # 读取的费米能保存
        self.PREDICTED_FREQS = []   # 读取的预测频率列表保存
        self.PREDICTED_M = []       # 读取的预测频率对应质量列表保存
        self.PREDICTED_CURV = []    # 读取的预测频率对应面积列表保存
        self.PREDICTED_ORBITTYPE = []   # 读取的预测频率对应电子/空穴类型概率列表保存
        self.PREDICTED_RUC_COORDS = []  # 读取的预测频率对应轨道中心点坐标列表保存（尚未解决该坐标的位置）
        self.PREDICTED_SLICE = []   # 读取的预测频率对应slice列表保存
        self.PREDICTED_ORBIT = []   # 读取的预测频率对应orbit列表保存
        self.PREDICTED_POINT = []   # 读取的预测频率对应point列表保存
        self.ORBIT_DATA = []        # 读取的预测频率对应轨道坐标列表保存    
        self.ANG_FLAG = False       # 记录是否读取过results_orbitoutlines_invAng.out文件
        self.LONG_DF = pd.DataFrame()   # 读取results_long.out后记录表
        self.LONG_F = ''            # 
        self.LONG_FLAG = False      # 记录是否读取results_long.out文件
        self.LINE_COLOR = saveVal((0., 0., 0.)) # 边界线颜色值保存并设定初始值
        self.INNER_COLOR = saveVal((0.97647059, 0.87843137, 0.5372549)) # 内面颜色值保存并设定初始值 
        self.OUTER_COLOR = saveVal((0.1254902, 0.14509804, 0.35294118)) # 外面颜色值保存并设定初始值
        self.TRAIT_COLOR = saveVal((0, 1, 0))   # 轨道颜色值保存并设定初始值
        self.MODE = 0              # mode模式：0 为 simple, 1 为 seperate
        self.BZ = 2                # 元胞显示方式：2 为 first bz, 3 为 primitive bz
        self.INTERPOL_RATIO = saveVal(1.)   # 插值精度值保存
        self.FERMI_ENERGY = saveVal(0.)     # 费米能值保存
        self.EPS = saveVal(1.)     # 密度聚类精度值保存
        self.SECTION_V1 = saveVal(0.)   # 切片坐标x值保存
        self.SECTION_V2 = saveVal(0.)   # 切片坐标y值保存
        self.SECTION_V3 = saveVal(0.)   # 切片坐标z值保存
        self.BZ_NUMBER1 = saveVal(1.)   # 扩胞x方向值保存
        self.BZ_NUMBER2 = saveVal(1.)   # 扩胞y方向值保存
        self.BZ_NUMBER3 = saveVal(1.)   # 扩胞z方向值保存
        self.ROTATE1 = saveVal(0.)  # 旋转角度x方向值保存
        self.ROTATE2 = saveVal(0.)  # 旋转角度y方向值保存
        self.ROTATE3 = saveVal(0.)  # 旋转角度z方向值保存
        self.SHOW_SLICE = saveVal(False)    # 记录是否显示切片
        self.INNER_OUTER = saveVal(False)   # 记录是否显示内外面区分功能
        self.LINE_WIDTH = saveVal(1.)   # 边界线宽值保存
        self.TRAIT_WIDTH = saveVal(1.)  # 轨道线宽值保存
        self.AXES = saveVal(False)  # 记录是否显示轴标记
        self.TRAIT = saveVal(False) # 记录是否显示轨道
        self.SELECT_ROW = saveVal([])   # 记录所选择显示的轨道
        self.PARA_DICT = {}         # 保存所有参数字典
        self.THETA = 0.             # theta角度
        self.PHI = 0.               # phi角度
        self.BXSF_FERMI_ENERGY = 0. # 记录费米面文件中的费米能
        self.LONG_FERMI_ENERGY = 0. # 记录skeaf文件中的费米能
        self.RATIO_CHANGED = False  # 记录插值精度是否改变
        self.FULL_TRAIT = saveVal(False)    # 记录是否显示完整轨道
        self.FULL_FS_FULL_TRAIT = saveVal(False)    # 记录是否显示完整轨道+费米面
        self.LINE = saveVal(True)   # 记录是否显示边界线
        self.FS_OP = 1              # 费米面透明度值保存
        self.SLICE_OP = 0.4         # 切片透明度值保存
        self.TRAIT_OP = 1           # 轨道透明度值保存
        self.LINE_OP = 1            # 边界线透明度值保存
        self.BG_COLOR = saveVal((1, 1, 1))  # 背景颜色值保存
        self.SLICE_COLOR = saveVal((0, 0, 0))   # 切片颜色值保存
        self.MAGNETIC_FIELD = False # 记录磁场方向是否确定
        self.MESH = np.array([0, 0, 0]) # 保存原始插值点数目
        self.BXSF_FERMI = 0.
        self.xB1 = 0.
        self.xB2 = 0.
        self.xB3 = 1.
        self.gv_mag = False
        self.gv_vector = False
        self.gv_avg = False

        '''链接函数'''
        self.mayavi_widget.camera.connect(self.recCamera)
        # 选择相应频率可显示极值选取与极值形状图
        self.predicted_frequencies_tabel_view.itemClicked.connect(self.showPic)
        # 存储极值选取数据
        self.extreme_value_set_parameters_button.clicked.connect(self.showSaveExtremeData)
        # 存储极值轨道数据
        self.orbit_outline_set_parameters_button.clicked.connect(self.showSaveOrbitData)

        # mode设计
        self.mode_button_group.setId(self.simple_button, 0)
        self.mode_button_group.setId(self.seperate_button, 1)
        self.mode_button_group.buttonClicked[int].connect(self.showMode)
        self.simple_button.setChecked(True) # default为simple

        # 插值设计
        self.interpol_ratio_line.setPlaceholderText('1.0')
        self.interpol_ratio_line.editingFinished.connect(self.showRatio)

        self.fermi_energy_line.textChanged.connect(lambda: self.showChanged(self.fermi_energy_line.text(), self.FERMI_ENERGY))

        self.first_bz_button.setChecked(True)
        self.bz_button_group.setId(self.first_bz_button, 2)
        self.bz_button_group.setId(self.primitive_bz_button, 3)
        self.bz_button_group.buttonClicked[int].connect(self.showBZMode)

        self.show_slice_checkbox.setDisabled(False)
        self.show_slice_checkbox.stateChanged.connect(self.showSlice)

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
        
        self.rotate_line1.editingFinished.connect(self.showCamera)
        self.rotate_line2.editingFinished.connect(self.showCamera)
        self.rotate_line3.editingFinished.connect(self.showCamera)
        
        self.line_checkbox.stateChanged.connect(self.showLine)
        self.line_checkbox.setChecked(True)
        
        self.line_color_button.clicked.connect(self.showLineColor)
        
        self.inner_outer_checkbox.setChecked(False)
        self.inner_color_button.setDisabled(True)
        self.outer_color_button.setDisabled(True)
        self.inner_color_button.clicked.connect(self.showInnerColor)
        self.outer_color_button.clicked.connect(self.showOuterColor)
        self.inner_outer_checkbox.stateChanged.connect(self.showInnerOuter)

        self.line_width_line.setPlaceholderText(str(self.LINE_WIDTH.getVal()))
        self.line_width_line.textChanged.connect(self.showLineChange)

        self.trait_select_combobox.setDisabled(True)
        self.trait_select_combobox.activated.connect(self.show_select_trait)

        self.trait_check_box.setDisabled(True)
        self.trait_check_box.stateChanged.connect(self.showTrait)

        self.trait_width_line.setPlaceholderText(str(self.TRAIT_WIDTH.getVal()))
        self.trait_width_line.setReadOnly(True)
        self.trait_width_line.textChanged.connect(self.showTraitChange)

        self.trait_color_button.setDisabled(True)
        self.trait_color_button.clicked.connect(self.showTraitColor)

        self.axes_checkbox.stateChanged.connect(self.showAxes)

        self.full_trait_checkbox.stateChanged.connect(self.showFullTrait)
        self.full_trait_checkbox.setDisabled(True)

        self.full_fs_full_trait_checkbox.setDisabled(True)
        self.full_fs_full_trait_checkbox.stateChanged.connect(self.showFullFsFullTrait)

        self.more_dialog = more()
        self.more_button.setDisabled(True)
        self.more_button.clicked.connect(self.showMore)

        self.update_button.clicked.connect(self.myUpdate)

        self.actionImport_bxsf_file.triggered.connect(self.showBxsf)
        self.actionImport_results_long_out.triggered.connect(self.showLong)
        self.actionImport_results_orbitoutlines_invAng_out.triggered.connect(self.showAng)
        self.actionImport_results_orbitoutlines_invau_out.triggered.connect(self.showAu)
        self.actionImport_all_results_file.triggered.connect(self.showSkeaf)
        # self.actionImport_config_in.triggered.connect(self.showNormalFile)

        self.import_all_results_file = all_file_dialog()
        self.actionSettng_parameters.triggered.connect(self.showNormalCalc)
        self.actionSaving_results.triggered.connect(self.showSave)
        self.actionB_direction.triggered.connect(self.showBText)
        self.actionB_direction.setDisabled(True)
        self.menuGroup_Velocity.setDisabled(True)
        self.actionvector.setChecked(False)
        self.actioncalculate.triggered.connect(self.showCalcGv)
        self.actionavg.triggered.connect(self.showAvg)
        self.actionmag.triggered.connect(self.showGv)
        self.actionvector.triggered.connect(self.showGv)
        self.actionslice.triggered.connect(self.showWinSlice)
    
    # 显示slice窗口
    def showWinSlice(self):
        self.mySliceWin = mySlice()
        self.mySliceWin.show()
        self.mySliceWin._midSignal.connect(self.synSlice)
    
    # 同步slice窗口截面在实际费米面中的位置
    def synSlice(self, mid_points):
        mid_points = np.array(mid_points)
        self.first_bz_button.setChecked(True)
        self.THETA = self.mySliceWin.THETA
        self.PHI = self.mySliceWin.PHI
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.synSliceMayavi(self.PARA_DICT, mid_points)
    
    # 计算群速度
    def showCalcGv(self):
        self.actionvector.setChecked(False)
        self.actionmag.setChecked(False)
        self.actionavg.setChecked(False)
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.gvCalcMayavi(self.BCELL, self.EBANDS)
    
    # 显示群速度平均值（在mayavi窗口中）
    def showAvg(self):
        self.gv_avg = self.actionavg.isChecked()
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.gvMayavi(self.PARA_DICT)
    
    # 显示群速度图像
    def showGv(self):
        self.gv_mag = self.actionmag.isChecked()
        self.gv_vector = self.actionvector.isChecked()
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.gvMayavi(self.PARA_DICT)

    # 显示磁场方向设定窗口
    def showBText(self):
        self.B_direction = B_direction(self.BCELL)
        self.B_direction.show()
        self.B_direction._angle_signal.connect(self.showB)
        self.B_direction._close_signal.connect(self.endB)
    
    # 结束磁场方向设定窗口
    def endB(self):
        self.mayavi_widget.magSelectEndMayavi()
    
    # 显示磁场设定窗口所设定的磁场方向（蓝色箭头）
    def showB(self, theta, phi):
        self.THETA = theta
        self.PHI = phi
        mag = self.ruc2sc()[2]
        rel_mag = np.dot(mag, np.linalg.inv(self.BCELL))
        self.B_direction.b1_line.setText('{:.2f}'.format(rel_mag[0]))
        self.B_direction.b2_line.setText('{:.2f}'.format(rel_mag[1]))
        self.B_direction.b3_line.setText('{:.2f}'.format(rel_mag[2]))
        self.mayavi_widget.magSelectMayavi(self.BCELL, mag)
    
    # 保存计算结果函数
    def showSave(self):
        saving_files = ['config.in', 'results_long.out', 'results_short.out', 'results_orbitoutlines_invAng.out', 'results_orbitoutlines_invau.out']
        fname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Open directory', DEFAULT_PATH)
        if fname != '':
            try:
                for file in saving_files:
                    print(file)
                    shutil.copy('./'+CALC_FILE+'/'+file, os.path.join(fname, file))
            except:
                    QtWidgets.QMessageBox.critical(self, "Error", "Problems with results file, may be unfinished calculating.")
 
    # 普通计算窗口
    def showNormalCalc(self):
        self.NormalCalc = normal_calculate()
        self.NormalCalc.show()
        self.NormalCalc._bxsfSignal.connect(self.postDeal)
    
    # 普通计算窗口结果导入
    def postDeal(self, f):
        self.readBxsf(f)
        self.readLong('./'+CALC_FILE+'/results_long.out')
        self.readAng('./'+CALC_FILE+'/results_orbitoutlines_invAng.out')
          
    '''
    读取bxsf文件
    '''
    def readBxsf(self, f):
        try:
            self.BCELL, self.EBANDS, efermi = self.interpol(str(f), 1.0)
            self.BXSF_FERMI_ENERGY = efermi
            self.FERMI_ENERGY.savVal(efermi)
            self.fermi_energy_line.setText(str(self.FERMI_ENERGY.getVal()))
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.calcFsMayavi(self.BCELL, self.EBANDS, self.FERMI_ENERGY.getVal(), self.PARA_DICT)
            self.mayavi_widget.seperateCalcFsMayavi()
            self.mayavi_widget.fsMayavi(self.PARA_DICT)
            self.BXSF_FILE_PATH = f
            self.more_button.setDisabled(False)
            self.actionB_direction.setDisabled(False)
            self.menuGroup_Velocity.setDisabled(False)
            self.actionvector.setChecked(False)
            self.actionmag.setChecked(False)
            self.actionavg.setChecked(False)
        except:
            QtWidgets.QMessageBox.critical(self, "Error", "The format or data in your bxsf file is not correct, please check the sample file!")

    def showBxsf(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', DEFAULT_PATH)

        if fname[0] != '':
            self.readBxsf(fname[0])
            return fname[0]
            
        return ''         
        
    
    '''
    一次性打开所有文件
    '''
    def showSkeaf(self):

        self.import_all_results_file.show()
        self.import_all_results_file._signal.connect(self.getFilePath)

    '''
    more子窗口打开函数
    '''
    def showMore(self):
        if self.BXSF_FILE_PATH != '':
            self.more_dialog.fs_opacity_line.setReadOnly(False)
            self.more_dialog.fs_opacity_slider.setDisabled(False)
            if self.line_checkbox.isChecked():
                self.more_dialog.line_opacity_line.setReadOnly(False)
                self.more_dialog.line_opacity_slider.setDisabled(False)
            else:
                self.more_dialog.line_opacity_line.setReadOnly(True)
                self.more_dialog.line_opacity_slider.setDisabled(True)
            if self.trait_check_box.isChecked():
                self.more_dialog.trait_opacity_line.setReadOnly(False)
                self.more_dialog.trait_opacity_slider.setDisabled(False)
            else:
                self.more_dialog.trait_opacity_line.setReadOnly(True)
                self.more_dialog.trait_opacity_slider.setDisabled(True)
            if self.show_slice_checkbox.isChecked():
                self.more_dialog.slice_opacity_line.setReadOnly(False)
                self.more_dialog.slice_opacity_slider.setDisabled(False)
                self.more_dialog.slice_color_button.setDisabled(False)
            else:
                self.more_dialog.slice_opacity_line.setReadOnly(True)
                self.more_dialog.slice_opacity_slider.setDisabled(True)
                self.more_dialog.slice_color_button.setDisabled(True)
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
        
        self.more_dialog.show()
        # 子窗口传递到主窗口信息后保存的变量值
        self.more_dialog._BG_COLOR_signal.connect(self.showBgColor)
        self.more_dialog._fs_opacity_signal.connect(self.showFsOpacity)
        self.more_dialog._line_opacity_signal.connect(self.showLineOpacity)
        self.more_dialog._magnetic_field_signal.connect(self.showMag)
        self.more_dialog._slice_opacity_signal.connect(self.showSliceOpacity)
        self.more_dialog._trait_opacity_signal.connect(self.showTraitOpacity)
        self.more_dialog._SLICE_COLOR_signal.connect(self.showSliceColor)
    
    '''more窗口传递的信号函数'''
    def showBgColor(self, bg_color):
        self.BG_COLOR = saveVal(bg_color)
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.bgMayavi(self.PARA_DICT)
    
    def showFsOpacity(self, fs_Opacity):
        self.FS_OP = fs_Opacity
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.fsMayavi(self.PARA_DICT)
    
    def showLineOpacity(self, line_Opacity):
        self.LINE_OP = line_Opacity
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.fsMayavi(self.PARA_DICT)
    
    def showMag(self, mag):
        self.MAGNETIC_FIELD = mag
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.magMayavi(self.PARA_DICT)
    
    def showTraitOpacity(self, trait_opacity):
        self.TRAIT_OP = trait_opacity
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.traitMayavi(self.PARA_DICT)
    
    def showSliceOpacity(self, slice_opacity):
        self.SLICE_OP = slice_opacity
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.traitMayavi(self.PARA_DICT)

    def showSliceColor(self, slice_color):
        self.SLICE_COLOR.savVal(slice_color)
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.traitMayavi(self.PARA_DICT)

    '''
    传递all_file_dialog子窗口参数
    bxsf, long, ang, au为子窗口通过_signal传递而来的变量
    '''
    def getFilePath(self, bxsf, long, ang, au):
        self.BXSF_FILE_PATH = bxsf
        self.LONG_FILE_PATH = long
        self.ANG_FILE_PATH = ang
        self.AU_FILE_PATH = au

        if self.BXSF_FILE_PATH != '':
            try:
                self.readBxsf(self.BXSF_FILE_PATH)
            except:
                QtWidgets.QMessageBox.critical(self, "Error", "The format or data in your bxsf file is not correct, please check the sample file!")

        if self.LONG_FILE_PATH != '':
            try:
                self.LONG_DF, self.LONG_F = self.read_long_file(self.LONG_FILE_PATH)
            except:
                QtWidgets.QMessageBox.critical(self, "Error", "The format or data in your long file is not correct, please check the sample file!")

        if self.AU_FILE_PATH != '':
            try:
                self.ANG_FLAG = False
                self.read_au_Ang_file(self.AU_FILE_PATH)
                self.trait_check_box.setChecked(False)
                self.trait_check_box.setDisabled(False)
                if self.BXSF_FILE_PATH != '':
                    self.PARA_DICT = self.read_para_dict()
                    self.mayavi_widget.traitCalcMayavi(self.BCELL, self.PARA_DICT)
                    self.mayavi_widget.seperateCalcTraitMayavi(self.PARA_DICT)
            except:
                QtWidgets.QMessageBox.critical(self, "Error", "The format or data in your orbit file is not correct, please check the sample file!")
 
        if self.ANG_FILE_PATH != '':
            try:
                self.ANG_FLAG = True
                self.read_au_Ang_file(self.ANG_FILE_PATH)
                self.trait_check_box.setChecked(False)
                self.trait_check_box.setDisabled(False)
                if self.BXSF_FILE_PATH != '':
                    self.PARA_DICT = self.read_para_dict()
                    self.mayavi_widget.traitCalcMayavi(self.BCELL, self.PARA_DICT)
                    self.mayavi_widget.seperateCalcTraitMayavi(self.PARA_DICT)
            except:
                QtWidgets.QMessageBox.critical(self, "Error", "The format or data in your orbit file is not correct, please check the sample file!")
 
    '''
    打开任意文件（主要用于调试还未实现的功能）
    '''
    def showNormalFile(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', DEFAULT_PATH)

        if fname[0]:
            return fname[0]
        else:
            return ''
    
    '''
    打开long文件
    '''
    def readLong(self, f):
        self.LONG_DF, self.LONG_F = self.read_long_file(f)
        self.LONG_FILE_PATH = f

    def showLong(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', DEFAULT_PATH)

        if fname[0]:
            try:
                self.readLong(fname[0])
                return fname[0]
            except:
                QtWidgets.QMessageBox.critical(self, "Error", "The format or data in your long file is not correct, please check the sample file!")

        return ''

    '''
    打开results_orbitoutlines_invau.out文件
    '''
    def readAu(self, f):
        self.ANG_FLAG = False
        self.read_au_Ang_file(f)
        self.AU_FILE_PATH = f
        self.trait_check_box.setChecked(False)
        self.trait_check_box.setDisabled(False)
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.traitCalcMayavi(self.BCELL, self.PARA_DICT)
            self.mayavi_widget.seperateCalcTraitMayavi(self.PARA_DICT)   

    def showAu(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', DEFAULT_PATH)

        if fname[0]:
            try:
                self.readAu(fname[0])
                return fname[0]
            except:
                QtWidgets.QMessageBox.critical(self, "Error", "The format or data in your orbit file is not correct, please check the sample file!")
        
        return ''
    
    '''
    打开results_orbitoutlines_invAng.out文件
    '''
    def readAng(self, f):
        self.ANG_FLAG = True
        self.read_au_Ang_file(f)
        self.ANG_FILE_PATH = f
        self.trait_check_box.setChecked(False)
        self.trait_check_box.setDisabled(False)
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.traitCalcMayavi(self.BCELL, self.PARA_DICT)
            self.mayavi_widget.seperateCalcTraitMayavi(self.PARA_DICT)                 

    def showAng(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', DEFAULT_PATH)

        if fname[0]:
            try:
                self.readAng(fname[0])
                return fname[0]
            except:
                QtWidgets.QMessageBox.critical(self, "Error", "The format or data in your orbit file is not correct, please check the sample file!")
        
        return ''

    '''
    根据选择显示matplotlib组件对应图片
    Item：调用时自动传递的所选择行数的变量
    '''
    def showPic(self, Item=None):

        if Item == None:
            return
        
        else:
            row = Item.row()
            if (self.AU_FILE_PATH != '') | (self.ANG_FILE_PATH != ''):
                self.orbit_outline_pic_widget.drawOrbit(row, self.ORBIT_DATA)

            if (self.LONG_FILE_PATH != ''):
                self.extreme_value_pic_widget.drawExtreme(self.PREDICTED_SLICE[row], self.PREDICTED_ORBIT[row], self.LONG_DF)

    '''
    颜色选取函数
    '''
    def showLineColor(self):
        self.showColor(self.line_color_frame, self.LINE_COLOR)
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.edgeMayavi(self.PARA_DICT)

    def showInnerColor(self):
        self.showColor(self.inner_color_frame, self.INNER_COLOR)
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.fsMayavi(self.PARA_DICT)
    
    def showOuterColor(self):
        self.showColor(self.outer_color_frame, self.OUTER_COLOR)
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.fsMayavi(self.PARA_DICT)

    def showTraitColor(self):
        self.showColor(self.trait_color_frame, self.TRAIT_COLOR)
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.traitMayavi(self.PARA_DICT)

    def showColor(self, frm, colorVal):

        col = QtWidgets.QColorDialog.getColor()

        if col.isValid():
            frm.setStyleSheet("QWidget { background-color: %s}" % col.name())
            colorVal.savVal(col.getRgbF())

    '''
    值改变函数
    '''
    def showLineChange(self):
        self.showChanged(self.line_width_line.text(), self.LINE_WIDTH)
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.edgeMayavi(self.PARA_DICT)
    
    def showTraitChange(self):
        self.showChanged(self.trait_width_line.text(), self.TRAIT_WIDTH)
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.traitMayavi(self.PARA_DICT)
    
    def showChanged(self, text, textVal):

        if isFloat(text):
            textVal.savVal(float(text))
    
    '''
    插值改变与保存函数
    '''
    def showRatio(self):

        text = self.interpol_ratio_line.text()
        if isFloat(text):
            if float(text) != self.INTERPOL_RATIO:
                self.INTERPOL_RATIO.savVal(float(text))
                self.RATIO_CHANGED = True
                if self.BXSF_FILE_PATH != '':
                    self.PARA_DICT = self.read_para_dict()
                    self.mayavi_widget.calcFsMayavi(self.BCELL, self.EBANDS, self.FERMI_ENERGY.getVal(), self.PARA_DICT)
                    self.mayavi_widget.seperateCalcFsMayavi()
                    self.mayavi_widget.fsMayavi(self.PARA_DICT)


    '''
    参数基础button的显示与记录函数
    id：所选择的button
    '''
    def showBZMode(self, id):
        self.BZ = id
        self.PARA_DICT = self.read_para_dict()
        if self.BXSF_FILE_PATH != '':
            self.mayavi_widget.fsMayavi(self.PARA_DICT)
    
    def showMode(self, id):
        self.MODE = id
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.fsMayavi(self.PARA_DICT)

    '''
    参数基础checkbox的显示与记录函数
    '''
    def showTrait(self):
        self.TRAIT.savVal(self.trait_check_box.isChecked())
        if self.trait_check_box.isChecked():
            self.trait_width_line.setReadOnly(False)
            self.full_trait_checkbox.setDisabled(False)
            self.trait_color_button.setDisabled(False)
            self.trait_select_combobox.setDisabled(False)
            self.full_fs_full_trait_checkbox.setDisabled(False)
            if self.LONG_FILE_PATH != '':
                self.fermi_energy_line.setText(str(self.LONG_FERMI_ENERGY))
            self.fermi_energy_line.setReadOnly(True)
        else:
            self.full_trait_checkbox.setDisabled(True)
            self.full_fs_full_trait_checkbox.setDisabled(True)
            self.trait_color_button.setDisabled(True)
            self.trait_select_combobox.setDisabled(True)
            self.fermi_energy_line.setReadOnly(False)
    
    # 显示完整的费米面与轨道函数
    def showFullFsFullTrait(self):
        self.FULL_FS_FULL_TRAIT.savVal(self.full_fs_full_trait_checkbox.isChecked())
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
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.traitCalcMayavi(self.BCELL, self.PARA_DICT)
            self.PARA_DICT['line_color'] = (0,0,1)
            self.mayavi_widget.calcFsMayavi(self.BCELL, self.EBANDS, self.FERMI_ENERGY.getVal(), self.PARA_DICT)            
            self.mayavi_widget.seperateCalcFsMayavi()
            self.mayavi_widget.seperateCalcTraitMayavi(self.PARA_DICT)
            self.mayavi_widget.fsMayavi(self.PARA_DICT)            

    # 显示完整轨道函数
    def showFullTrait(self):
        self.FULL_TRAIT.savVal(self.full_trait_checkbox.isChecked())
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.traitMayavi(self.PARA_DICT)
    
    # 显示费米面内外颜色
    def showInnerOuter(self):
        self.INNER_OUTER.savVal(self.inner_outer_checkbox.isChecked())
        # 区分内外功能开启才可编辑选择颜色功能
        if self.inner_outer_checkbox.isChecked():
            self.inner_color_button.setDisabled(False)
            self.outer_color_button.setDisabled(False)
        else:
            self.inner_color_button.setDisabled(True)
            self.outer_color_button.setDisabled(True)
        if self.BXSF_FILE_PATH != '':
            self.PARA_DICT = self.read_para_dict()
            self.mayavi_widget.fsMayavi(self.PARA_DICT)

    # 显示切片位置
    def showSlice(self):
        self.SHOW_SLICE.savVal(self.show_slice_checkbox.isChecked())
        # 显示切片功能开启才可输入坐标
        if self.show_slice_checkbox.isChecked():
            self.section_v_line1.setReadOnly(False)
            self.section_v_line2.setReadOnly(False)
            self.section_v_line3.setReadOnly(False)
        else:
            self.section_v_line1.setReadOnly(True)
            self.section_v_line2.setReadOnly(True)
            self.section_v_line3.setReadOnly(True)
    
    # 显示边界
    def showLine(self):
        self.LINE.savVal(self.line_checkbox.isChecked())
         # 边界线显示功能开启才能编辑边界线特征
        if self.line_checkbox.isChecked():
            self.line_width_line.setReadOnly(False)
            self.line_color_button.setDisabled(False)
        else:
            self.line_width_line.setReadOnly(True)
            self.line_color_button.setDisabled(True)
        self.PARA_DICT = self.read_para_dict()
        if self.BXSF_FILE_PATH != '':
            self.mayavi_widget.edgeMayavi(self.PARA_DICT)
        
    # 显示轴
    def showAxes(self):
        self.AXES.savVal(self.axes_checkbox.isChecked())
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.axesMayavi(self.PARA_DICT)

    # 显示目前的图像位置信息
    def recCamera(self, v1, v2, v3):
        self.rotate_line1.setText('{:.2f}'.format(v1))
        self.rotate_line2.setText('{:.2f}'.format(v2))
        self.rotate_line3.setText('{:.2f}'.format(v3))

    # 更新所设定的图像位置信息
    def showCamera(self):
        self.showChanged(self.rotate_line1.text(), self.ROTATE1)
        self.showChanged(self.rotate_line2.text(), self.ROTATE2)
        self.showChanged(self.rotate_line3.text(), self.ROTATE3)
        self.mayavi_widget.updateCameraMayavi(self.ROTATE1.getVal(), self.ROTATE2.getVal(), self.ROTATE3.getVal())
    
    # 显示所选轨道
    def show_select_trait(self):
        self.PARA_DICT = self.read_para_dict()
        self.mayavi_widget.traitMayavi(self.PARA_DICT)

    '''
    输出极值选取数据
    '''
    def showSaveExtremeData(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Open directory', DEFAULT_PATH)

        if len(dirname) != 0:
            self.extreme_value_pic_widget.saveExtremeData(dirname)
    
    '''
    输出极值轨道数据
    '''
    def showSaveOrbitData(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Open directory', DEFAULT_PATH)

        if len(dirname) != 0:
            self.orbit_outline_pic_widget.saveOrbitData(dirname, self.PREDICTED_FREQS)

    '''
    传递参数函数
    '''
    def read_para_dict(self):
        b1, b2, b3 = np.linalg.norm(self.BCELL, axis=1)
        para_dict = {}
        para_dict['mode'] = self.mode_button_group.checkedId()
        para_dict['bz_mode'] = self.bz_button_group.checkedId()
        if self.RATIO_CHANGED & (self.BXSF_FILE_PATH != ''):
            self.BCELL, self.EBANDS, self.BXSF_FERMI_ENERGY = self.interpol(self.BXSF_FILE_PATH, self.INTERPOL_RATIO.getVal())
            self.RATIO_CHANGED = False
        para_dict['interpol_ratio'] = self.INTERPOL_RATIO.getVal() * self.MESH
        if self.trait_check_box.isChecked() & (self.LONG_FILE_PATH != ''):
            self.fermi_energy_line.setText(str(self.LONG_FERMI_ENERGY))
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
        para_dict['show_slice'] = self.show_slice_checkbox.isChecked()
        para_dict['inner_outer'] = self.inner_outer_checkbox.isChecked()
        para_dict['line_width'] = self.LINE_WIDTH.getVal() * (b1 / 200 / 1.5)
        para_dict['trait_width'] = self.TRAIT_WIDTH.getVal() * (b1 / 200 / 1.5)
        para_dict['trait_color'] = self.TRAIT_COLOR.getVal()[0:3]
        para_dict['axes'] = self.axes_checkbox.isChecked()
        para_dict['trait'] = self.trait_check_box.isChecked()
        if (self.ANG_FILE_PATH != '') | (self.AU_FILE_PATH != '') | (self.LONG_FILE_PATH != ''):
            # 使用默认磁场方向
            para_dict['mag_h'] = self.ruc2sc()[2]
            para_dict['is_mag_h'] = True
        else:
            para_dict['mag_h'] = self.ruc2sc()[2]
            para_dict['is_mag_h'] = False
        para_dict['bxsf_file_path'] = self.BXSF_FILE_PATH
        para_dict['orbit_data'] = self.ORBIT_DATA
        para_dict['select_row'] = np.array(self.trait_select_combobox.getCheckItem()) - 1
        para_dict['ang_flag'] = self.ANG_FLAG
        para_dict['freqs'] = self.PREDICTED_FREQS
        para_dict['full_trait'] = self.full_trait_checkbox.isChecked()
        para_dict['full_fs_full_trait'] = self.full_fs_full_trait_checkbox.isChecked()
        para_dict['line'] = self.line_checkbox.isChecked()
        para_dict['fs_opacity'] = self.FS_OP
        para_dict['slice_opacity'] = self.SLICE_OP
        para_dict['trait_opacity'] = self.TRAIT_OP
        para_dict['line_opacity'] = self.LINE_OP
        para_dict['slice_color'] = self.SLICE_COLOR.getVal()
        para_dict['background_color'] = self.BG_COLOR.getVal()
        para_dict['magnetic_field'] = self.MAGNETIC_FIELD
        para_dict['vector'] = self.gv_vector
        para_dict['mag'] = self.gv_mag
        para_dict['avg_mag'] = self.gv_avg
        para_dict['sc'] = self.ruc2sc()
        para_dict['bcell'] = self.BCELL
        return para_dict

    '''
    update更新函数（针对于扩胞等）
    '''
    def myUpdate(self):
        self.PARA_DICT = self.read_para_dict()
        # 调用mayavi组件函数进行处理
        if self.BXSF_FILE_PATH != '':
            self.mayavi_widget.calcFsMayavi(self.BCELL, self.EBANDS, self.FERMI_ENERGY.getVal(), self.PARA_DICT)
            self.mayavi_widget.seperateCalcFsMayavi()
            if (self.AU_FILE_PATH != '') or (self.ANG_FILE_PATH != ''):
                self.mayavi_widget.traitCalcMayavi(self.BCELL, self.PARA_DICT)
                self.mayavi_widget.seperateCalcTraitMayavi(self.PARA_DICT)
            self.mayavi_widget.fsMayavi(self.PARA_DICT)


    '''
    读取au/Ang轨道坐标文件
    file_path：文件路径
    '''
    def read_au_Ang_file(self, file_path):
        with open(file_path, 'r') as filereader:
            f = filereader.read()

        # 处理没有极值频率的情况（此时文件为空）
        try:                     
            orbit_data = re.findall(r'kx[\sE\-\.\dkyz\+]+', f)
            orbit_data  = list(map(lambda data: list(map(lambda x: list(map(float, x.split())), data.split('\n')[1:-1])), orbit_data))
            self.ORBIT_DATA = orbit_data
            self.PREDICTED_POINT = list(map(lambda x: int(x), re.findall(r'Points\s+\=\s+(\d+)', f)))
            orbit_pred_freqs = list(map(lambda x: float(x), re.findall(r'Freq\(kT, average of all copies\)\s+\=\s+(\d+\.\d+)', f)))
            
            if self.LONG_FILE_PATH == '':
                self.PREDICTED_FREQS = orbit_pred_freqs
                self.THETA = float(re.findall(r'Theta\(deg\)\s+\=\s+(\d+\.\d+)', f)[0])
                self.PHI = float(re.findall(r'Phi\(deg\)\s+\=\s+(\d+\.\d+)', f)[0])
                self.PREDICTED_SLICE = list(map(lambda x: int(x), re.findall(r'Slice\s+\=\s+(\d+)', f)))
                self.PREDICTED_ORBIT = list(map(lambda x: int(x), re.findall(r'Orbit\s+\#\s+\=\s+(\d+)', f)))
                
                self.predicted_frequencies_tabel_view.clear()
                self.predicted_frequencies_tabel_view.setRowCount(0)
                self.predicted_frequencies_tabel_view.setColumnCount(4)
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
            # long 与 orbit不匹配的情况
            else:
                if self.PREDICTED_FREQS != orbit_pred_freqs:
                    self.PREDICTED_FREQS = orbit_pred_freqs
                    self.LONG_FILE_PATH = ''
                    QtWidgets.QMessageBox.critical(self, "Error", "long file and orbit file do not match, please reload the long file!")
        except:
            QtWidgets.QMessageBox.warning(self, "Warning", "au/Ang file has problems, maybe because of empty results.")

    '''
    插值函数
    file_path：bxsf文件路径
    mesh_num：原始的'''
    def interpol(self, file_path, mesh_ratio=1.0):
        with open(file_path, 'r') as filereader:
            bxsf_orig = filereader.read()
        self.MESH = np.fromiter(map(lambda x: int(x), list(re.findall('BANDGRID_3D_BANDS\n \d+[\s]+(\d+) (\d+) (\d+)', bxsf_orig)[0])), dtype=int)
        mesh_num = np.fromiter(self.MESH * mesh_ratio, dtype=int)
        bxsf = bxsf_orig.split('\n')
        bcell = np.array(list(map(lambda x: float(x), ' '.join(bxsf[9:12]).split()))).reshape(3,3)  ## 这个读取不准确，实际上与kVectors_xyz一致

        ## 读取所有采样点的能量值
        fermi_ebands3d_uc = ' '.join(bxsf[13:-3]).split()
        fermi_ebands3d_uc = np.array(list(map(lambda x: float(x), fermi_ebands3d_uc)))

        ebands = fermi_ebands3d_uc

        ## 插值处理过程
        ## xi为插值后的坐标位置
        nx = np.linspace(0, 1, mesh_num[0])
        ny = np.linspace(0, 1, mesh_num[1])
        nz = np.linspace(0, 1, mesh_num[2])
        X,Y,Z = np.meshgrid(nx,ny,nz)
        xi = np.concatenate((X[:,:,:,None], Y[:,:,:,None], Z[:,:,:,None]),axis=-1).reshape(-1,3)

        ebands = np.array(list(map(lambda x: float(x), ebands)))
        ebands = ebands.reshape(self.MESH[0], self.MESH[1], self.MESH[2])

        points = (np.linspace(0, 1, self.MESH[0]), np.linspace(0, 1, self.MESH[1]), np.linspace(0, 1, self.MESH[2]))
        values = ebands
        # xi = np.asarray(xi)
        res = interpn(points, values, xi, method='linear', fill_value=True, bounds_error=False) 
        # 此处使用的插值方法为linear，可以调整其他方法，但计算速度会不一样，详情参考scipy.interpolate.interpn的method介绍

        ebands = res

        ## 获取画画费米面所需参数
        ebands = np.array(list(map(lambda x: float(x), ebands)))
        efermi = float(re.findall(r'[\d\.]+', bxsf[1])[0])
        ebands = ebands.reshape(mesh_num[0], mesh_num[1], mesh_num[2])

        return bcell, ebands, efermi
    
    '''
    读取results_long.out文件的信息
    '''
    def read_long_file(self, file_path):
        with open(file_path, 'r') as filereader:
            f = filereader.read()

        self.LONG_FERMI_ENERGY = float(re.findall(r'Fermi energy:\s+(\d+.\d+)', f)[0])

        ## 处理没有orbit的slice
        raw_data = re.findall(r'No orbits|Orbit[\s\d:,a-zA-KM-Z^\(\)\-\_\.\=\*/+]+', f)
        raw_data = raw_data[0:-1]

        print('done')
        ## 处理全部slice
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
            df = pd.concat([df, pd.DataFrame(raw_slice)])

        ## 记录磁场方向
        self.THETA = float(re.findall(r'Theta\s+\=\s+(\d+.\d+)', f)[0])
        self.PHI = float(re.findall(r'Phi\s+\=\s+(\d+.\d+)', f)[0])
        
        try:
            ## 记录预测频率
            predicted_area = re.findall(r'Predicted dHvA frequencies:[\sa-zA-Z\d\.\+\-\=\(\)\:\*\,\_\/\^\%]+', f)[0]
            ## 处理long与orbit不匹配的情况
            long_pre_freqs = list(map(lambda x: float(x), re.findall(r'Freq.\s+\=\s+(\-?\d+\.\d+)', predicted_area)))
            if (self.PREDICTED_FREQS != []) & (self.PREDICTED_FREQS != long_pre_freqs):
                self.PREDICTED_FREQS = long_pre_freqs
                self.ANG_FILE_PATH = ''
                self.AU_FILE_PATH = ''
                QtWidgets.QMessageBox.critical(self, "Error", "long file and orbit file do not match, please reload the orbit file!")
            else:
                self.PREDICTED_FREQS = long_pre_freqs
                self.PREDICTED_M = list(map(lambda x: float(x), re.findall(r'm\*\s+\=\s+(\-?\d+\.\d+)', predicted_area)))
                self.PREDICTED_CURV = list(map(lambda x: float(x), re.findall(r'Curv\.\s+\=\s+(\-?\d+\.\d+)', predicted_area)))
                self.PREDICTED_ORBITTYPE = list(map(lambda x: float(x), re.findall(r'orbit \(1=e,-1=h\):\s+(\-?\d\.\d)', predicted_area)))
                ruc_coords_raw = re.findall(r'RUC avg coords: \(\s+(\-?\d.\d+)\+\/\-\s+[\d\.Na]+,\s+(\-?\d\.\d+)\+\/\-\s+[\d\.Na]+,\s+(\-?\d\.\d+)', predicted_area)
                self.PREDICTED_RUC_COORDS = []
                for i in ruc_coords_raw:
                    self.PREDICTED_RUC_COORDS.append(list(map(lambda x: float(x), list(i))))
                self.PREDICTED_SLICE = list(map(lambda x: int(x), re.findall(r'slice\s+(\d+)', predicted_area)))
                self.PREDICTED_ORBIT = list(map(lambda x: int(x), re.findall(r'orbit\s+(\d+)', predicted_area)))
                
                self.predicted_frequencies_tabel_view.clear()
                self.predicted_frequencies_tabel_view.setRowCount(0)
                self.predicted_frequencies_tabel_view.setColumnCount(7)
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
        except:
            QtWidgets.QMessageBox.warning(self, "Warning", "long file has problems, maybe because of empty results.")       
        return df, f
    
    '''
    计算磁场方向
    '''
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
        return sc

__author__ = "Guanzhang Liu"
__maintainer__ = "Guanzhang Liu"
__email__ = "MG21340029@smail.nju.edu.cn"
__date__ = "October 23, 2022"

'''
normal_calculate_ds.py为ui_normal_calculate.py界面对应功能的设计文件
'''

from ui_py_files.ui_normal_calculate import *
from PyQt5 import QtWidgets
from PyQt5.QtCore import QProcess, pyqtSignal
from time import sleep
from share import *
import re
import shutil
import os

'''normal_calculate类继承自QWidget属性，并承接Ui_NormalCalculate的界面'''    
class normal_calculate(QtWidgets.QWidget, Ui_NormalCalculate):
    _bxsfSignal = pyqtSignal(str)   # 传递导入文件信号

    def __init__(self, parent=None):
        super(normal_calculate, self).__init__(parent)
        self.setupUi(self)

        self.p = None
        self.start_button.pressed.connect(self.start_process)
        self.run_content.setReadOnly(True)
        self.bxsf_button.pressed.connect(self.showBxsf)
        self.inter_line.setText(str(120))
        self.theta_line.setText(str(0))
        self.phi_line.setText(str(0))
        self.min_ext_fs_line.setText(str(0.01))
        self.max_fc_diff_line.setText(str(0.01))
        self.max_dist_line.setText(str(0.8))
        self.fermi_energy = 0.
        self.interpolation = 120
        self.theta = 0.
        self.phi = 0.
        self.min_ext_fs = 0.01
        self.max_fac_diff = 0.01
        self.max_dist = 0.8
        self.bxsf_label.setText('')

    # 将终端信息更新到窗口
    def message(self, s):
        self.run_content.append(s)
    
    # 检查参数设定是否正确
    def checkParam(self, s, param_name):
        try:
            float(s)
            return s
        except:
            QtWidgets.QMessageBox.critical(self, "Error", "Input for "+param_name+" is not correct!")
            raise ValueError('Input should be number!')

    # 检查计算环境是否合理并写入config.in文件
    def checkConfig(self):
        bxsf_file = [x for x in os.listdir() if x != 'a.out']
        if bxsf_file != []:
            for i in bxsf_file:
                os.remove(i)
        shutil.copy(self.bxsf_label.text(), './'+self.bxsf_label.text().split('/')[-1])
        self.fermi_energy = self.checkParam(self.fermi_line.text(), 'fermi_energy')
        self.interpolation = self.checkParam(self.inter_line.text(), 'interpolation')
        self.theta = self.checkParam(self.theta_line.text(), 'theta')
        self.phi = self.checkParam(self.phi_line.text(), 'phi')
        self.min_ext_fs = self.checkParam(self.min_ext_fs_line.text(), 'Minumum extremal Fs freq.(kT)')
        self.max_fac_diff = self.checkParam(self.max_fc_diff_line.text(), 'Maximum factional diff. between orbit freqs.')
        self.max_dist = self.checkParam(self.max_dist_line.text(), 'Maximum distance between orbit avg. coords.')
        with open('./config.in', 'w') as filewriter:
            filewriter.writelines(self.bxsf_label.text().split('/')[-1]+'\n')
            filewriter.writelines(self.fermi_energy+'\n')
            filewriter.writelines(self.interpolation+'\n')
            filewriter.writelines(self.theta+'\n')
            filewriter.writelines(self.phi+'\nn\n')
            filewriter.writelines(self.min_ext_fs+'\n')
            filewriter.writelines(self.max_fac_diff+'\n')
            filewriter.writelines(self.max_dist+'\nn\n')
            filewriter.writelines(self.theta+'\n'+self.theta+'\n'+self.phi+'\n'+self.phi+'\n'+'1')

    # 开始运行
    def start_process(self):
        try:
            os.chdir('./'+CALC_FILE)
            self.checkConfig()
            if self.p is None:
                self.message("Start Calculating")
                self.p = QProcess()
                self.p.readyReadStandardOutput.connect(self.handle_stdout)
                self.p.readyReadStandardError.connect(self.handle_stderr)
                self.p.stateChanged.connect(self.handle_state)
                self.p.finished.connect(self.process_finished)
                self.p.start("./a.out", [])
                sleep(1)
                self.p.write(u"y\n".encode('utf-8'))
                sleep(1)
                self.p.write(u"y\n".encode('utf-8'))
            os.chdir('..')
        except:
            # pass
            QtWidgets.QMessageBox.critical(self, "Error", "Please recheck your parameters!")
 
    # 运行的错误信息
    def handle_stderr(self):
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.message(stderr)
    
    # 运行的终端信息更新
    def handle_stdout(self):
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.message(stdout)
    
    # 运行状态更新
    def handle_state(self, state):
        states = {
            QProcess.NotRunning: "Not running",
            QProcess.Starting: "Starting",
            QProcess.Running: "Running"
        }
        state_name = states[state]
        self.message(f"State changed: {state_name}")
    
    # 结束运行提醒，并询问是否导入结果文件分析
    def process_finished(self):
        self.message("Process finished")
        self.p = None
        reply = QtWidgets.QMessageBox.question(self, '', 'Do you want to import result files to software?',QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if reply == QtWidgets.QMessageBox.Yes:
            self._bxsfSignal.emit(self.bxsf_label.text())
        else:
            pass
    
    # 选择bxsf文件
    def showBxsf(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', DEFAULT_PATH)

        if fname[0] != '':
            try:
                with open(fname[0], 'r') as filereader:
                    bxsf_orig = filereader.read()
                bxsf = bxsf_orig.split('\n')
                efermi = re.findall(r'[\d\.]+', bxsf[1])[0]
                self.fermi_line.setText(efermi)
                self.bxsf_label.setText(fname[0])
            except:
                QtWidgets.QMessageBox.critical(self, "Error", "The format or data in your bxsf file is not correct, please check the sample file!")


    
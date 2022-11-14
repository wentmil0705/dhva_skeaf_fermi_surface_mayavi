#!/usr/bin
__author__ = "Guanzhang Liu"
__maintainer__ = "Guanzhang Liu"
__email__ = "MG21340029@smail.nju.edu.cn"
__date__ = "October 23, 2022"

'''本文件为打开软件所用，只需要在配置完成环境后，终端cd到demo.py所在位置，并运行 python demo.py即可'''

from ds_files.mainWindow_ds import myWindow
from PyQt5.QtWidgets import QApplication, QMainWindow
import sys


if __name__ == '__main__':

    app = QApplication([])
    skeaf = myWindow()
    skeaf.show()
    sys.exit(app.exec_())  
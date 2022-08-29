#!/usr/bin

from mainWindow_ds import myWindow
from PyQt5 import sip
from PyQt5.QtWidgets import QApplication, QMainWindow
import sys


if __name__ == '__main__':

    app = QApplication([])
    skeaf = myWindow()
    skeaf.show()
    sys.exit(app.exec_())  
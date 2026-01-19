import sys
import random
from PySide6 import QtGui, QtCore, QtGui, QtWidgets

icon_path = "M57.png"

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName("ASTROPHOTO")
    app.setWindowIcon(QtGui.QIcon(icon_path))

    # parser = QtCore.QCommandLineParser()
    # parser.addPositionalArgument("", "Project directory", "[path]")
    # parser.process(app)

    app.exec()




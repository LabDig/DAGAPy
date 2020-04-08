from PyQt5 import QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class MplCanvas(FigureCanvas):
     def __init__(self):
         # setup Matplotlib Figure and Axis
         self.fig = Figure(tight_layout=True) #fit in widget
         self.ax = self.fig.add_subplot(111)
         # initialization of the canvas
         FigureCanvas.__init__(self, self.fig)
         # we define the widget as expandable
         FigureCanvas.setSizePolicy(self,
                                    QtWidgets.QSizePolicy.Expanding,
                                    QtWidgets.QSizePolicy.Expanding)
         # notify the system of updated policy
         FigureCanvas.updateGeometry(self)
         
class mplwidget(QtWidgets.QWidget):
     def __init__(self, parent = None):
         
         # initialization of Qt MainWindow widget
         QtWidgets.QWidget.__init__(self, parent)
         # set the canvas to the Matplotlib widget
         self.canvas = MplCanvas()
         # create a vertical box layout
         self.vbl = QtWidgets.QVBoxLayout()
         # add mpl widget to vertical box
         self.vbl.addWidget(self.canvas)
         # set the layout to th vertical box
         self.setLayout(self.vbl)

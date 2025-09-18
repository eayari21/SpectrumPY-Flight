#!/opt/anaconda3/envs/idex/bin/python
# -*- coding: utf-8 -*-

title = """

//  ================================================================================
//  ||                                                                            ||
//  ||              idex_quicklook                                                ||
//  ||              ------------------------------------------------------        ||
//  ||                           I D E X  Q U I C K L O O K                       ||
//  ||              ------------------------------------------------------        ||
//  ||                                                                            ||
//  ||                __author__      = Ethan Ayari                               ||
//  ||                IMPACT/LASP, CU Boulder                                     ||
//  ||                                                                            ||
//  ||                For: IDEX INTEGRATION AND TEST                              ||
//  ||                                                                            ||
//  ||                2022                                                        ||
//  ||                                                                            ||
//  ||                                                                            ||
//  ||                Works with Python 3.10.4                                    ||
//  ||                                                                            ||
//  ================================================================================


Import an HDF5 file and vizualize its content

"""
print(title)

# %%DEPENDENCIES
import datetime
import time
import shutil
import pandas as pd
import h5py
import sys
import os
import matplotlib
import time
import argparse

import numpy as np
import qtawesome as qta
from matplotlib.backends.backend_qtagg import (FigureCanvasQTAgg,
                                            NavigationToolbar2QT as
                                            NavigationToolbar)

from PyQt6.QtCore import QSize 
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QStatusBar,
    QComboBox,
    QToolBar,
    QFileDialog,
    QVBoxLayout,
    QListWidget,
    QInputDialog,
    QMessageBox,
    QWidget,
    QLabel,
    QCheckBox,
    QPushButton,
)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# %%PRETTY PLOTS
matplotlib.use('Agg')
plt.style.use("seaborn-pastel")

fn = None
eventnum = None

# %%DO OUR BEST TO SET THE BASE DIRECTORY
try:
    basePath = __file__
except NameError:
    try:
        basePath = sys.argv[0]
        os.chdir(os.path.join(sys.argv[0]))
    except FileNotFoundError:
        basePath = None


# logging.debug("Initial basepath attempt complete. Basepath = %s"%(basePath))
basePath = os.path.abspath(basePath)

if getattr(sys, 'frozen', False):
    if hasattr(sys, "_MEIPASS"):
        application_path = os.path.join(sys._MEIPASS)
else:
    application_path = os.path.dirname(basePath)
# logging.debug("All basepath attempts complete. application_path = %s"%(application_path))

# %%SPECIFY WINDOWS ENVIRONMENT
try:
    from ctypes import windll # Only accessible on Windows OS
    myappid = "IMPACT.SpectrumPy.0.1"
    windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

except ImportError:
    pass


# %%GLOBAL EVENTNUM/AID VARIABLE
eventnum = str(1) # || Set it to event one as a start.

# %%CUSTOM TICK FORMATTING (https://stackoverflow.com/questions/25750170/show-decimal-places-and-scientific-notation-on-the-axis)
class MathTextSciFormatter(mticker.Formatter):
    def __init__(self, fmt="%1.2e"):
        self.fmt = fmt
    def __call__(self, x, pos=None):
        s = self.fmt % x
        decimal_point = '.'
        positive_sign = '+'
        tup = s.split('e')
        significand = tup[0].rstrip(decimal_point)
        sign = tup[1][0].replace(positive_sign, '')
        exponent = tup[1][1:].lstrip('0')
        if exponent:
            exponent = '10^{%s%s}' % (sign, exponent)
        if significand and exponent:
            s =  r'%s{\times}%s' % (significand, exponent)
        else:
            s =  r'%s%s' % (significand, exponent)
        return "${}$".format(s)
    
# %%SET UP INTERACTIVE PLOT
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=18.5, height=10.5, dpi=100):
        """A class to set up the interactive matplotlib plot. Uses the global
        traceNumber and channelNames variables to decide which data to display.

        Args:
           width (float):  The width of the figure.
           height (float):  The height of the figure.
           dpi (float):  Controls the quality of the interactive plot. Stands
           for "dots per inch".

        Kwargs:
           parent (FigureCanvasQTAgg): The parent canvas of the current plot.
           This is useful for opening new windows

        Returns:
           None

        Raises:
           None
           """

        global fn, eventnum
        AID = eventnum
        nrows = 6  # Make this general
        self.fig, self.ax = plt.subplots(nrows=nrows, sharex=False)
        self.fig.set_size_inches(width, height)
        self.fig.suptitle(f"Filename = {os.path.basename(fn)}, Event {eventnum}", font="Times New Roman", fontsize=30, fontweight='bold')


        # self.fig.tight_layout()

        super().__init__(self.fig)


# %%SET UP QT WINDOW OBJECT
# ||Subclass QMainWindow to customize SpectrumPY's main window
class MainWindow(QMainWindow):

    def __init__(self, filename=None, eventnumber=None, parent=None):
        """A class to set up the main pyQT user interface.

        Args:
           None

        Kwargs:
           None

        Returns:
           None

        Raises:
           None
           """
        super(MainWindow, self).__init__(parent)


        global fn, eventnum

        self.renderPlot = False
        self.setWindowTitle("SpectrumPY Waveform View")
        self.setWindowIcon(QIcon("./Plots/IMAP.ico"))  # Set window icon


        # || Get HDF5 file
        # Prompt user to select file if not provided
        if filename is None:
            filename = QFileDialog.getOpenFileName(
                self,
                "Open File",
                os.path.join(os.getcwd(), "HDF5"),
                "HDF5 Files (*.h5)",
            )[0]
        
        hdf_file = filename
        fn = filename
        self.dataset = read_hdf(hdf_file)  # || Store it in the canonical h5py format

        # Set the initial event number
        if eventnumber is None:
            print("No eventnumber chosen.")
            eventnumber = 1
            eventnum = 1
        self.eventnumber = str(eventnumber)
        eventnum = self.eventnumber
        print(f"eventnum = {eventnum}")



        # dlg = QMessageBox(self)
        # dlg.setWindowTitle("Save channel settings")
        # dlg.setText(("Would you like to save these channel names for"
        #             " future use?"))
        # dlg.setStandardButtons(QMessageBox.StandardButton.Yes |
        #                         QMessageBox.StandardButton.No)
        # dlg.setIcon(QMessageBox.Icon.Question)
        # button = dlg.exec()


        self.v_layout = QVBoxLayout()  # || Set the layout of the main window

        self._createMenuBar()
        self._createActions()

        self.sc = MplCanvas(self, width=18.5, height=10.5, dpi=100)  # || Add in the figure object
        self.toolbar = NavigationToolbar(self.sc, self)  # || Nav toolbar
        self.setCentralWidget(self.sc)  # || Make the figure the main object of the window.


        # || Build the main plotting window
        plt.xlabel(r"Time [$\mu$s]", font="Times New Roman", fontsize=20, fontweight='bold', labelpad=4)
        # for eventnum in self.dataset.key
        # s():  # || Iterate through each event
        # Format with 2 decimal places
        # plt.gca().yaxis.set_major_formatter(MathTextSciFormatter("%1.2e"))

        axdex = 0

        for count, label in enumerate(self.dataset[eventnum].keys()):  # || Iterate through each channel
           label = "Some Metadata Analysis"
           
        #    data = self.dataset[f"/{eventnum}/{label}/Ion Grid Rise Time (s)"]
        #    data = self.dataset[f"/{eventnum}/{label}/Ion Grid Amplitude (pC)"]
        #    data = self.dataset[f"/{eventnum}/{label}/Ion Grid Trigger (s)"]
           
        #   "Target Rise Time (s)", "Target amplitude (pC)", "Target Trigger (s)"
           
           if "Time" not in label and "Meta" not in label and "Analysis" not in label:

                # print(count)
                # print(label)
                # plt.gca().yaxis.set_major_formatter(MathTextSciFormatter("%1.2e"))
                axdex += 1
                self.sc.ax[axdex].tick_params(axis='both', which='major', labelsize=12)
                self.sc.ax[axdex].tick_params(axis='both', which='minor', labelsize=12)
                # self.sc.ax[axdex].yaxis.set_major_formatter(mticker.LogFormatterSciNotation(base=10,minor_thresholds=(10,10)))
                self.sc.ax[axdex].yaxis.set_major_formatter(MathTextSciFormatter("%1.2e"))  # || Custom axes class


                # self.sc.ax[axdex].ticklabel_format(style='sci', axis='both', scilimits=(0,0), useMathText=True)
                # self.sc.ax[axdex].set_xlabel("Time [s]", fontsize=16, labelpad=40)
                # self.sc.ax[axdex].set_ylabel(label, font="Times New Roman", fontweight="bold", fontsize=15, labelpad=70)
                self.sc.ax[axdex].grid(True)
                data = self.dataset[f"/{eventnum}/{label}"]
                # print(f"Length of data {axdex} = {len(data)}")

                # || For events with one time array
                # self.sc.ax[axdex].plot(self.dataset[f"{eventnum}/Time"][:], data[:], c='b', lw=1.0)
                # || Time is determined based on channel name
                # if("TOF" in label):  # || Time-of-flight channels sample at 130 MHz
                #     self.sc.ax[axdex].plot(self.dataset[f"{eventnum}/Time (high sampling)"][:], data[:], c='b', lw=1.0)
                # else:  # || Analog channels sample at 4 MHz
                #     self.sc.ax[axdex].plot(self.dataset[f"{eventnum}/Time (low sampling)"][:], data[:], c='b', lw=1.0)
        self.sc.ax[0].plot(self.dataset[f"{eventnum}/Time (high sampling)"][:], self.dataset[f"{eventnum}/TOF L"][:])
        self.sc.ax[0].set_ylabel("TOF L", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)
        self.sc.ax[1].plot(self.dataset[f"{eventnum}/Time (high sampling)"][:], self.dataset[f"{eventnum}/TOF M"][:])
        self.sc.ax[1].set_ylabel("TOF M", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)
        self.sc.ax[2].plot(self.dataset[f"{eventnum}/Time (high sampling)"][:], self.dataset[f"{eventnum}/TOF H"][:])
        self.sc.ax[2].set_ylabel("TOF H", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)

        self.sc.ax[3].plot(self.dataset[f"{eventnum}/Time (low sampling)"][:], self.dataset[f"{eventnum}/Ion Grid"][:])
        # self.sc.ax[3].plot(self.dataset[f"{eventnum}/Analysis/Ion GridFitTime"][:], self.dataset[f"{eventnum}/Analysis/Ion GridFitResult"][:])
        self.sc.ax[3].set_ylabel("Ion Grid", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)

        self.sc.ax[4].plot(self.dataset[f"{eventnum}/Time (low sampling)"][:], self.dataset[f"{eventnum}/Target L"][:])
        # self.sc.ax[4].plot(self.dataset[f"{eventnum}/Analysis/Target LFitTime"][:], self.dataset[f"{eventnum}/Analysis/Target LFitResult"][:])
        self.sc.ax[4].set_ylabel("Target L", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)

        self.sc.ax[5].plot(self.dataset[f"{eventnum}/Time (low sampling)"][:], self.dataset[f"{eventnum}/Target H"][:])
        # self.sc.ax[5].plot(self.dataset[f"{eventnum}/Analysis/Target HFitTime"][:], self.dataset[f"{eventnum}/Analysis/Target HFitResult"][:])
        self.sc.ax[5].set_ylabel("Target H", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)




        # || Create shot selection tool
        shots = np.linspace(1, len(self.dataset.keys()), len(self.dataset.keys()))
        self.shot_list = QComboBox()  # || Selector tool for shots
        self.shot_list.addItems([str(int(shot)) for shot in shots])
        # Set the default selection to the event number
        if eventnumber is not None:
            index = eventnumber - 1  # Convert 1-based to 0-based index
            if 0 <= index < len(shots):  # Ensure index is within bounds
                self.shot_list.setCurrentIndex(index)
        self.shot_list.currentIndexChanged.connect(self.event_changed)
        self.shot_list.setFixedWidth(600)



        self.v_layout.addStretch()
        self.v_layout.addWidget(self.toolbar)
        self.v_layout.addWidget(self.shot_list)

        # self.v_layout.addStretch()
        self.sc.setLayout(self.v_layout)
        # self.v_layout.addStretch()


        self.show()

    # %%HANDLER FUNCTION TO PROMPT USER TO SAVE AND EXIT SAFELY
    def closeEvent(self, event):
 
        dlg = QMessageBox(self)
        dlg.setWindowTitle("Close window")
        dlg.setText(("Are you sure you want to quit SpectrumPY?"))
        dlg.setStandardButtons(QMessageBox.StandardButton.Yes |
                               QMessageBox.StandardButton.No)
        dlg.setIcon(QMessageBox.Icon.Question)
        button = dlg.exec()

        if button == QMessageBox.StandardButton.Yes:  # || Update h5 file
            pass

        else:
            event.ignore()
    
    # %%SELECTION FUNCTION FOR EVENTNUM TOOL
    def event_changed(self, index):
        print(f"Changing to event {index+1}")
        global eventnum
        eventnum = str(index+1)
        self.sc.fig.suptitle(f"Filename = {os.path.basename(fn)}, Event {eventnum}", font="Times New Roman", fontsize=30, fontweight='bold')

        for ax in self.sc.ax:  # || Get rid of existing plots
            ax.clear()

        # || Build the main plotting window
        plt.xlabel(r"Time [$\mu$s]", font="Times New Roman", fontsize=20, fontweight='bold', labelpad=4)
        # for eventnum in self.dataset.keys():  # || Iterate through each event
        axdex = 0

        for count, label in enumerate(self.dataset[eventnum].keys()):  # || Iterate through each channel
            if any(keyword in label for keyword in ["Time", "Meta", "Analysis", "Mass", "SpiceData"]):
                
                self.sc.ax[axdex].tick_params(axis='both', which='major', labelsize=12)
                self.sc.ax[axdex].tick_params(axis='both', which='minor', labelsize=12)
                self.sc.ax[axdex].yaxis.set_major_formatter(MathTextSciFormatter("%1.2e"))  # || Custom axes class

                # self.sc.ax[axdex].ticklabel_format(style='sci', axis='both', scilimits=(0,0), useMathText=True)
                # self.sc.ax[axdex].set_xlabel("Time [s]", fontsize=16, labelpad=40)
                self.sc.ax[axdex].set_ylabel(label, font="Times New Roman", fontweight="bold", fontsize=12, labelpad=70)
                self.sc.ax[axdex].grid(True)
                data = self.dataset[f"/{eventnum}/{label}"]
                axdex += 1
                print(label, axdex)
                # print(f"Length of data {axdex} = {len(data)}")
                # self.sc.ax[axdex].plot(self.dataset[f"{eventnum}/Time"][:], data[:], c='b', lw=1.0)

                # || Time is determined based on channel name
                # if("TOF" in label):  # || Time-of-flight channels sample at 130 MHz
                #     self.sc.ax[axdex].plot(self.dataset[f"{eventnum}/Time (high sampling)"][:], data[:], c='b', lw=1.0)
                # else:  # || Analog channels sample at 4 MHz
                #     self.sc.ax[axdex].plot(self.dataset[f"{eventnum}/Time (low sampling)"][:], data[:], c='b', lw=1.0)
        # plt.tight_layout()
        self.sc.ax[0].plot(self.dataset[f"{eventnum}/Time (high sampling)"][:], self.dataset[f"{eventnum}/TOF L"][:])
        self.sc.ax[0].set_ylabel("TOF L", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)

        self.sc.ax[1].plot(self.dataset[f"{eventnum}/Time (high sampling)"][:], self.dataset[f"{eventnum}/TOF M"][:])
        self.sc.ax[1].set_ylabel("TOF M", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)

        self.sc.ax[2].plot(self.dataset[f"{eventnum}/Time (high sampling)"][:], self.dataset[f"{eventnum}/TOF H"][:])
        self.sc.ax[2].set_ylabel("TOF H", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)

        self.sc.ax[3].plot(self.dataset[f"{eventnum}/Time (low sampling)"][:], self.dataset[f"{eventnum}/Ion Grid"][:])
        # self.sc.ax[3].plot(self.dataset[f"{eventnum}/Analysis/Ion GridFitTime"][:], self.dataset[f"{eventnum}/Analysis/Ion GridFitResult"][:])
        self.sc.ax[3].set_ylabel("Ion Grid", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)

        self.sc.ax[4].plot(self.dataset[f"{eventnum}/Time (low sampling)"][:], self.dataset[f"{eventnum}/Target L"][:])
        # self.sc.ax[4].plot(self.dataset[f"{eventnum}/Analysis/Target LFitTime"][:], self.dataset[f"{eventnum}/Analysis/Target LFitResult"][:])
        self.sc.ax[4].set_ylabel("Target L", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)

        self.sc.ax[5].plot(self.dataset[f"{eventnum}/Time (low sampling)"][:], self.dataset[f"{eventnum}/Target H"][:])
        # self.sc.ax[5].plot(self.dataset[f"{eventnum}/Analysis/Target HFitTime"][:], self.dataset[f"{eventnum}/Analysis/Target HFitResult"][:])
        self.sc.ax[5].set_ylabel("Target H", font="Times New Roman", fontsize=12, fontweight='bold', labelpad=70)
        plt.draw()

# %%CREATE MENU BAR AND FILE DROP DOWN OPTIONS
    def _createMenuBar(self):
        """This function sets up the toolbar options located above the
        interactive plot. For now, these options echo the menu options (file,
        edit, view, etc.)

        Args:
           None

        Kwargs:
           None

        Returns:
           None

        Raises:
           None
           """
        toolbar = QToolBar("My main toolbar")
        toolbar.setIconSize(QSize(16, 16))
        self.addToolBar(toolbar)
        scopecsvIcon = qta.icon("mdi6.application-export")
        scopeImport = qta.icon("mdi.file-import")
        traceupIcon = qta.icon("ei.arrow-up")
        tracedownIcon = qta.icon("ei.arrow-down")

        trace_up = QAction(traceupIcon, "&View Next Trace",
                                self)
        trace_up.setStatusTip("Import SpectrumPY Data (*.h5)")
        trace_up.triggered.connect(self.upTrace)
        # trace_up.setCheckable(True)
        toolbar.addAction(trace_up)
        toolbar.addSeparator()

        trace_down = QAction(tracedownIcon, "&View Previous Trace",
                                self)
        trace_down.setStatusTip("Import Scope Data (*.trc)")
        trace_down.triggered.connect(self.downTrace)
        # trace_down.setCheckable(True)
        toolbar.addAction(trace_down)
        toolbar.addSeparator()

        button_action = QAction(scopeImport, "&Import SpectrumPY Data (*.h5)",
                                self)
        trace_up.setStatusTip("Import SpectrumPY Data (*.h5)")
        button_action.triggered.connect(self.importScopeData)
        # button_action.setCheckable(True)
        toolbar.addAction(button_action)

        toolbar.addSeparator()


        button_action2 = QAction(scopecsvIcon, "&Export Scope Data to .csv", self)
        button_action2.setStatusTip("Export Dataset (*.csv)")
        button_action2.triggered.connect(self.exportScopeData)
        # button_action2.setCheckable(True)
        toolbar.addAction(button_action2)

        hIcon = qta.icon("fa.h-square")
        HDF_Action = QAction(hIcon, "&Export Scope Data to HDF5", self)
        HDF_Action.setStatusTip("Export Dataset (*.h5)")
        HDF_Action.triggered.connect(self.exportToHDF5)
        toolbar.addAction(HDF_Action)
        

        sqlwriteIcon = qta.icon("mdi6.database-export")
        sql_write_action = QAction(sqlwriteIcon, "&Export SQL Data to .csv", self)
        sql_write_action.triggered.connect(self.writeSQL)

        traceIcon = qta.icon("fa5s.list-ol")
        listTraces = QAction(traceIcon, "&Choose Trace", self)
        listTraces.triggered.connect(self.chooseTrace)

        helpIcon = qta.icon("mdi.help-circle-outline")
        getHelp = QAction(helpIcon, "&Read the Docs", self)
        getHelp.triggered.connect(self.helpPage)

        channelChangeIcon = qta.icon("fa5s.digital-tachograph")
        changeChannels = QAction(channelChangeIcon, "&Change Channels", self)
        changeChannels.triggered.connect(self.changeChannel)

        QDIcon = qta.icon("mdi.speedometer")
        fitQD = QAction(QDIcon, "&Fit QD Waveform", self)
        fitQD.triggered.connect(self.fitQD)

        IonIcon = qta.icon("mdi.square-wave")
        fitIon = QAction(IonIcon, "&Fit Ion Grid Signal", self)
        fitIon.triggered.connect(self.fitIon)

        TargetIcon = qta.icon("ph.target")
        fitTarget = QAction(TargetIcon, "&Fit Target Grid Signal", self)
        fitTarget.triggered.connect(self.fitTarget)

        SQLIcon = qta.icon("fa.database")
        viewSQL = QAction(SQLIcon, "&View SQL Data", self)
        viewSQL.triggered.connect(self.viewSQLWindow)

        toolbar.addAction(listTraces)
        toolbar.addSeparator()

        toolbar.addAction(changeChannels)
        toolbar.addSeparator()

        toolbar.addAction(viewSQL)
        toolbar.addSeparator()

        toolbar.addAction(sql_write_action)
        toolbar.addSeparator()

        toolbar.addAction(fitIon)
        toolbar.addSeparator()

        toolbar.addAction(fitTarget)
        toolbar.addSeparator()

        toolbar.addAction(fitQD)
        toolbar.addSeparator()

        toolbar.addAction(getHelp)
        toolbar.addSeparator()


        self.setStatusBar(QStatusBar(self))

        menu = self.menuBar()

        file_menu = menu.addMenu("&File")
        # file_menu.addAction(button_action)
        import_submenu = file_menu.addMenu("Import Data")
        import_submenu.addAction(button_action)
        file_menu.addSeparator()
        # file_menu.addAction(button_action2)
        export_submenu = file_menu.addMenu("Export Data")
        export_submenu.addAction(button_action2)
        export_submenu.addAction(sql_write_action)
        export_submenu.addAction(HDF_Action)

        view_menu = menu.addMenu("&View")
        view_menu.addAction(listTraces)
        view_menu.addAction(changeChannels)
        view_menu.addAction(trace_up)
        view_menu.addAction(trace_down)

        run_menu = menu.addMenu("&Run")
        run_menu.addAction(fitQD)
        run_menu.addAction(fitIon)
        run_menu.addAction(fitTarget)
        run_menu.addAction(viewSQL)

        help_menu = menu.addMenu("&Help")
        help_menu.addAction(getHelp)

# %%CREATE TOOLBAR ACTIONS
    def _createActions(self):
        """Initializes the actions associated with the toolbar options.

        Args:
           None

        Kwargs:
           None

        Returns:
           None

        Raises:
           None
           """
        # Creating action using the first constructor
        self.newAction = QAction(self)
        self.newAction.setText("&New")
        # Creating actions using the second constructor
        self.openAction = QAction("&Open...", self)
        self.saveAction = QAction("&Save", self)
        self.exitAction = QAction("&Exit", self)
        self.copyAction = QAction("&Copy", self)
        self.pasteAction = QAction("&Paste", self)
        self.cutAction = QAction("C&ut", self)
        self.helpContentAction = QAction("&Help Content", self)
        self.aboutAction = QAction("&About", self)

# %%WRITE TIMES AND AMPS TO CSV
    def writeSQL(self, s):
        pass



# %%WRITE TIMES AND AMPS TO CSV
    def exportScopeData(self, s):
        pass


# %%WRITE TIMES AND AMPS TO CSV
    def exportToHDF5(self, s):
        pass
            




# %%CHANGE SHOT BEING VIEWED
    def upTrace(self, s):
        pass

# %%CHANGE SHOT BEING VIEWED
    def downTrace(self, s):
        pass

# %%CHANGE SHOT BEING VIEWED
    def chooseTrace(self, s):
        pass

# %%CHANGE SHOT BEING VIEWED
    def changeChannel(self, s):
        pass

# %%BRING UP HTML DOCS
    def helpPage(self, s):
        """A helper function that gets help from the QT main window.

        Args:
           None

        Kwargs:
           s (str): The choice of trace provided by the user.

        Returns:
           None

        Raises:
           None
           """
        pass
        # import webbrowser
        # filename = dname + '../_build/index.html'
        # webbrowser.open('file://' + os.path.realpath(filename))
        # # webbrowser.open_new_tab("../_build/index.html")

# %%FIT EXISTING TARGET WAVEFORMS
    def fitTarget(self, s):
        pass


# %%FIT EXISTING ION GRID WAVEFORMS
    def fitIon(self, s):
        pass

        
        


# %%FIT EXISTING QD WAVEFORMS
    def fitQD(self, s):
        pass

# %%BRING UP MATCHED SQL DATA
    def viewSQLWindow(self, s):
        pass


# %%UPDATE PLOT WITH CHOICE OF TRACE FILE
    def updatePlot(self, s):
        pass
        

# %%OPEN A FILE DIALOG TO IMPORT SCOPE DATA
    def importScopeData(self, s):
        pass


# %%SET UP QT POPUP WINDOW OBJECT
class ChannelChoosingWindow(QWidget):
    """
    This "window" is a child. If it has no parent, it
    will appear as a free-floating window as we want.
    """

    def __init__(self, MainWindow):
        super(ChannelChoosingWindow, self).__init__()
        # channelNames = 
        
        # super(ChannelChoosingWindow, self).__init__(parent)
        self.complete = False
        self.parent = MainWindow
        self.setWindowTitle(("SpectrumPY FM Channel Choosing Window"))

        # Set the label text for the user
        lb = QLabel("Please choose the channels you would like to display:",
                    self)
        lb.setGeometry(20, 20, 350, 20)
        lb.move(20, 20)

        self.label = QLabel('Nothing Selected')
        self.label.move(20, 55)

        self.submitButton = QPushButton("Plot Chosen Data")
        self.submitButton.clicked.connect(self.Submit_Plot)

        # Set the vertical Qt Layout
        vbox = QVBoxLayout()
        vbox.addWidget(lb)
        vbox.addWidget(self.label)
        # Create as many checkboxes as channels
        ySpacing = 100
        cB = {}
        for name in channelNames:
            cB[name] = QCheckBox('{}'.format(name), self)
            cB[name].move(20, ySpacing)
            cB[name].toggled.connect(self.Selected_Value)
            ySpacing += 20
            vbox.addWidget(cB[name])

        vbox.addWidget(self.submitButton)
        self.submitButton.setGeometry(500, 700, 60, 20)
        self.setLayout(vbox)
        self.setGeometry(60, 60, 700, 700)
        self.lblText = ''

        self.show()

    # %%DEFINE FUNCTION TO READ THE USER'S INPUT
    def Selected_Value(self):
        channelNames = None
        displayDex = None
        numDisplay = None

        # print("Selected: ", self.sender().isChecked(),
        #      "  Name: ", self.sender().text())

        if(self.sender().isChecked()):
            self.lblText += self.sender().text()
            numDisplay += 1
            displayDex.append(channelNames.index(self.sender().text()))
        else:
            self.lblText.replace(self.sender().text(), '')
            numDisplay -= 1
            displayDex.remove(channelNames.index(self.sender().text()))

        if self.lblText :
            self.label.setText('You have selected \n' + self.lblText)

        print(displayDex)

# %%DEFINE FUNCTION TO READ THE USER'S INPUT
    def Submit_Plot(self):
        pass

    # %%HANDLER FUNCTION TO TRIGGER MAIN WINDOW UPON CLOSE
    def closeEvent(self, event):
        pass



# %%TRACE FILE DISPLAY
def read_hdf(filename):
    print(filename)
    hdf5_file = h5py.File(filename, 'r')
    # print(list(hdf5_file.keys()))
    # dataset = hdf5_file['your_dataset']
    return hdf5_file

# %%TRACE FILE DISPLAY
def displayTRC(times, amps, sc):
    pass

# %%EXECUTABLE CODE BELOW
if __name__ == "__main__":
    # Set commandline arguments
    parser = argparse.ArgumentParser(description="Run the SpectrumPY application.")
    parser.add_argument(
        "--filename", nargs="?", default=None, help="Path to the HDF5 file."
    )
    parser.add_argument(
        "--eventnumber", nargs="?", type=int, default=None, help="Event number to display."
    )
    args = parser.parse_args()

    # Print each argument and its value
    for arg, value in vars(args).items():
        print(f"{arg}: {value}")


    # Start application
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    eventnum = args.eventnumber
    fn = args.filename

    w = MainWindow(filename=fn, eventnumber=eventnum)
    w.show()
    app.exec()

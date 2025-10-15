import sys
import os
import subprocess
import h5py
import importlib.util
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QComboBox, QPushButton, QSplitter, QLabel, QMessageBox
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas

# Define the path to the file and the module name
file_path = './IDEX-quicklook.py'  # Adjust the path if necessary
module_name = 'IDEX_quicklook'  # Module name without hyphen

# Load the module dynamically
spec = importlib.util.spec_from_file_location(module_name, file_path)
module = importlib.util.module_from_spec(spec)
sys.modules[module_name] = module
spec.loader.exec_module(module)

# Now you can use the MainWindow class from the module
from IDEX_quicklook import MainWindow

normal_size = 2.5
last_hovered_point = None  # Track the last hovered point
app = None

class HDF5PlotterApp(QWidget):
    def __init__(self, hdf5_folder):
        super().__init__()
        self.hdf5_folder = hdf5_folder

        # Load datasets from HDF5 file
        self.datasets, self.dataset_sources = self.load_datasets()

        # Initialize UI
        self.init_ui()

    def load_datasets(self):
        datasets = {}
        dataset_sources = {}  # To track which file each point comes from

        # Iterate over all .h5 files in the specified folder
        for filename in os.listdir(self.hdf5_folder):
            if not filename.endswith(".h5"):
                continue  # Skip non-HDF5 files

            file_path = os.path.join(self.hdf5_folder, filename)
            print(f"Reading filepath {file_path}")
            with h5py.File(file_path, "r") as file:
                for group in file:
                    analysis_group = f"{group}/Analysis"
                    if analysis_group not in file:
                        continue
                    for name, dataset in file[analysis_group].items():
                        if "FitParams" in name or name == "TOF H":
                            continue  # Exclude unwanted datasets
                        dataset_key = name
                        if dataset_key not in datasets:
                            datasets[dataset_key] = []
                            dataset_sources[dataset_key] = []

                        # Check if the dataset is scalar or an array
                        if dataset.shape == ():  # Scalar dataset
                            datasets[dataset_key].append(dataset[()])
                            dataset_sources[dataset_key].append(filename)
                        else:  # Array dataset
                            datasets[dataset_key].extend(dataset[:])
                            dataset_sources[dataset_key].extend([filename] * len(dataset))

        # Convert lists to numpy arrays for easier manipulation
        for key in datasets:
            datasets[key] = np.array(datasets[key])
            dataset_sources[key] = np.array(dataset_sources[key])
        print(len(datasets))
        return datasets, dataset_sources


    def init_ui(self):
        layout = QVBoxLayout()

        # Create the splitter to divide the layout into control and plot sections
        splitter = QSplitter()

        # Create the controls panel
        control_panel = QWidget()
        control_layout = QVBoxLayout()

        # Dropdowns for selecting x, y, and color
        self.x_dropdown1 = QComboBox()
        self.y_dropdown1 = QComboBox()
        self.color_dropdown1 = QComboBox()

        self.x_dropdown2 = QComboBox()
        self.y_dropdown2 = QComboBox()
        self.color_dropdown2 = QComboBox()

        dataset_keys = list(self.datasets.keys())
        self.x_dropdown1.addItems(dataset_keys)
        self.y_dropdown1.addItems(dataset_keys)
        self.color_dropdown1.addItems(dataset_keys)
        self.x_dropdown2.addItems(dataset_keys)
        self.y_dropdown2.addItems(dataset_keys)
        self.color_dropdown2.addItems(dataset_keys)

        control_layout.addWidget(QLabel("Select X-axis 1:"))
        control_layout.addWidget(self.x_dropdown1)
        control_layout.addWidget(QLabel("Select Y-axis 1:"))
        control_layout.addWidget(self.y_dropdown1)
        control_layout.addWidget(QLabel("Select Marker Color 1:"))
        control_layout.addWidget(self.color_dropdown1)

        control_layout.addWidget(QLabel("Select X-axis 2:"))
        control_layout.addWidget(self.x_dropdown2)
        control_layout.addWidget(QLabel("Select Y-axis 2:"))
        control_layout.addWidget(self.y_dropdown2)
        control_layout.addWidget(QLabel("Select Marker Color 2:"))
        control_layout.addWidget(self.color_dropdown2)

        # Plot button
        self.plot_button = QPushButton("Plot")
        self.plot_button.clicked.connect(self.plot_data)
        control_layout.addWidget(self.plot_button)
        control_panel.setLayout(control_layout)


        # Create the plot panel
        plot_panel = QWidget()
        plot_layout = QVBoxLayout()

        # Matplotlib canvas
        self.figure = plt.figure(figsize=(8, 10), constrained_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        plot_layout.addWidget(self.canvas)

        # Add the navigation toolbar
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        plot_panel.setLayout(plot_layout)


        # Add the panels to the splitter
        splitter.addWidget(control_panel)
        splitter.addWidget(plot_panel)

        # Set the splitter as the central widget
        layout.addWidget(splitter)

        self.setLayout(layout)
        # plt.subplots_adjust(hspace=0.3)  # Increase space between subplots
        # self.figure.tight_layout()  # Ensure everything fits well within the figure


    def plot_data(self):
        # Get selected keys for both plots
        x_key1 = self.x_dropdown1.currentText()
        y_key1 = self.y_dropdown1.currentText()
        color_key1 = self.color_dropdown1.currentText()

        x_key2 = self.x_dropdown2.currentText()
        y_key2 = self.y_dropdown2.currentText()
        color_key2 = self.color_dropdown2.currentText()

        # Retrieve data
        x_data1 = self.datasets[x_key1]
        y_data1 = self.datasets[y_key1]
        color_data1 = self.datasets[color_key1]

        x_data2 = self.datasets[x_key2]
        y_data2 = self.datasets[y_key2]
        color_data2 = self.datasets[color_key2]

        # Clear previous plots
        self.figure.clear()
        self.figure.set_constrained_layout(True)

        # Create gridspec for vertically stacked joint plots
        gs = self.figure.add_gridspec(
            2, 1, height_ratios=[3, 3], hspace=1.0, wspace=0.3
        )

        # First joint plot (top)
        ax_joint1 = self.figure.add_subplot(gs[0, 0])
        # ax_marg_x1 = self.figure.add_subplot(gs[0, 0], sharex=ax_joint1)
        # ax_marg_y1 = self.figure.add_subplot(gs[1, 1], sharey=ax_joint1)

        # Second joint plot (bottom)
        ax_joint2 = self.figure.add_subplot(gs[1, 0])
        # ax_marg_x2 = self.figure.add_subplot(gs[3, 0], sharex=ax_joint2)
        # ax_marg_y2 = self.figure.add_subplot(gs[4, 1], sharey=ax_joint2)

        # Plot first joint plot
        scatter1 = ax_joint1.scatter(np.abs(x_data1), np.abs(y_data1), c=np.abs(color_data1), s=normal_size, cmap="rainbow", alpha=0.7, picker=True)
        self.figure.colorbar(scatter1, ax=ax_joint1, label=color_key1)

        # sns.histplot(x=x_data1, ax=ax_marg_x1, kde=True, color="blue")
        # sns.histplot(y=y_data1, ax=ax_marg_y1, kde=True, color="blue", orientation="horizontal")

        ax_joint1.set_xlabel(x_key1)
        ax_joint1.set_ylabel(y_key1)
        ax_joint1.set_title(f"{x_key1} vs {y_key1} (colored by {color_key1})")

        # Plot second joint plot
        scatter2 = ax_joint2.scatter(np.abs(x_data2), np.abs(y_data2), c=np.abs(color_data2), s=normal_size, cmap="rainbow", alpha=0.7, picker=True)
        self.figure.colorbar(scatter2, ax=ax_joint2, label=color_key2)

        # sns.histplot(x=x_data2, ax=ax_marg_x2, kde=True, color="green")
        # sns.histplot(y=y_data2, ax=ax_marg_y2, kde=True, color="green", orientation="horizontal")

        ax_joint2.set_xlabel(x_key2)
        ax_joint2.set_ylabel(y_key2)
        ax_joint2.set_title(f"{x_key2} vs {y_key2} (colored by {color_key2})")

        # Hide overlapping tick labels
        # plt.setp(ax_marg_x1.get_xticklabels(), visible=False)
        # plt.setp(ax_marg_y1.get_yticklabels(), visible=False)
        # plt.setp(ax_marg_x2.get_xticklabels(), visible=False)
        # plt.setp(ax_marg_y2.get_yticklabels(), visible=False)

        # Apply tight layout
        # self.figure.tight_layout()

        def on_click(event):
            print("Click")
            if event.artist not in [scatter1, scatter2]:
                return

            ind = event.ind  # Get indices of clicked points
            closest_point = None
            min_distance = float('inf')  # Start with an infinitely large distance

            if event.artist == scatter1:
                for i in ind:
                    x, y = x_data1[i], y_data1[i]
                    filename = self.dataset_sources[x_key1][i]

                    # Calculate the Euclidean distance from the clicked point to the data point
                    click_x, click_y = event.mouseevent.xdata, event.mouseevent.ydata
                    distance = np.sqrt((x - click_x)**2 + (y - click_y)**2)

                    # Find the closest point (smallest distance)
                    if distance < min_distance:
                        min_distance = distance
                        closest_point = (x, y, filename, i)

                # Only print the closest point
                if closest_point:
                    x, y, filename, i = closest_point

                    # Find the group name in the HDF5 file
                    with h5py.File(self.hdf5_folder + filename, "r") as file:
                        group_names = list(file.keys())
                        num_groups = len(group_names)
                        local_index = i % num_groups
                        group_name = group_names[local_index]
                        print(f"Filename: {filename}, Event: Scatter 1, X = {x}, Y = {y}, Index in File = {group_name}")
                        # Start application
                        print("Starting application")
                        eventnum = group_name

                        # Start a new window for each event
                        subprocess.Popen(['python', 'IDEX-quicklook.py', "--filename", self.hdf5_folder+filename, "--eventnumber", str(eventnum)])



            elif event.artist == scatter2:
                for i in ind:
                    x, y = x_data2[i], y_data2[i]
                    filename = self.dataset_sources[x_key2][i]

                    # Calculate the Euclidean distance from the clicked point to the data point
                    click_x, click_y = event.mouseevent.xdata, event.mouseevent.ydata
                    distance = np.sqrt((x - click_x)**2 + (y - click_y)**2)

                    # Find the closest point (smallest distance)
                    if distance < min_distance:
                        min_distance = distance
                        closest_point = (x, y, filename, i)

                # Only print the closest point
                if closest_point:
                    x, y, filename, i = closest_point
                    # Find the group name in the HDF5 file
                    with h5py.File(self.hdf5_folder + filename, "r") as file:
                        group_names = list(file.keys())
                        num_groups = len(group_names)
                        local_index = i % num_groups
                        group_name = group_names[local_index]
                        eventnum = group_name

                        print(f"Filename: {filename}, Event: Scatter 2, X = {x}, Y = {y}, Index in File = {group_name}")
                        print("Starting application")
                        
                        # Start a new window for each event
                        subprocess.Popen(['python', 'IDEX-quicklook.py', "--filename", self.hdf5_folder+filename, "--eventnumber", str(eventnum)])
        def on_hover(event):
            # print("Hovering")
            global last_hovered_point
            # print(f"last_hovered_point = {last_hovered_point}")


            # Check if the event is inside the axes and whether it's in scatter1 or scatter2
            if event.inaxes != scatter1.axes and event.inaxes != scatter2.axes:
                return
            
            closest_point = None
            min_distance = float('inf')  # Start with an infinitely large distance

            # Initialize _sizes properly for scatter1 and scatter2
            scatter1._sizes = np.ones(len(x_data1)) * normal_size  # Ensure it matches the number of data points
            scatter2._sizes = np.ones(len(x_data2)) * normal_size  # Ensure it matches the number of data points

            # Check for hover over scatter1
            if event.inaxes == scatter1.axes:
                # print("Scatter 1 hover")
                for i in range(len(x_data1)):
                    x, y = x_data1[i], y_data1[i]
                    click_x, click_y = event.xdata, event.ydata
                    distance = np.sqrt((x - click_x)**2 + (y - click_y)**2)

                    if distance < min_distance:
                        min_distance = distance
                        closest_point = i

                # Enlarge the closest point by 20% if it's a new point
                if closest_point != last_hovered_point and closest_point is not None and min_distance < 1000:
                    # print(f"closest_hovered_point = {closest_point}")

                    if last_hovered_point is not None:  # Reset the size of the previous point
                        scatter1._sizes[last_hovered_point] = normal_size
                    scatter1._sizes[closest_point] = normal_size * 3  # Enlarge the new closest point
                    self.canvas.draw()

                    last_hovered_point = closest_point  # Update the last hovered point

            # Check for hover over scatter2
            elif event.inaxes == scatter2.axes:
                # print("Scatter 2 hover")
                for i in range(len(x_data2)):
                    x, y = x_data2[i], y_data2[i]
                    click_x, click_y = event.xdata, event.ydata
                    distance = np.sqrt((x - click_x)**2 + (y - click_y)**2)

                    if distance < min_distance:
                        min_distance = distance
                        closest_point = i

                # Enlarge the closest point by 20% if it's a new point
                if closest_point != last_hovered_point and closest_point is not None and min_distance < 1000:
                    if last_hovered_point is not None:  # Reset the size of the previous point
                        scatter2._sizes[last_hovered_point] = normal_size
                    scatter2._sizes[closest_point] = normal_size * 3  # Enlarge the new closest point
                    self.canvas.draw()

                    last_hovered_point = closest_point  # Update the last hovered point


        # Reset the point size when the mouse moves away
        def on_move(event):
            # print("Moving")
            global last_hovered_point
            # print(f"last_hovered_point = {last_hovered_point}")

            if event.inaxes not in [scatter1.axes, scatter2.axes]:
                return

            # if last_hovered_point is not None:
            #     # Reset size of the last hovered point when the mouse moves
            #     # Reset all sizes to normal
            #     if event.inaxes == scatter1.axes:
            #         print("Scatter 1 motion")
            #         scatter1._sizes = np.ones_like(scatter1._sizes) * normal_size  # Reset to original size
            #         self.canvas.draw()

            #     elif event.inaxes == scatter2.axes:
            #         print("Scatter 2 motion")
            #         scatter2._sizes = np.ones_like(scatter2._sizes) * normal_size  # Reset to original size
            #         self.canvas.draw()

            # last_hovered_point = None

        # Connect the motion events to the figure
        self.canvas.mpl_connect('motion_notify_event', on_hover)
        self.canvas.mpl_connect('motion_notify_event', on_move)

        # Connect the pick event to the handler
        self.canvas.mpl_connect('pick_event', on_click)

        # Redraw the canvas
        self.canvas.draw()




if __name__ == "__main__":
    app = QApplication(sys.argv)

   # Path to the folder containing HDF5 files
    hdf5_folder = "HDF5/"

    # Launch application
    main_window = HDF5PlotterApp(hdf5_folder)
    main_window.setWindowTitle("HDF5 Data Plotter")
    main_window.resize(800, 1000)
    main_window.show()

    app.exec()

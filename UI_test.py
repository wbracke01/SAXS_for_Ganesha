import sys
import os
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton,
    QHBoxLayout, QLabel, QFileDialog, QComboBox, QTextEdit, QMessageBox
)
from PySide6.QtCore import Qt

# Import your existing functions
from Form_fit import fit_file, ito_contrast, estimate_volfrac, log_results
import matplotlib.pyplot as plt


class SAXSGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SAXS Form Factor Fitting")

        # Variables
        self.data_path = None
        self.save_folder = None
        self.distribution = "volume"
        self.params = None

        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # File selection row
        file_layout = QHBoxLayout()
        self.file_label = QLabel("No data file selected")
        btn_select_file = QPushButton("Select Data File")
        btn_select_file.clicked.connect(self.select_file)
        file_layout.addWidget(self.file_label)
        file_layout.addWidget(btn_select_file)
        layout.addLayout(file_layout)

        # Save folder selection
        save_layout = QHBoxLayout()
        self.save_label = QLabel("No save folder selected")
        btn_select_save = QPushButton("Select Save Folder")
        btn_select_save.clicked.connect(self.select_save_folder)
        save_layout.addWidget(self.save_label)
        save_layout.addWidget(btn_select_save)
        layout.addLayout(save_layout)

        # Distribution type
        dist_layout = QHBoxLayout()
        dist_layout.addWidget(QLabel("Distribution:"))
        self.dist_combo = QComboBox()
        self.dist_combo.addItems(["volume", "number"])
        self.dist_combo.currentTextChanged.connect(self.set_distribution)
        dist_layout.addWidget(self.dist_combo)
        layout.addLayout(dist_layout)

        # Run button
        btn_run = QPushButton("Run Fit")
        btn_run.clicked.connect(self.run_fit)
        layout.addWidget(btn_run)

        # Output display
        self.output_box = QTextEdit()
        self.output_box.setReadOnly(True)
        layout.addWidget(self.output_box)

    def select_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Data File", "", "Data Files (*.csv *.npy)"
        )
        if path:
            self.data_path = path
            self.file_label.setText(os.path.basename(path))

    def select_save_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Save Folder")
        if folder:
            self.save_folder = folder
            self.save_label.setText(folder)

    def set_distribution(self, value):
        self.distribution = value

    def run_fit(self):
        if not self.data_path:
            QMessageBox.warning(self, "Error", "Please select a data file first.")
            return

        try:
            params = fit_file(self.data_path, self.distribution, self.save_folder)
            self.params = params

            rho_squared = ito_contrast(0.05, 'hexane')
            vol_frac, mass_conc = estimate_volfrac(params, rho_squared, self.save_folder)

            # Log results
            samplename = os.path.basename(self.data_path).split('.')[0]
            if self.save_folder:
                log_results(samplename, params, vol_frac, mass_conc, self.save_folder)

            # Show output
            self.output_box.append(f"File: {self.data_path}")
            self.output_box.append(f"Distribution: {self.distribution}")
            self.output_box.append(f"Params: {params}")
            self.output_box.append(f"Contrast: {rho_squared} x 10^20 / cm^4")
            self.output_box.append(f"Volume Fraction: {vol_frac*100:.3f} %")
            self.output_box.append(f"Mass Concentration: {mass_conc:.3f} mg/mL")
            self.output_box.append("-"*40)

            plt.show()

        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SAXSGUI()
    window.resize(600, 400)
    window.show()
    sys.exit(app.exec())

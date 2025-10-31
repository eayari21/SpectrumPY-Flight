"""Interactive GUI for the dust velocity and mass estimator."""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Callable, Optional, Sequence

import numpy as np

try:  # pragma: no cover - optional at runtime
    import h5py
except Exception:  # pragma: no cover - graceful fallback when unavailable
    h5py = None  # type: ignore[assignment]

try:  # pragma: no cover - prefer PySide6 when available
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import (
        QApplication,
        QComboBox,
        QDoubleSpinBox,
        QFileDialog,
        QGridLayout,
        QGroupBox,
        QLabel,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QScrollArea,
        QVBoxLayout,
        QWidget,
    )
except Exception:  # pragma: no cover - fallback to PyQt6
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QApplication,
        QComboBox,
        QDoubleSpinBox,
        QFileDialog,
        QGridLayout,
        QGroupBox,
        QLabel,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QScrollArea,
        QVBoxLayout,
        QWidget,
    )

from lookup.dust_estimator import (
    CoefficientTable,
    CurveParameters,
    compute_mass_from_charge,
    compute_velocity_from_collection_efficiency,
    compute_velocity_from_rise_time,
    load_default_tables,
    select_velocity,
)


def _render_formula(text: str) -> str:
    return f"<span style='font-size:16px;'>{text}</span>"


def _render_info(text: str) -> str:
    return f"<div style='font-size:13px; color:#343a40;'>{text}</div>"


VELOCITY_FORMULA_TEMPLATE = (
    "<b>v({variable})</b> = <b>A</b>·{variable}<sup><b>a₁</b></sup> · "
    "[1 + ({variable}/<b>v_b</b>)<sup><b>k</b></sup>]<sup>(<b>a₂</b> − <b>a₁</b>)/<b>k</b></sup> · "
    "[1 + ({variable}/<b>v_c</b>)<sup><b>k</b></sup>]<sup>(<b>a₃</b> − <b>a₂</b>)/<b>k</b></sup>"
)

YIELD_FORMULA_HTML = _render_formula(
    "<b>Y(v)</b> = <b>A</b>·<b>v</b><sup><b>a₁</b></sup> · "
    "[1 + (<b>v</b>/<b>v_b</b>)<sup><b>k</b></sup>]<sup>(<b>a₂</b> − <b>a₁</b>)/<b>k</b></sup> · "
    "[1 + (<b>v</b>/<b>v_c</b>)<sup><b>k</b></sup>]<sup>(<b>a₃</b> − <b>a₂</b>)/<b>k</b></sup>"
)

MASS_FORMULA_HTML = _render_formula("<b>m</b> = <b>Q</b> / <b>Y</b>(<b>v</b>)")

PARAMETER_DEFINITIONS_HTML = _render_info(
    """
    <ul style='margin:6px 0 0 18px;'>
      <li><b>A</b> – scaling coefficient for the curve.</li>
      <li><b>a₁</b>, <b>a₂</b>, <b>a₃</b> – power-law exponents for the three regimes.</li>
      <li><b>v_b</b>, <b>v_c</b> – breakpoint velocities (km/s) controlling the transitions.</li>
      <li><b>k</b> – transition sharpness controlling how quickly the curve bends.</li>
    </ul>
    """
)


class ParameterDisplay(QWidget):
    """Widget showing the numerical values of a coefficient set."""

    PARAMETER_KEYS = ("A", "a1", "a2", "a3", "vb", "vc", "k")

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        layout = QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setHorizontalSpacing(18)
        layout.setVerticalSpacing(6)
        self._value_labels: dict[str, QLabel] = {}
        for row, key in enumerate(self.PARAMETER_KEYS):
            label = QLabel(f"{key}:", self)
            label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            value = QLabel("—", self)
            value.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            value.setStyleSheet("font-weight: 500;")
            layout.addWidget(label, row, 0)
            layout.addWidget(value, row, 1)
            self._value_labels[key] = value

    def set_parameters(self, params: Optional[CurveParameters]) -> None:
        if params is None:
            for label in self._value_labels.values():
                label.setText("—")
            return
        self._value_labels["A"].setText(f"{params.A:.6g}")
        self._value_labels["a1"].setText(f"{params.a1:.6g}")
        self._value_labels["a2"].setText(f"{params.a2:.6g}")
        self._value_labels["a3"].setText(f"{params.a3:.6g}")
        self._value_labels["vb"].setText(f"{params.vb:.6g}")
        self._value_labels["vc"].setText(f"{params.vc:.6g}")
        self._value_labels["k"].setText(f"{params.k:.6g}")


class DustEstimatorWindow(QMainWindow):
    """Main application window that combines inputs, formulas, and results."""

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        *,
        h5: Optional["h5py.File"] = None,
        event_name: Optional[str] = None,
        event_names: Optional[Sequence[str]] = None,
        on_event_changed: Optional[Callable[[str], None]] = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Dust Velocity and Mass Estimator")
        self.resize(960, 760)

        if h5py is None:
            self._h5 = None
        else:
            self._h5 = h5
        self._event_names = list(event_names or [])
        if event_name is not None:
            self._event_name = event_name
            if event_name not in self._event_names:
                self._event_names.append(event_name)
        elif self._event_names:
            self._event_name = self._event_names[0]
        else:
            self._event_name = None
        self._on_event_changed = on_event_changed
        self._event_combo: Optional[QComboBox] = None
        self._file_label: Optional[QLabel] = None
        self._current_event_label: Optional[QLabel] = None

        (self._rise_table, self._ratio_table, self._yield_table) = load_default_tables()

        central = QWidget(self)
        outer_layout = QVBoxLayout(central)
        outer_layout.setContentsMargins(18, 18, 18, 18)
        outer_layout.setSpacing(14)

        title_label = QLabel("Impact Parameters")
        title_label.setStyleSheet("font-size: 22px; font-weight: 600;")
        outer_layout.addWidget(title_label)

        event_section = self._build_event_section()
        if event_section is not None:
            outer_layout.addWidget(event_section)

        description = QLabel(
            "Enter the target high-gain amplitude, ion grid amplitude, and target H rise time. "
            "The tool maps rise time and collection efficiency to a velocity using the "
            "calibration curves, then applies the yield curve to estimate the dust mass."
        )
        description.setWordWrap(True)
        description.setStyleSheet("color: #495057;")
        outer_layout.addWidget(description)

        input_group = self._build_input_group()
        outer_layout.addWidget(input_group)

        result_group = self._build_result_group()
        outer_layout.addWidget(result_group)

        formulas_widget = self._build_formula_panel()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(formulas_widget)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        outer_layout.addWidget(scroll, stretch=1)

        self.setCentralWidget(central)

        self._update_comboboxes()
        self._update_parameter_sections()
        self._update_estimates()
        self._load_event_data()

    # ---- UI construction -------------------------------------------------
    def _build_event_section(self) -> Optional[QGroupBox]:
        if self._h5 is None and not self._event_names and self._event_name is None:
            return None

        group = QGroupBox("Event data source")
        layout = QGridLayout(group)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setHorizontalSpacing(14)
        layout.setVerticalSpacing(10)

        row = 0
        if self._h5 is not None:
            filename = getattr(self._h5, "filename", None)
            display = "—"
            if isinstance(filename, bytes):
                try:
                    filename = filename.decode("utf-8")
                except Exception:
                    filename = filename.decode("utf-8", errors="ignore")
            if isinstance(filename, str) and filename:
                display = Path(filename).name or filename
            file_label = QLabel(display)
            file_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            layout.addWidget(QLabel("File:"), row, 0)
            layout.addWidget(file_label, row, 1)
            self._file_label = file_label
            row += 1

        if self._event_names:
            combo = QComboBox(group)
            combo.setEditable(False)
            self._event_combo = combo
            for name in self._event_names:
                combo.addItem(name)
            if self._event_name and self._event_name in self._event_names:
                combo.setCurrentIndex(self._event_names.index(self._event_name))
            combo.currentIndexChanged.connect(self._on_event_combo_changed)
            layout.addWidget(QLabel("Event:"), row, 0)
            layout.addWidget(combo, row, 1)
            row += 1
        elif self._event_name is not None:
            event_label = QLabel(self._event_name)
            event_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            layout.addWidget(QLabel("Event:"), row, 0)
            layout.addWidget(event_label, row, 1)
            self._current_event_label = event_label
            row += 1

        if row == 0:
            return None

        layout.setColumnStretch(1, 1)
        return group

    def _build_input_group(self) -> QGroupBox:
        group = QGroupBox("Inputs")
        layout = QGridLayout(group)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setHorizontalSpacing(14)
        layout.setVerticalSpacing(10)

        self.target_charge_spin = QDoubleSpinBox(group)
        self.target_charge_spin.setDecimals(12)
        self.target_charge_spin.setRange(0.0, 10.0)
        self.target_charge_spin.setSingleStep(1e-12)
        self.target_charge_spin.setSuffix(" C")
        self.target_charge_spin.setValue(1e-12)
        self.target_charge_spin.valueChanged.connect(self._update_estimates)

        self.ion_charge_spin = QDoubleSpinBox(group)
        self.ion_charge_spin.setDecimals(12)
        self.ion_charge_spin.setRange(0.0, 10.0)
        self.ion_charge_spin.setSingleStep(1e-12)
        self.ion_charge_spin.setSuffix(" C")
        self.ion_charge_spin.setValue(1e-12)
        self.ion_charge_spin.valueChanged.connect(self._update_estimates)

        self.rise_time_spin = QDoubleSpinBox(group)
        self.rise_time_spin.setDecimals(6)
        self.rise_time_spin.setRange(0.0, 1.0e6)
        self.rise_time_spin.setSingleStep(0.01)
        self.rise_time_spin.setSuffix(" µs")
        self.rise_time_spin.setValue(1.0)
        self.rise_time_spin.valueChanged.connect(self._update_estimates)

        self.min_velocity_spin = QDoubleSpinBox(group)
        self.min_velocity_spin.setDecimals(3)
        self.min_velocity_spin.setRange(0.1, 1000.0)
        self.min_velocity_spin.setSingleStep(0.1)
        self.min_velocity_spin.setSuffix(" km/s")
        self.min_velocity_spin.setValue(1.0)
        self.min_velocity_spin.valueChanged.connect(self._update_estimates)

        self.max_velocity_spin = QDoubleSpinBox(group)
        self.max_velocity_spin.setDecimals(3)
        self.max_velocity_spin.setRange(1.0, 2000.0)
        self.max_velocity_spin.setSingleStep(0.5)
        self.max_velocity_spin.setSuffix(" km/s")
        self.max_velocity_spin.setValue(100.0)
        self.max_velocity_spin.valueChanged.connect(self._update_estimates)

        layout.addWidget(QLabel("Target H amplitude (charge):"), 0, 0)
        layout.addWidget(self.target_charge_spin, 0, 1)
        layout.addWidget(QLabel("Ion grid amplitude (charge):"), 1, 0)
        layout.addWidget(self.ion_charge_spin, 1, 1)
        layout.addWidget(QLabel("Target H rise time:"), 2, 0)
        layout.addWidget(self.rise_time_spin, 2, 1)
        layout.addWidget(QLabel("Velocity clamp – minimum:"), 3, 0)
        layout.addWidget(self.min_velocity_spin, 3, 1)
        layout.addWidget(QLabel("Velocity clamp – maximum:"), 4, 0)
        layout.addWidget(self.max_velocity_spin, 4, 1)

        ratio_label = QLabel("Collection efficiency (η = Ion/Target):")
        ratio_label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.collection_ratio_value = QLabel("—")
        self.collection_ratio_value.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.collection_ratio_value.setStyleSheet("font-weight: 500;")
        layout.addWidget(ratio_label, 5, 0)
        layout.addWidget(self.collection_ratio_value, 5, 1)

        self.rise_combo = QComboBox(group)
        self.rise_combo.currentIndexChanged.connect(self._on_curve_changed)
        self.rise_load_button = QPushButton("Load rise-time CSV…", group)
        self.rise_load_button.clicked.connect(lambda: self._load_table("rise"))

        self.ratio_combo = QComboBox(group)
        self.ratio_combo.currentIndexChanged.connect(self._on_curve_changed)
        self.ratio_load_button = QPushButton("Load collection-efficiency CSV…", group)
        self.ratio_load_button.clicked.connect(lambda: self._load_table("ratio"))

        self.yield_combo = QComboBox(group)
        self.yield_combo.currentIndexChanged.connect(self._on_curve_changed)
        self.yield_load_button = QPushButton("Load yield CSV…", group)
        self.yield_load_button.clicked.connect(lambda: self._load_table("yield"))

        layout.addWidget(QLabel("Rise-time coefficients:"), 6, 0)
        layout.addWidget(self.rise_combo, 6, 1)
        layout.addWidget(self.rise_load_button, 6, 2)

        layout.addWidget(QLabel("Collection-efficiency coefficients:"), 7, 0)
        layout.addWidget(self.ratio_combo, 7, 1)
        layout.addWidget(self.ratio_load_button, 7, 2)

        layout.addWidget(QLabel("Yield coefficients:"), 8, 0)
        layout.addWidget(self.yield_combo, 8, 1)
        layout.addWidget(self.yield_load_button, 8, 2)

        layout.setColumnStretch(1, 1)
        return group

    def _build_result_group(self) -> QGroupBox:
        group = QGroupBox("Estimates")
        layout = QGridLayout(group)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setHorizontalSpacing(18)
        layout.setVerticalSpacing(8)

        def _add_row(row: int, label_text: str) -> QLabel:
            label = QLabel(label_text)
            label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            value = QLabel("—")
            value.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            value.setStyleSheet("font-weight: 600; font-size: 14px;")
            layout.addWidget(label, row, 0)
            layout.addWidget(value, row, 1)
            return value

        row = 0
        self.rise_velocity_value = _add_row(row, "Velocity from rise time:")
        row += 1
        self.ratio_velocity_value = _add_row(row, "Velocity from collection efficiency:")
        row += 1
        self.selected_velocity_value = _add_row(row, "Selected velocity (km/s):")
        row += 1
        self.velocity_source_value = _add_row(row, "Selection source:")
        row += 1
        self.charge_yield_value = _add_row(row, "Yield Y(v) (C/kg):")
        row += 1
        self.mass_value = _add_row(row, "Estimated mass (kg):")
        row += 1

        stored_header = QLabel("Stored values from file:")
        stored_header.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        stored_header.setStyleSheet("font-weight: 600; color: #495057;")
        layout.addWidget(stored_header, row, 0, 1, 2)
        row += 1

        self.stored_charge_value = _add_row(row, "Impact charge (pC):")
        row += 1
        self.stored_velocity_value = _add_row(row, "Velocity (km/s):")
        row += 1
        self.stored_mass_value = _add_row(row, "Dust mass (kg):")
        row += 1
        self.stored_yield_value = _add_row(row, "Charge yield (pC/kg):")
        row += 1
        self.stored_rise_velocity_value = _add_row(row, "Velocity from rise (km/s):")
        row += 1
        self.stored_ratio_velocity_value = _add_row(row, "Velocity from collection efficiency (km/s):")
        row += 1
        self.stored_velocity_source_value = _add_row(row, "Velocity source:")

        return group

    def _clear_stored_values(self) -> None:
        for label in (
            self.stored_charge_value,
            self.stored_velocity_value,
            self.stored_mass_value,
            self.stored_yield_value,
            self.stored_rise_velocity_value,
            self.stored_ratio_velocity_value,
        ):
            label.setText("—")
        self.stored_velocity_source_value.setText("—")

    def _set_spin_from_dataset(self, spin: QDoubleSpinBox, value: Optional[float]) -> None:
        spin.blockSignals(True)
        if value is None or not math.isfinite(value):
            spin.setValue(spin.minimum())
        else:
            clamped = min(max(value, spin.minimum()), spin.maximum())
            spin.setValue(clamped)
        spin.blockSignals(False)

    def _read_scalar_dataset(self, dataset_path: str) -> Optional[float]:
        if h5py is None or self._h5 is None:
            return None
        try:
            dataset = self._h5[dataset_path]
        except Exception:
            return None
        if not isinstance(dataset, h5py.Dataset):
            return None
        try:
            data = np.asarray(dataset[()], dtype=float).ravel()
        except Exception:
            return None
        if data.size == 0:
            return None
        value = float(data[0])
        if not math.isfinite(value):
            return None
        return value

    def _read_string_dataset(self, dataset_path: str) -> Optional[str]:
        if h5py is None or self._h5 is None:
            return None
        try:
            dataset = self._h5[dataset_path]
        except Exception:
            return None
        if not isinstance(dataset, h5py.Dataset):
            return None
        data = dataset[()]
        if isinstance(data, str):
            return data or None
        if isinstance(data, bytes):
            try:
                decoded = data.decode("utf-8")
            except Exception:
                decoded = data.decode("utf-8", errors="ignore")
            return decoded or None
        try:
            array = np.asarray(data).ravel()
        except Exception:
            return None
        if array.size == 0:
            return None
        element = array[0]
        if isinstance(element, str):
            return element or None
        if isinstance(element, bytes):
            try:
                decoded = element.decode("utf-8")
            except Exception:
                decoded = element.decode("utf-8", errors="ignore")
            return decoded or None
        return None

    def _load_event_data(self) -> None:
        self._clear_stored_values()

        if self._current_event_label is not None and self._event_name is not None:
            self._current_event_label.setText(self._event_name)

        if h5py is None or self._h5 is None or self._event_name is None:
            self._set_spin_from_dataset(self.target_charge_spin, None)
            self._set_spin_from_dataset(self.ion_charge_spin, None)
            self._set_spin_from_dataset(self.rise_time_spin, None)
            self._update_estimates()
            return

        try:
            self._h5[self._event_name]
        except Exception:
            self._set_spin_from_dataset(self.target_charge_spin, None)
            self._set_spin_from_dataset(self.ion_charge_spin, None)
            self._set_spin_from_dataset(self.rise_time_spin, None)
            self._update_estimates()
            return

        charge_c = self._read_scalar_dataset(f"/{self._event_name}/Analysis/Target H Impact Charge")
        ion_charge_c = self._read_scalar_dataset(f"/{self._event_name}/Analysis/Ion Grid Impact Charge")
        rise_time_us = self._read_scalar_dataset(f"/{self._event_name}/Analysis/Target H 10-90 Risetime")
        velocity = self._read_scalar_dataset(f"/{self._event_name}/Analysis/Target H Velocity Estimate")
        mass = self._read_scalar_dataset(f"/{self._event_name}/Analysis/Target H Dust Mass Estimate")
        yield_pc_per_kg = self._read_scalar_dataset(f"/{self._event_name}/Analysis/Target H Charge Yield Estimate")
        rise_velocity = self._read_scalar_dataset(f"/{self._event_name}/Analysis/Target H Velocity From Rise")
        ratio_velocity = self._read_scalar_dataset(f"/{self._event_name}/Analysis/Target H Velocity From Ratio")
        velocity_source = self._read_string_dataset(f"/{self._event_name}/Analysis/Target H Velocity Source")

        self._set_spin_from_dataset(self.target_charge_spin, charge_c)
        self._set_spin_from_dataset(self.ion_charge_spin, ion_charge_c)
        self._set_spin_from_dataset(self.rise_time_spin, rise_time_us)

        self._set_value_label(
            self.stored_charge_value,
            charge_c * 1.0e12 if charge_c is not None and math.isfinite(charge_c) else None,
            "pC",
        )
        self._set_value_label(self.stored_velocity_value, velocity, "km/s")
        self._set_value_label(self.stored_mass_value, mass, "kg")
        self._set_value_label(self.stored_yield_value, yield_pc_per_kg, "pC/kg")
        self._set_value_label(self.stored_rise_velocity_value, rise_velocity, "km/s")
        self._set_value_label(self.stored_ratio_velocity_value, ratio_velocity, "km/s")
        if velocity_source:
            self.stored_velocity_source_value.setText(velocity_source)
        else:
            self.stored_velocity_source_value.setText("—")

        self._update_estimates()

    def _on_event_combo_changed(self, index: int) -> None:  # noqa: ARG002 - required signature
        if self._event_combo is None or index < 0 or index >= self._event_combo.count():
            return
        name = self._event_combo.itemText(index)
        if not name:
            return
        if name == self._event_name:
            return
        self._event_name = name
        self._load_event_data()
        if self._on_event_changed is not None:
            try:
                self._on_event_changed(name)
            except Exception:
                pass

    def set_current_event(self, event_name: Optional[str]) -> None:
        if not event_name:
            return
        self._event_name = event_name
        if self._event_combo is not None:
            index = self._event_combo.findText(event_name)
            if index == -1:
                self._event_combo.addItem(event_name)
                self._event_names.append(event_name)
                index = self._event_combo.count() - 1
            if self._event_combo.currentIndex() != index:
                self._event_combo.blockSignals(True)
                self._event_combo.setCurrentIndex(index)
                self._event_combo.blockSignals(False)
        elif self._current_event_label is not None:
            self._current_event_label.setText(event_name)
        if event_name not in self._event_names:
            self._event_names.append(event_name)
        self._load_event_data()

    def _build_formula_panel(self) -> QWidget:
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(18)

        self.rise_formula_group, self.rise_parameters = self._build_formula_group(
            title="Velocity from Target H Rise Time",
            formula_html=_render_formula(
                VELOCITY_FORMULA_TEMPLATE.format(variable="<b>t_r</b>")
            ),
            variable_description=_render_info(
                "<b>t_r</b> is the Target H rise time in microseconds; the curve returns velocity in km/s."
            ),
        )
        layout.addWidget(self.rise_formula_group)

        self.ratio_formula_group, self.ratio_parameters = self._build_formula_group(
            title="Velocity from Collection Efficiency",
            formula_html=_render_formula(
                VELOCITY_FORMULA_TEMPLATE.format(variable="<b>η</b>")
            ),
            variable_description=_render_info(
                "<b>η</b> is the collection efficiency ratio Ion Grid / Target H (dimensionless)."
            ),
        )
        layout.addWidget(self.ratio_formula_group)

        yield_group = QGroupBox("Yield Curve and Mass Estimate")
        yield_layout = QVBoxLayout(yield_group)
        yield_layout.setContentsMargins(12, 12, 12, 12)
        yield_layout.setSpacing(10)

        yield_formula_label = QLabel()
        yield_formula_label.setTextFormat(Qt.TextFormat.RichText)
        yield_formula_label.setWordWrap(True)
        yield_formula_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        yield_formula_label.setText(YIELD_FORMULA_HTML)
        yield_layout.addWidget(yield_formula_label)

        yield_description = QLabel(
            _render_info(
                "<b>v</b> is the selected velocity (km/s). <b>Y(v)</b> returns the target charge yield in C/kg."
            )
        )
        yield_description.setTextFormat(Qt.TextFormat.RichText)
        yield_description.setWordWrap(True)
        yield_layout.addWidget(yield_description)

        self.yield_parameters = ParameterDisplay()
        yield_layout.addWidget(self.yield_parameters)

        mass_formula_label = QLabel()
        mass_formula_label.setTextFormat(Qt.TextFormat.RichText)
        mass_formula_label.setWordWrap(True)
        mass_formula_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        mass_formula_label.setText(MASS_FORMULA_HTML)
        yield_layout.addWidget(mass_formula_label)

        mass_description = QLabel(
            _render_info("<b>Q</b> is the Target H charge amplitude in coulombs.")
        )
        mass_description.setTextFormat(Qt.TextFormat.RichText)
        mass_description.setWordWrap(True)
        yield_layout.addWidget(mass_description)

        yield_layout.addWidget(self._build_parameter_definition_label())

        layout.addWidget(yield_group)
        layout.addStretch(1)
        return container

    def _build_formula_group(
        self,
        *,
        title: str,
        formula_html: str,
        variable_description: str,
    ) -> tuple[QGroupBox, ParameterDisplay]:
        group = QGroupBox(title)
        layout = QVBoxLayout(group)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)

        formula_label = QLabel()
        formula_label.setTextFormat(Qt.TextFormat.RichText)
        formula_label.setWordWrap(True)
        formula_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        formula_label.setText(formula_html)
        layout.addWidget(formula_label)

        variable_label = QLabel(variable_description)
        variable_label.setTextFormat(Qt.TextFormat.RichText)
        variable_label.setWordWrap(True)
        layout.addWidget(variable_label)

        params = ParameterDisplay()
        layout.addWidget(params)

        layout.addWidget(self._build_parameter_definition_label())
        return group, params

    def _build_parameter_definition_label(self) -> QLabel:
        label = QLabel(PARAMETER_DEFINITIONS_HTML)
        label.setTextFormat(Qt.TextFormat.RichText)
        label.setWordWrap(True)
        return label

    # ---- Table loading and selection ------------------------------------
    def _update_comboboxes(self) -> None:
        self._populate_combo(self.rise_combo, self._rise_table)
        self._populate_combo(self.yield_combo, self._yield_table)
        self._populate_combo(self.ratio_combo, self._ratio_table)
        self.ratio_combo.setEnabled(self._ratio_table is not None)
        self.ratio_load_button.setEnabled(True)

    def _populate_combo(self, combo: QComboBox, table: Optional[CoefficientTable]) -> None:
        combo.blockSignals(True)
        combo.clear()
        if table is None:
            combo.addItem("Unavailable", None)
            combo.setEnabled(False)
        else:
            for label in table.labels():
                display = label.replace("_", " ").title() if label != label.upper() else label
                combo.addItem(display, userData=label)
            combo.setEnabled(True)
            combo.setCurrentIndex(0)
        combo.blockSignals(False)

    def _load_table(self, which: str) -> None:
        dialog_caption = {
            "rise": "Select rise-time coefficient CSV",
            "ratio": "Select collection-efficiency coefficient CSV",
            "yield": "Select yield coefficient CSV",
        }.get(which, "Select coefficient CSV")
        path, _ = QFileDialog.getOpenFileName(self, dialog_caption, str(Path.cwd()), "CSV Files (*.csv)")
        if not path:
            return
        try:
            table = CoefficientTable.from_csv(path, default_label=Path(path).stem or "custom")
        except Exception as exc:  # pragma: no cover - user-driven
            QMessageBox.critical(self, "Load Error", f"Unable to load coefficients:\n{exc}")
            return
        if which == "rise":
            self._rise_table = table
        elif which == "ratio":
            self._ratio_table = table
        elif which == "yield":
            self._yield_table = table
        self._update_comboboxes()
        self._update_parameter_sections()
        self._update_estimates()

    def _on_curve_changed(self, index: int) -> None:  # noqa: ARG002
        self._update_parameter_sections()
        self._update_estimates()

    def _current_parameters(self, combo: QComboBox, table: Optional[CoefficientTable]) -> Optional[CurveParameters]:
        if table is None:
            return None
        label = combo.currentData()
        if not isinstance(label, str):
            return None
        try:
            return table.get(label)
        except KeyError:
            return None

    def _update_parameter_sections(self) -> None:
        rise_params = self._current_parameters(self.rise_combo, self._rise_table)
        self.rise_parameters.set_parameters(rise_params)

        ratio_params = self._current_parameters(self.ratio_combo, self._ratio_table)
        self.ratio_parameters.set_parameters(ratio_params)
        self.ratio_formula_group.setVisible(self._ratio_table is not None)

        yield_params = self._current_parameters(self.yield_combo, self._yield_table)
        self.yield_parameters.set_parameters(yield_params)

    # ---- Calculations ----------------------------------------------------
    def _update_estimates(self) -> None:
        charge = float(self.target_charge_spin.value())
        ion_charge = float(self.ion_charge_spin.value())
        ratio = math.nan
        if charge > 0.0 and ion_charge > 0.0:
            ratio = ion_charge / charge
            self.collection_ratio_value.setText(f"{ratio:.6g}")
        else:
            self.collection_ratio_value.setText("—")

        rise_time = float(self.rise_time_spin.value())
        if rise_time <= 0.0:
            rise_time = math.nan

        rise_params = self._current_parameters(self.rise_combo, self._rise_table)
        ratio_params = self._current_parameters(self.ratio_combo, self._ratio_table)
        yield_params = self._current_parameters(self.yield_combo, self._yield_table)

        rise_velocity = self._compute_velocity(rise_time, rise_params, "rise")
        ratio_velocity = self._compute_velocity(ratio, ratio_params, "ratio")

        self._set_value_label(self.rise_velocity_value, rise_velocity, "km/s")
        self._set_value_label(self.ratio_velocity_value, ratio_velocity, "km/s")

        selected_velocity = None
        source = "—"
        min_velocity = float(self.min_velocity_spin.value())
        max_velocity = float(self.max_velocity_spin.value())
        if min_velocity >= max_velocity:
            max_velocity = min_velocity + 1.0
            self.max_velocity_spin.blockSignals(True)
            self.max_velocity_spin.setValue(max_velocity)
            self.max_velocity_spin.blockSignals(False)
        try:
            selected_velocity, source = select_velocity(
                rise_velocity,
                ratio_velocity,
                (min_velocity, max_velocity),
            )
        except ValueError:
            selected_velocity = None
            source = "no valid velocity"
        self._set_value_label(self.selected_velocity_value, selected_velocity, "km/s")
        self.velocity_source_value.setText(source)

        if yield_params is None or selected_velocity is None or not math.isfinite(selected_velocity) or selected_velocity <= 0:
            self.charge_yield_value.setText("—")
            self.mass_value.setText("—")
            return

        try:
            mass, yield_value = compute_mass_from_charge(charge, selected_velocity, yield_params)
        except Exception:
            self.charge_yield_value.setText("—")
            self.mass_value.setText("—")
            return

        self._set_value_label(self.charge_yield_value, yield_value, "C/kg")
        self._set_value_label(self.mass_value, mass, "kg")

    def _compute_velocity(
        self,
        value: float,
        params: Optional[CurveParameters],
        mode: str,
    ) -> Optional[float]:
        if params is None or not math.isfinite(value) or value <= 0.0:
            return None
        try:
            if mode == "rise":
                return compute_velocity_from_rise_time(value, params)
            return compute_velocity_from_collection_efficiency(value, params)
        except Exception:
            return None

    def _set_value_label(self, label: QLabel, value: Optional[float], unit: str) -> None:
        if value is None or not math.isfinite(value):
            label.setText("—")
        else:
            label.setText(f"{value:.6g} {unit}")


def launch_dust_estimator_window(
    *,
    parent: Optional[QWidget] = None,
    h5: Optional["h5py.File"] = None,
    event_name: Optional[str] = None,
    event_names: Optional[Sequence[str]] = None,
    on_event_changed: Optional[Callable[[str], None]] = None,
) -> DustEstimatorWindow:
    """Create a dust estimator window for embedding in other Qt applications."""

    return DustEstimatorWindow(
        parent=parent,
        h5=h5,
        event_name=event_name,
        event_names=event_names,
        on_event_changed=on_event_changed,
    )


def main() -> int:
    app = QApplication(sys.argv)
    try:
        window = launch_dust_estimator_window()
    except Exception as exc:
        QMessageBox.critical(None, "Impact Parameters Error", f"Unable to start the impact parameters tool:\n{exc}")
        return 1
    window.show()
    return app.exec()


if __name__ == "__main__":  # pragma: no cover - GUI entry point
    sys.exit(main())

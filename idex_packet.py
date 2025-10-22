#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

"""IDEX packet tools for generating higher level data products.

This module ingests Level-0 (L0) telemetry packets stored in the
``Data/`` directory and extracts waveform and housekeeping information for
the Instrument for Dust Experiment (IDEX). The primary workflow is to
transform raw telemetry into Level-1A (L1A) through Level-2C (L2C)
products by decoding waveforms, fitting instrument responses, and writing
standardized HDF5/CDF outputs. The original implementation targeted
Python 3.8.10 and relies on the LASP packet definitions along with a
collection of numerical utilities.

The file retains extensive functionality that is used by several science
pipelines. Documentation strings were refreshed to clarify intent without
altering established behaviour.

__author__ = "Ethan Ayari & Gavin Medley, Institute for Modeling Plasmas,\
 Atmospheres and Cosmic Dust"
"""

# || Python libraries
import argparse
import os
import socket
import bitstring
import h5py
import shutil
import struct
import matplotlib.pyplot as plt
from collections import defaultdict
from pathlib import Path
from typing import Dict, Tuple

from idex_variable_definitions import load_variable_definitions
# try:
#     plt.style.use("seaborn-pastel")
# except:
#     plt.style.use("seaborn-v0_8-pastel")
plt.style.use("seaborn-pastel")
import numpy as np

try:
    import cupy as cp  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    cp = None

from datetime import datetime, timedelta, timezone

from scipy.optimize import curve_fit
from scipy.signal import detrend, butter, filtfilt, find_peaks
from scipy.special import erfc


# || LASP software
from lasp_packets import xtcedef  # Gavin Medley's xtce UML implementation
from lasp_packets import parser  # Gavin Medley's constant bitstream implementation
from rice_decode import idex_rice_Decode
from time2mass import time2mass
import cdflib.cdfwrite as cdfwrite
import cdflib.cdfread as cdfread

# %%IDEX ION GRID FUNCTION DEFINITON
def IDEXIonGrid(x, P0, P1, P4, P5, P6):
    """Evaluate the analytic ion grid response model for the target signal.

    Parameters
    ----------
    x
        Time axis in microseconds.
    P0, P1, P4, P5, P6
        Model coefficients representing the impact time, offset, amplitude,
        rise time, and decay time respectively.
    """

    return P1 + np.heaviside(x-P0, 0) * ( P4 * (1.0 - np.exp(-(x-P0)/P5)) * np.exp( -(x-P0)/P6))

# Define the EMG function
def EMG(x, amplitude, mu, sigma, lam):
    """Return the exponentially modified Gaussian evaluated at ``x``.

    Parameters
    ----------
    x
        Time coordinate.
    amplitude, mu, sigma, lam
        Exponentially modified Gaussian parameters (area, mean, standard
        deviation, and exponential decay rate).
    """

    sigma = abs(sigma) if abs(sigma) > 1.0e-15 else 1.0e-15
    lam = abs(lam) if abs(lam) > 1.0e-15 else 1.0e-15
    prefactor = (lam * amplitude) / 2
    exponent = np.exp((lam / 2) * (2 * mu + lam * sigma**2 - 2 * x))
    erfc_part = erfc((mu + lam * sigma**2 - x) / (np.sqrt(2) * sigma))
    return prefactor * exponent * erfc_part


def EMG_CDF(x, mu, sigma, lam):
    """Return the cumulative distribution function of a unit EMG."""

    sigma = abs(sigma) if abs(sigma) > 1.0e-15 else 1.0e-15
    lam = abs(lam) if abs(lam) > 1.0e-15 else 1.0e-15
    z = (x - mu) / sigma
    first = 0.5 * erfc(-z / np.sqrt(2.0))
    exponent = np.exp(-lam * (x - mu) + 0.5 * (lam * sigma) ** 2)
    second = 0.5 * erfc(-(z - lam * sigma) / np.sqrt(2.0))
    return np.clip(first - exponent * second, 0.0, 1.0)

# Function to calculate the area under the EMG fit curve
def calculate_area_under_emg(x_slice, param):
    """Integrate the exponentially modified Gaussian across ``x_slice``.

    Parameters
    ----------
    x_slice
        Domain over which to integrate the fitted EMG profile.
    param
        Optimized EMG parameters. If the fit fails the value ``0`` is
        returned and the integral is skipped.

    Returns
    -------
    float
        Estimated area under the curve for the selected slice.
    """

    if isinstance(param, (list, tuple, np.ndarray)) and len(param) >= 4:
        amplitude, mu, sigma, lam = param[:4]
        if len(x_slice) == 0:
            return 0.0
        start = float(np.min(x_slice))
        end = float(np.max(x_slice))
        if not np.isfinite(start) or not np.isfinite(end) or end <= start:
            return 0.0
        cdf_start = EMG_CDF(start, mu, sigma, lam)
        cdf_end = EMG_CDF(end, mu, sigma, lam)
        area = float(amplitude) * float(cdf_end - cdf_start)
        return float(area)
    return 0.0

# Helper utilities for interacting with the HDF5 hierarchy
def _ensure_group(hdf_obj, path: str) -> str:
    """Return ``path`` stripped of a leading slash and create groups as needed."""

    dataset_name = path.lstrip('/')
    group_path = os.path.dirname(dataset_name)
    if group_path:
        hdf_obj.require_group(group_path)
    return dataset_name


def create_dataset_if_not_exists(hdf5_obj, dataset_path, data, **kwargs):
    """Create ``dataset_path`` within ``hdf5_obj`` if it is not present."""

    dataset_name = _ensure_group(hdf5_obj, dataset_path)
    if dataset_name in hdf5_obj:
        return hdf5_obj[dataset_name]
    return hdf5_obj.create_dataset(dataset_name, data=data, **kwargs)


def create_or_replace_dataset(hdf5_obj, dataset_path, data, **kwargs):
    """Create ``dataset_path`` and replace existing data if necessary."""

    dataset_name = _ensure_group(hdf5_obj, dataset_path)

    if dataset_name in hdf5_obj:
        del hdf5_obj[dataset_name]
    return hdf5_obj.create_dataset(dataset_name, data=data, **kwargs)

# Fit routine for EMG
def FitEMG(time, amplitude):
    """Fit an exponentially modified Gaussian to the provided waveform.

    Parameters
    ----------
    time
        Sample axis for the waveform.
    amplitude
        Measured detector response in data numbers.

    Returns
    -------
    tuple
        The fitted parameters, covariance matrix, peak amplitude estimate,
        and the evaluated model at the original ``time`` points.
    """

    x = np.asarray(time)
    y = np.asarray(amplitude)

    # || Initial Guess for the parameters of the EMG
    mu_guess = x[np.argmax(y)]  # Initial guess for the mean
    sigma_guess = np.std(x) / 10  # Initial guess for standard deviation
    lam_guess = 1 / max(x[-1] - x[0], 1.0e-6)  # Initial guess for decay rate
    area_guess = np.trapz(np.clip(y, 0.0, None), x)
    if not np.isfinite(area_guess) or area_guess <= 0:
        area_guess = max(np.max(np.clip(y, 0.0, None)) * max(x[-1] - x[0], 1.0e-6), 1.0e-6)

    p0 = [area_guess, mu_guess, sigma_guess, lam_guess]  # Initial parameter guesses

    # Fit the data using curve_fit
    try:
        param, param_cov = curve_fit(
            lambda t, area, mu, sigma, lam: EMG(t, area, mu, sigma, lam),
            x,
            y,
            p0=p0,
            maxfev=100_000,
        )

        # Generate the fitted curve
        area_fit, mu_fit, sigma_fit, lam_fit = param
        area_fit = abs(float(area_fit))
        sigma_fit = abs(float(sigma_fit))
        lam_fit = abs(float(lam_fit))
        param = np.array([area_fit, mu_fit, sigma_fit, lam_fit])
        result = EMG(x, area_fit, mu_fit, sigma_fit, lam_fit)
        sig_amp = max(result) - np.mean(y)

        return param, param_cov, sig_amp, result
    except RuntimeError as e:
        print(f"Fit failed: {e}")
        return 0, 0, 0, 0

# %%Target Signal Fitting Routine %% #

# || Very noisy due to "microphonics", so we will:
# || 1) Remove a linear baseline (y = a*x + b), and 
# || 2) Remove a sinusoidal background (y = c*sin(d*x + e)

def FitTargetSignal(time, targetAmp):
    """Fit the ion grid target signal after background subtraction.

    Parameters
    ----------
    time
        Sample axis centred on the trigger.
    targetAmp
        Raw target amplitude used for event-level charge estimation.

    Returns
    -------
    tuple
        Best-fit parameters, covariance matrix, and signal amplitude.
    """

    x = np.asarray(time)
    y = np.asarray(targetAmp)

    # || Select only raw noise (where we know the signal is not)
    mask = np.logical_and(time >= -7, time <= -5)
    print(mask)

    baselineraw = y[mask]
    baselinedomain = time[mask]

    # || Remove Linear Background

    try:
        # slopeguess = (baselineraw[len(baselineraw)-1] - baselineraw[0]) / (baselinedomain[len(baselinedomain-1)] - baselinedomain[0])
        slopeguess = 0
        linparam, lin_cov = curve_fit(LinearFit, baselinedomain, baselineraw, p0=[slopeguess,0], maxfev=100_000)
        linearbase = LinearFit(time, linparam[0], linparam[1])
        # y -= linearbase
        y = detrend(y)
    except Exception:
        print(f"Linear background not found.")
        linearbase = None


    # || Remove Sine Wave Background

    try:
        baselinedelined = y[basedex]
        sineparam, sine_cov = curve_fit(SineFit, baselinedomain, baselinedelined, p0=[max(baselinedelined), 7000, 45], maxfev=100_000)
        sinebase = SineFit(time, sineparam[0], sineparam[1], sineparam[2])
        y -= sinebase
        y = butter_lowpass_filter(y, time)


    except Exception:
        print(f"Sinusoidal background not found.")
        sinebase=None


    # print(f"idx = {idx}")
    x = time[mask]
    y = y[mask]
    pre = -2.0 # Before image charge
    # print(f"*** y is of type {type(y)}")
    # Assuming x and y are numpy arrays and `pre` is the threshold
    yBaseline = np.where(x < pre, y, np.nan)  # Set elements of y to nan where x >= pre
    yImage = y[(x >= pre) & (x < 0.0)]
    ionError = np.std(yBaseline)
    ionMean = yBaseline.mean()
    yBaseline -= ionMean
    parameters = []

    ionTime = np.array([float(x) for x in x])
    ionAmp = np.array([float(y) for y in y])

    print(f"Calculating Target Fit.")
    # || Initial Guess for the parameters of the ion grid signal

    t0 = 0.0                         # P[0] time of impact
    c = 0.                           # P[1] Constant offset
    # b = np.abs(min(yImage))     # P[2] Image amplitude
    # b = .01

    # s = 4.e-6                        # P[3] Image pulse width
    A  = max(ionAmp)  # np.abs(min(ionAmp) - max(ionAmp))     # P[4] amplitude (v)
    # A = .05

    t1 = .371                         # P[5] rise  time (s)
    t2 = .371                          # P[6] discharge time (s)
    c = 0



    filtered_time = ionTime
    filtered_signal = ionAmp
    fit_result = None

    try:
        param, param_cov = curve_fit(IDEXIonGrid, ionTime, ionAmp, p0=[t0, c, A, t1, t2], maxfev=100_000)
        result = IDEXIonGrid(ionTime, param[0], param[1], param[2], param[3], param[4])
        sig_amp = max(result) - yBaseline.mean()
        fit_result = result
    except RuntimeError as exc:
        print(f"Ion grid fit failed: {exc}")
        param, param_cov = None, None
        sig_amp = 0.0

    return (param, param_cov, sig_amp, filtered_time, filtered_signal, fit_result)


# ||
# ||
# || Generator object from LASP packets
# || to read in the data
class IDEXEvent:
    """Container for decoded IDEX packets and derived metadata.

    The class ingests a binary telemetry stream, typically sourced from the
    ``Data/`` directory, and exposes convenience helpers to plot and write
    Level-0 derived products for subsequent L1A through L2C processing.
    """

    def __init__(self, filename: str):
        """Parse a binary telemetry file and populate waveform containers.

        Parameters
        ----------
        filename
            Path to the raw IDEX L0 binary packet stream produced by the
            instrument. The parser consumes the stream, extracts per-event
            metadata, and decodes science waveforms in preparation for
            higher-level processing.
        """
        # TODO: CHge location of xml definition
        idex_xtce = 'idex_combined_science_definition.xml'
        idex_definition = xtcedef.XtcePacketDefinition(xtce_document=idex_xtce)
        # assert isinstance(idex_definition, xtcedef.XtcePacketDefinition)


        idex_packet_file = filename
        print(f"Reading in data file {idex_packet_file}")
        idex_binary_data = bitstring.ConstBitStream(filename=idex_packet_file)
        print("Data import completed, writing packet structures.")

        idex_parser = parser.PacketParser(idex_definition)
        idex_packet_generator = idex_parser.generator(idex_binary_data,
                                                    # skip_header_bits=64,
                                                    skip_header_bits=32,  # For sciData
                                                    show_progress=True,
                                                    yield_unrecognized_packet_errors=True)
    

        print("Packet structures written.")
        idex_binary_data.pos = 0
        idex_packet_generator = idex_parser.generator(idex_binary_data)
        self.data = {}
        self.header = {}
        self.analysis_flags = defaultdict(lambda: {
            'failed_fits': [],
            'saturated_channels': [],
            'notes': []
        })
        self.metadata_units: Dict[Tuple[int, str], str] = {}
        self._definitions_catalog = load_variable_definitions()
        evtnum = 0
        for pkt in idex_packet_generator:
            print(evtnum)
            if 'IDX__SCI0TYPE' in pkt.data:
                # print(evtnum)
                if pkt.data['IDX__SCI0TYPE'].raw_value == 1:
                    evtnum += 1
                    print(pkt.data)
                    # Iterate over all items in pkt.data and store them in the header
                    for key, item in pkt.data.items():
                        self.header[(evtnum, key)] = item.derived_value
                        print(f"{key} = {self.header[(evtnum, key)]}")
                    print(f"^*****Event header {evtnum}******^")
                    # sciEvtnum = bin(pkt.data['IDX__SCI0EVTNUM'].derived_value).replace('b', '')

                    # print(f"NBlocks = binary: {bin(pkt.data['IDX__TXHDRBLOCKS'].derived_value)} hex: {hex(pkt.data['IDX__TXHDRBLOCKS'].derived_value)}")
                    
                    # nBlocks = bin(pkt.data['IDX__TXHDRBLOCKS'].derived_value).replace('b', '')
                    # Extract the 17-22-bit integer (usually 8)
                    self.lspretrigblocks = (pkt.data['IDX__TXHDRBLOCKS'].derived_value >> 16) &  0b1111
                    # Extract the next 4-bit integer (usually 8)
                    self.lsposttrigblocks = (pkt.data['IDX__TXHDRBLOCKS'].derived_value >> 12) & 0b1111
                    # Extract the next 6 bits integer (usually 32)
                    self.hspretrigblocks = (pkt.data['IDX__TXHDRBLOCKS'].derived_value >> 6) & 0b111111
                    # Extract the first 6 bits (usually 32)
                    self.hsposttrigblocks = (pkt.data['IDX__TXHDRBLOCKS'].derived_value) & 0b111111

                    print("HS pre trig sampling blocks: ", self.hspretrigblocks)
                    print("LS pre trig sampling blocks: ", self.lspretrigblocks)
                    print("HS post trig sampling blocks: ", self.hsposttrigblocks)
                    print("LS post trig sampling blocks: ", self.lsposttrigblocks)
                    print(f"IDX__TXHDRHVPSHKCH01 = {pkt.data['IDX__TXHDRHVPSHKCH01'].derived_value}")
                    # Extract raw DN value for Voltage reading of Detector on HVPS Board (ADC CHnel 0)
                    self.header[(evtnum, 'detector_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH01'].derived_value) & 0b111111111111
                    print("Detector voltage = ", self.header[(evtnum, 'detector_voltage')])
                    # Extract raw DN value for Voltage reading of Sensor on HVPS Board (ADC CHnel 1)
                    self.header[(evtnum, 'sensor_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH01'].derived_value >> 16) & 0b111111111111
                    print("Sensor voltage = ", self.header[(evtnum, 'sensor_voltage')])
                    # HVPS Board signal "Target Voltage" (ADC CHnel 23)
                    self.header[(evtnum, 'target_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH23'].derived_value) & 0b111111111111
                    print("Target voltage = ", self.header[(evtnum, 'target_voltage')])
                    # HVPS Board signal "Reflectron Voltage" (ADC CHnel 23)
                    self.header[(evtnum, 'reflectron_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH23'].derived_value >> 16) & 0b111111111111
                    print("Reflectron voltage = ", self.header[(evtnum, 'reflectron_voltage')])
                    # HVPS Board signal "Rejection Voltage" (ADC CHnel 45)
                    self.header[(evtnum, 'rejection_voltage')] = (pkt.data['IDX__TXHDRHVPSHKCH45'].derived_value) & 0b111111111111
                    print("Rejection voltage = ", self.header[(evtnum, 'rejection_voltage')])
                    # HVPS Board signal "Current for the HVPS sensor" (ADC CHnel 45)
                    self.header[(evtnum, 'current_hvps_sensor')] = (pkt.data['IDX__TXHDRHVPSHKCH45'].derived_value >> 16) & 0b111111111111
                    print("Current for HVPS sensor = ", self.header[(evtnum, 'current_hvps_sensor')])
                    # HVPS Board signal "Positive current for the HVPS sensor" (ADC CHnel 67)
                    self.header[(evtnum, 'positive_current_hvps')] = (pkt.data['IDX__TXHDRHVPSHKCH67'].derived_value) & 0b111111111111
                    print("Positive current for HVPS sensor = ", self.header[(evtnum, 'positive_current_hvps')])
                    # HVPS Board signal "Negative current for the HVPS sensor" (ADC CHnel 67)
                    self.header[(evtnum, 'negative_current_hvps')] = (pkt.data['IDX__TXHDRHVPSHKCH67'].derived_value >> 16) & 0b111111111111
                    print("Negative current for HVPS sensor = ", self.header[(evtnum, 'negative_current_hvps')])
                    # LVPS Board signal "Voltage of +3.3V reference" (ADC CHnel 01)
                    self.header[(evtnum, 'voltage_3V3_ref')] = (pkt.data['IDX__TXHDRLVHK0CH01'].derived_value) & 0b111111111111
                    print("Voltage +3.3V reference = ", self.header[(evtnum, 'voltage_3V3_ref')])
                    # LVPS Board signal "Voltage of +3.3V operational reference" (ADC CHnel 01)
                    self.header[(evtnum, 'voltage_3V3_op_ref')] = (pkt.data['IDX__TXHDRLVHK0CH01'].derived_value >> 16) & 0b111111111111
                    print("Voltage +3.3V operational reference = ", self.header[(evtnum, 'voltage_3V3_op_ref')])
                    # LVPS Board signal "Voltage on -6V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'voltage_neg6V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH23'].derived_value) & 0b111111111111
                    print("Voltage -6V bus = ", self.header[(evtnum, 'voltage_neg6V_bus')])
                    # LVPS Board signal "Voltage on +6V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'voltage_pos6V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH23'].derived_value >> 16) & 0b111111111111
                    print("Voltage +6V bus = ", self.header[(evtnum, 'voltage_pos6V_bus')])
                    # LVPS Board signal "Voltage on +16V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'voltage_pos16V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH45'].derived_value) & 0b111111111111
                    print("Voltage +16V bus = ", self.header[(evtnum, 'voltage_pos16V_bus')])
                    # LVPS Board signal "Voltage on +3.3V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'voltage_pos3V3_bus')] = (pkt.data['IDX__TXHDRLVHK0CH45'].derived_value >> 16) & 0b111111111111
                    print("Voltage +3.3V bus = ", self.header[(evtnum, 'voltage_pos3V3_bus')])
                    # LVPS Board signal "Voltage on -5V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'voltage_neg5V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH67'].derived_value) & 0b111111111111
                    print("Voltage -5V bus = ", self.header[(evtnum, 'voltage_neg5V_bus')])
                    # LVPS Board signal "Voltage on +5V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'voltage_pos5V_bus')] = (pkt.data['IDX__TXHDRLVHK0CH67'].derived_value >> 16) & 0b111111111111
                    print("Voltage +5V bus = ", self.header[(evtnum, 'voltage_pos5V_bus')])
                    # LVPS Board signal "Current on +3.3V bus" (ADC CHnel 01)
                    self.header[(evtnum, 'current_3V3_bus')] = (pkt.data['IDX__TXHDRLVHK1CH01'].derived_value) & 0b111111111111
                    print("Current +3.3V bus = ", self.header[(evtnum, 'current_3V3_bus')])
                    # LVPS Board signal "Current on +16V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'current_16V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH23'].derived_value >> 16) & 0b111111111111
                    print("Current +16V bus = ", self.header[(evtnum, 'current_16V_bus')])
                    # LVPS Board signal "Current on +6V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'current_6V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH23'].derived_value) & 0b111111111111
                    print("Current +6V bus = ", self.header[(evtnum, 'current_6V_bus')])
                    # LVPS Board signal "Current on -6V bus" (ADC CHnel 23)
                    self.header[(evtnum, 'current_neg6V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH23'].derived_value >> 16) & 0b111111111111
                    print("Current -6V bus = ", self.header[(evtnum, 'current_neg6V_bus')])
                    # LVPS Board signal "Current on +5V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'current_5V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH45'].derived_value) & 0b111111111111
                    print("Current +5V bus = ", self.header[(evtnum, 'current_5V_bus')])
                    # LVPS Board signal "Current on -5V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'current_neg5V_bus')] = (pkt.data['IDX__TXHDRLVHK1CH45'].derived_value >> 16) & 0b111111111111
                    print("Current -5V bus = ", self.header[(evtnum, 'current_neg5V_bus')])
                    # LVPS Board signal "Current on +2.5V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'current_2V5_bus')] = (pkt.data['IDX__TXHDRLVHK1CH67'].derived_value) & 0b111111111111
                    print("Current +2.5V bus = ", self.header[(evtnum, 'current_2V5_bus')])
                    # LVPS Board signal "Current on -2.5V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'current_neg2V5_bus')] = (pkt.data['IDX__TXHDRLVHK1CH67'].derived_value >> 16) & 0b111111111111
                    print("Current -2.5V bus = ", self.header[(evtnum, 'current_neg2V5_bus')])

                    # LVPS Board signal "Current on the 1V POL" (ADC CHnel 01)
                    self.header[(evtnum, 'current_1V_pol')] = (pkt.data['IDX__TXHDRPROCHKCH01'].derived_value) & 0b111111111111
                    print("Current on the 1V POL = ", self.header[(evtnum, 'current_1V_pol')])
                    # LVPS Board signal "Current on the 1.9V POL" (ADC CHnel 01)
                    self.header[(evtnum, 'current_1.9V_pol')] = (pkt.data['IDX__TXHDRPROCHKCH01'].derived_value >> 16) & 0b111111111111
                    print("Current on the 1.9V POL = ", self.header[(evtnum, 'current_1.9V_pol')])
                    # LVPS Board signal "ProcBd Temperature 1" (ADC CHnel 23)
                    self.header[(evtnum, 'temperature_1')] = (pkt.data['IDX__TXHDRPROCHKCH23'].derived_value) & 0b111111111111
                    print("ProcBd Temperature 1 = ", self.header[(evtnum, 'temperature_1')])
                    # LVPS Board signal "ProcBd Temperature 2" (ADC CHnel 23)
                    self.header[(evtnum, 'temperature_2')] = (pkt.data['IDX__TXHDRPROCHKCH23'].derived_value >> 16) & 0b111111111111
                    print("ProcBd Temperature 2 = ", self.header[(evtnum, 'temperature_2')])
                    # LVPS Board signal "Voltage on 1V bus" (ADC CHnel 45)
                    self.header[(evtnum, 'voltage_1V_bus')] = (pkt.data['IDX__TXHDRPROCHKCH45'].derived_value) & 0b111111111111
                    print("Voltage on 1V bus = ", self.header[(evtnum, 'voltage_1V_bus')])
                    # LVPS Board signal "FPGA Temperature" (ADC CHnel 45)
                    self.header[(evtnum, 'fpga_temperature')] = (pkt.data['IDX__TXHDRPROCHKCH45'].derived_value >> 16) & 0b111111111111
                    print("FPGA Temperature = ", self.header[(evtnum, 'fpga_temperature')])
                    # LVPS Board signal "Voltage on 1.9V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'voltage_1.9V_bus')] = (pkt.data['IDX__TXHDRPROCHKCH67'].derived_value) & 0b111111111111
                    print("Voltage on 1.9V bus = ", self.header[(evtnum, 'voltage_1.9V_bus')])
                    # LVPS Board signal "Voltage on 3.3V bus" (ADC CHnel 67)
                    self.header[(evtnum, 'voltage_3.3V_bus')] = (pkt.data['IDX__TXHDRPROCHKCH67'].derived_value >> 16) & 0b111111111111
                    print("Voltage on 3.3V bus = ", self.header[(evtnum, 'voltage_3.3V_bus')])





                    mapping_dict = {
                        'detector_voltage': 'Last measurement in raw DN for HVPS Board signal “Detector Voltage”',
                        'sensor_voltage': 'Last measurement in raw DN for HVPS Board signal “Sensor Voltage"',
                        'target_voltage': 'Last measurement in raw DN for HVPS Board signal “Target Voltage”',
                        'reflectron_voltage': 'Last measurement in raw DN for HVPS Board signal “Reflectron Voltage”',
                        'rejection_voltage': 'Last measurement in raw DN for HVPS Board signal “Rejection Voltage”',
                        'current_hvps_sensor': 'Last measurement in raw DN for HVPS Board signal “Detector Current”',
                        'positive_current_hvps': 'Last measurement in raw DN for HVPS Board signal “Sensor IP”',
                        'negative_current_hvps': 'Last measurement in raw DN for HVPS Board signal “Sensor IN”',
                        'voltage_3V3_ref': 'Last measurement in raw DN for LVPS Board signal “P3.3VREF_HK”',
                        'voltage_3V3_op_ref': 'Last measurement in raw DN for LVPS Board signal “P3.3VREF_OP”',
                        'voltage_neg6V_bus': 'Last measurement in raw DN for LVPS Board signal “N6V”',
                        'voltage_pos6V_bus': 'Last measurement in raw DN for LVPS Board signal “P6V”',
                        'voltage_pos16V_bus': 'Last measurement in raw DN for LVPS Board signal “P16V”',
                        'voltage_pos3V3_bus': 'Last measurement in raw DN for LVPS Board signal “P3.3V”',
                        'voltage_neg5V_bus': 'Last measurement in raw DN for LVPS Board signal “N5V”',
                        'voltage_pos5V_bus': 'Last measurement in raw DN for LVPS Board signal “P5V”',
                        'current_3V3_bus': 'Last measurement in raw DN for LVPS Board signal “P3.3_IMON”',
                        'current_16V_bus': 'Last measurement in raw DN for LVPS Board signal “P16V_IMON”',
                        'current_6V_bus': 'Last measurement in raw DN for LVPS Board signal “P6V_IMON”',
                        'current_neg6V_bus': 'Last measurement in raw DN for LVPS Board signal “N6V_IMON”',
                        'current_5V_bus': 'Last measurement in raw DN for LVPS Board signal “P5V_IMON”',
                        'current_neg5V_bus': 'Last measurement in raw DN for LVPS Board signal “N5V_IMON”',
                        'current_2V5_bus': 'Last measurement in raw DN for LVPS Board signal “P2.5V_IMON”',
                        'current_neg2V5_bus': 'Last measurement in raw DN for LVPS Board signal “N2.5V_IMON”',
                        'spare_signal': 'Last measurement in raw DN for LVPS Board signal “Spare”',
                        'current_1V_pol':'Last measurement in raw DN for Processor Board signal “1V POL Current”',
                        'current_1.9V_pol':'Last measurement in raw DN for Processor Board signal “1.9V POL Current”',
                        'temperature_1': 'Last measurement in raw DN for Processor Board signal “ProcBd Temp1”',
                        'temperature_2': 'Last measurement in raw DN for Processor Board signal “ProcBd Temp2”',
                        'voltage_1V_bus': 'Last measurement in raw DN for Processor Board signal “1V Voltage”',
                        'fpga_temperature': 'Last measurement in raw DN for Processor Board signal “FPGA Temp”',
                        'voltage_1.9V_bus': 'Last measurement in raw DN for Processor Board signal “1.9V Voltage”',
                        'voltage_3.3V_bus': 'Last measurement in raw DN for Processor Board signal “3.3V Voltage”'
                    }
                    catalog = self._definitions_catalog
                    print("\n\n ***** Applying engineering unit conversions ***** \n\n")
                    for header_key, var_note in mapping_dict.items():
                        definition = catalog.find_by_var_notes(var_note)
                        if definition is None:
                            print(f"No variable definition found for note: {var_note}")
                            continue
                        raw_key = (evtnum, header_key)
                        if raw_key not in self.header:
                            print(f"Missing raw value for {header_key}")
                            continue
                        raw_value = self.header[raw_key]
                        try:
                            converted_value = definition.evaluate(raw_value)
                        except Exception as exc:
                            print(f"Failed to evaluate polynomial for {header_key}: {exc}")
                            continue
                        self.header[raw_key] = converted_value
                        self.metadata_units[raw_key] = definition.units or ""
                        print(
                            f"Converted {header_key} ({raw_value}) -> {converted_value} {definition.units}"
                        )

                     # Account for HS trigger delay
                    self.TOFdelay = pkt.data['IDX__TXHDRSAMPDELAY'].derived_value  # Last two bits are padding
                    # Mask to extract 10-bit values
                    mask = 0b1111111111
                    self.lgdelay = (self.TOFdelay) & mask # First 10 bits (0-9)
                    self.mgdelay = (self.TOFdelay >> 10) & mask # Next 10 bits (10-19)
                    self.hgdelay = (self.TOFdelay >> 20) & mask # Next 10 bits (20-29)
                    print(f"High gain delay = {self.hgdelay} samples.")
                    print(f"Mid gain delay = {self.mgdelay} samples.")
                    print(f"Low gain delay = {self.lgdelay} samples.")
                    if(pkt.data['IDX__TXHDRLSTRIGMODE'].derived_value!='DIS'):  # If this was a LS (Target Low Gain) trigger (DIS=disabled)
                        print(f"Low sampling trigger mode = {pkt.data['IDX__TXHDRLSTRIGMODE'].derived_value}")
                        self.Triggerorigin = 'LS' 
                        print("Low sampling trigger mode enabled.")
                        # Check the first 25th-bit integer for a coincidence trigger
                        # coincidence = (trigmode >> 24) &  0b1
                        # if(coincidence==1):
                            # self.Triggermode = 'Coincidence'
                        # else:
                            # self.Triggermode = 'Threshold'
                    print("High sampling trigger mode = ", pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value, pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value, pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value)
                    # Mask for extracting 11-bit and 10-bit values
                    mask_11_bit = 0b11111111111  # 11-bit mask
                    mask_10_bit = 0b1111111111   # 10-bit mask
                    if(pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value!=0):
                        print("Low gain TOF trigger mode enabled.")
                        self.Triggerorigin = 'LG'
                        # Extract the first 11 bits (bits 21-31)
                        minsamples = pkt.data['IDX__TXHDRLGTRIGCTRL1'].derived_value & mask_11_bit
                        # Extract the second 11 bits (bits 10-20)
                        maxsamples = (pkt.data['IDX__TXHDRLGTRIGCTRL1'].derived_value >> 11) & mask_11_bit
                        # Extract the last 10 bits (bits 0-9)
                        self.header[(evtnum, 'TriggerLevel')] = 2.89e-4*((pkt.data['IDX__TXHDRLGTRIGCTRL1'].derived_value >> 22) & mask_10_bit)
                        print(f"Trigger level = {self.header[(evtnum, 'TriggerLevel')]}")
                        if(pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value==1):
                            print("Threshold trigger mode enabled for low gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "LGThreshold"
                        elif(pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value==2):
                            print("Single pulse mode enabled for low gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "LGSinglePulse"
                        else:
                            print("Double pulse mode enabled for low gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "LGDoublePulse"
                    if(pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value!=0):
                        print("Mid gain TOF trigger mode enabled.")
                        self.Triggerorigin = 'MG'
                        # Extract the first 11 bits (bits 21-31)
                        minsamples = pkt.data['IDX__TXHDRMGTRIGCTRL1'].derived_value & mask_11_bit
                        # Extract the second 11 bits (bits 10-20)
                        maxsamples = (pkt.data['IDX__TXHDRMGTRIGCTRL1'].derived_value  >> 11) & mask_11_bit
                        # Extract the last 10 bits (bits 0-9)
                        self.header[(evtnum, 'TriggerLevel')] = 2.89e-4*((pkt.data['IDX__TXHDRMGTRIGCTRL1'].derived_value >> 22) & mask_10_bit)
                        print(f"Trigger level = {self.header[(evtnum, 'TriggerLevel')]}")
                        if(pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value==1):
                            print("Threshold trigger mode enabled for mid gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "MGThreshold"
                        elif(pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value==2):
                            print("Single pulse mode enabled for mid gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "MGSinglePulse"
                        else:
                            print("Double pulse mode enabled for mid gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "MGDoublePulse"
                    if(pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value!=0):
                        print("High gain trigger mode enabled.")
                        self.Triggerorigin = 'HG'
                        # Extract the first 11 bits (bits 21-31)
                        minsamples = pkt.data['IDX__TXHDRHGTRIGCTRL1'].derived_value & mask_11_bit
                        # Extract the second 11 bits (bits 10-20)
                        maxsamples = (pkt.data['IDX__TXHDRHGTRIGCTRL1'].derived_value  >> 11) & mask_11_bit
                        # Extract the last 10 bits (bits 0-9)
                        self.header[(evtnum, 'TriggerLevel')] = 2.89e-4*((pkt.data['IDX__TXHDRHGTRIGCTRL1'].derived_value >> 22) & mask_10_bit)
                        print(f"For {pkt.data['IDX__TXHDRHGTRIGCTRL1'].derived_value}, HG Trigger level = {self.header[(evtnum, 'TriggerLevel')]}, sample settings = {minsamples}, {maxsamples}")
                        if(pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value==1):
                            print("Threshold trigger mode enabled for high gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "HGThreshold"
                        elif(pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value==2):
                            print("Single pulse mode enabled for high gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "HGSinglePulse"
                        else:
                            print("Double pulse mode enabled for high gain channel.")
                            self.header[(evtnum, 'TriggerMode')] = "HGDoublePulse"
                    print(f"AID = {pkt.data['IDX__SCI0AID'].derived_value}")  # Instrument event number
                    print(f"Event number = {pkt.data['IDX__SCI0EVTNUM'].raw_value}")  # Event number out of how many events constitute the file
                    # print(f"Time = {pkt.data['IDX__SCI0TIME32'].derived_value}")  # Time in 20 ns intervals

                    print(f"Rice compression enabled = {bool(pkt.data['IDX__SCI0COMP'].raw_value)}")
                    compressed = bool(pkt.data['IDX__SCI0COMP'].raw_value)  # If we need to decompress the data

                    # self.header[evtnum][f"TimeIntervals"] = pkt.data['IDX__SCI0TIME32'].derived_value  # Store the number of 20 us intervals in the respective CDF "Time" variables
                    self.header[(evtnum, 'Timestamp')] = pkt.data['SHCOARSE'].derived_value + 20*(10**(-6))*pkt.data['SHFINE'].derived_value # Use this as the CDF epoch
                    print(f"Timestamp = {self.header[(evtnum, 'Timestamp')]} seconds since epoch (Midnight January 1st, 2012)")
                    # Convert to MST (UTC-7)
                    utc_time = datetime(2010, 1, 1, tzinfo=timezone.utc) + timedelta(seconds=self.header[(evtnum, 'Timestamp')])
                    # mst_offset = timedelta(hours=-7)
                    # mst_time = utc_time + mst_offset
                    print(f"Trigger time = {utc_time}")
                    self.header[(evtnum, 'Timestamp')] = utc_time.timestamp()


                if pkt.data['IDX__SCI0TYPE'].raw_value in [2, 4, 8, 16, 32, 64]:
                    # print(self.data.keys())
                    if (evtnum, pkt.data['IDX__SCI0TYPE'].raw_value) not in self.data.keys():  # If this is a new entry,
                        self.data.update({(evtnum, pkt.data['IDX__SCI0TYPE'].raw_value): pkt.data['IDX__SCI0RAW'].raw_value})
                    else:
                        self.data[(evtnum, pkt.data['IDX__SCI0TYPE'].raw_value)] += pkt.data['IDX__SCI0RAW'].raw_value


        # Parse the waveforms according to the scitype present (high gain and low gain CHnels encode waveform data differently).
        i = 1
        for scitype, waveform in self.data.items():
            if(compressed):  # If we need to decompress the data
                        print(waveform)
                        compressedFile = "test_compressed.txt"
                        dataFile = open(compressedFile, "wb")
                        index = 0
                        # print(f"||===waveform = {waveform}===||")
                        while index < len(waveform):
                            # Get 4 bytes (32 bits) from the 'waveform' binary string
                            data = waveform[index: index + 32]
                            # Convert the binary string to bytes using 'int' and 'to_bytes'
                            uint32 = int(data, 2).to_bytes(4, byteorder='big')
                            # Write the bytes to the file
                            dataFile.write(uint32)
                            index = index + 32
                        dataFile.close()
                        # print(waveform)
                        decompressor = RiceGolombDecompressor(waveform)
                        
                        if(scitype[1] < 12):  # LS
                            nsamples = 8*(self.lspretrigblocks + self.lsposttrigblocks)
                            # copy = gpt_rice_Decode(waveform, True, nsamples)
                            # copy = rice_Decode(compressedFile, f"test.txt", False, nsamples)
                            copy = decompressor.decompress(10)
                            waveform = copy
                        else:  # HS
                            nsamples = 512*(self.hspretrigblocks + self.hsposttrigblocks) # - pkt.data['IDX__TXHDRSAMPDELAY']
                            # copy = gpt_rice_Decode(waveform, True, nsamples)
                            # copy = rice_Decode(compressedFile, f"test.txt", True, nsamples)
                            copy = decompressor.decompress(12)
                            waveform = copy

            self.data[scitype] = parse_waveform_data(waveform, scitype[1])
        
        names = {2: "TOF H", 4: "TOF L", 8: "TOF M", 16: "Target L", 32: "Target H", 64: "Ion Grid"}
        datastore = {}
        for scitype, waveform in self.data.items():
            datastore[(scitype[0], names[scitype[1]])] = self.data[(scitype[0], scitype[1])]
        self.data = datastore
        self.numevents = evtnum
        # print(self.data.keys())

    # ||
    # ||
    # || Gather all of the events
    # || and plot them
    def plot_all_data(self, packets, fname: str):
        """Create per-event overview plots for all decoded channels.

        Parameters
        ----------
        packets
            Mapping of ``(event_number, channel_name)`` keys to decoded
            waveform arrays.
        fname
            Original filename used for naming the plot output directory.
        """

        fname = os.path.split(fname)[-1]
        # Create a folder to store the plots
        PlotFolder = os.path.join(os.getcwd(), f"Plots/{fname}")
        if os.path.exists(PlotFolder):  # If it exists, remove it
            shutil.rmtree(PlotFolder)
        os.makedirs(PlotFolder)

        # print("Number of packet items = ", len(packets.items()))
        def _build_stats_text(values):
            return (
                f"Min = {np.min(values)} [dN]\n"
                f"Avg = {np.mean(values):4.2f} [dN]\n"
                f"Std = {np.std(values):4.2f} [dN]\n"
                f"Max = {np.max(values)} [dN]"
            )

        fig, ax = plt.subplots(nrows=6)  # Make this general
        fig.set_size_inches(18.5, 10.5)
        fig.subplots_adjust(top=0.9, bottom=0.08, left=0.08, right=0.78, hspace=0.4)
        for i, (k, v) in enumerate(packets.items()):  # k[0] = Event number, k[1] = CHnel name, v=waveform data
            # fig = plt.figure(figsize=(17,12)) 
            # print(i%6)
            i=i%6  # We take modulo 6 so it is the same for each event
            x = np.linspace(0, len(v), len(v))  # Number of samples
            # Scale number of samples by ~4 ns (high rate) or ~250 ns (low rate) to get to time.
            if(i<=2):
                x *= 1/260  # HS
                self.hstime = x
            else:
                x *= 1/4.0625  # LS
                self.lstime = x

            # print("array length = ", len(v))

            
            print(f"Length of Channel {k[1]} = {len(v)}")
            # ax[i].fill_between(x, v, color='r')
            ax[i].set_xlabel(r"Time [$\mu$ s]", font="Times New Roman", fontsize=30, fontweight='bold')
            if(i<=2):
                self.lstriggertime = 8*(1/4.0625)*(self.lspretrigblocks+1) - (1/260)*self.hgdelay
                print(f"High sampling trigger time = {self.lstriggertime} microseconds.")
                self.hstime = self.hstime - self.lstriggertime
                ax[i].axvline(min(self.hstime)+self.lstriggertime, c="red", lw=2)

            else:
                self.hstriggertime = 512*(1/260)*(self.hspretrigblocks+1)  #  - (1/260)*self.hgdelay
                print(f"Low sampling trigger time = {self.hstriggertime} microseconds.")
                self.lstime = self.lstime-self.hstriggertime
                ax[i].axvline(min(self.lstime)+self.hstriggertime, c="red", lw=2)
                
            plt.suptitle(f"{fname} Event {k[0]}", font="Times New Roman", fontsize=30, fontweight='bold')
            # plt.tight_layout()

            if i==5:  #  End of the event, lets free up some memory
                ax[0].plot(self.hstime, packets[(k[0], "TOF L")])
                ax[0].set_ylabel("TOF L", font="Times New Roman", fontsize=15, fontweight='bold')
                text = _build_stats_text(packets[(k[0], "TOF L")])
                ax[0].text(
                    1.02,
                    0.98,
                    text,
                    fontsize=12,
                    va="top",
                    ha="left",
                    transform=ax[0].transAxes,
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                )
                # ax[0].set_xlim([0, 31.5])
                
                ax[1].plot(self.hstime, packets[(k[0], "TOF M")])
                ax[1].set_ylabel("TOF M", font="Times New Roman", fontsize=15, fontweight='bold')
                text = _build_stats_text(packets[(k[0], "TOF M")])
                ax[1].text(
                    1.02,
                    0.98,
                    text,
                    fontsize=12,
                    va="top",
                    ha="left",
                    transform=ax[1].transAxes,
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                )
                # ax[1].set_xlim([0, 31.5])
                
                ax[2].plot(self.hstime, packets[(k[0], "TOF H")])
                ax[2].set_ylabel("TOF H", font="Times New Roman", fontsize=15, fontweight='bold')
                text = _build_stats_text(packets[(k[0], "TOF H")])
                ax[2].text(
                    1.02,
                    0.98,
                    text,
                    fontsize=12,
                    va="top",
                    ha="left",
                    transform=ax[2].transAxes,
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                )
                # ax[2].set_xlim([0, 31.5])

                ax[3].plot(self.lstime, packets[(k[0], "Ion Grid")])
                ax[3].set_ylabel("Ion Grid", font="Times New Roman", fontsize=15, fontweight='bold')
                text = _build_stats_text(packets[(k[0], "Ion Grid")])
                ax[3].text(
                    1.02,
                    0.98,
                    text,
                    fontsize=12,
                    va="top",
                    ha="left",
                    transform=ax[3].transAxes,
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                )
                # ax[3].set_xlim([0, 126.5])
                
                if(self.header[(k[0], 'Timestamp')] < 494_733_600):  # If we are before September 27th, 2023 then we use the old definitions
                
                    ax[4].plot(self.lstime, packets[(k[0], "Target L")])
                    ax[4].set_ylabel("Target LG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = _build_stats_text(packets[(k[0], "Target L")])
                    ax[4].text(
                        1.02,
                        0.98,
                        text,
                        fontsize=12,
                        va="top",
                        ha="left",
                        transform=ax[4].transAxes,
                        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                    )
                    # ax[4].set_xlim([0, 126.5])
                    
                    ax[5].plot(self.lstime, packets[(k[0], "Target H")])
                    ax[5].set_ylabel("Target HG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = _build_stats_text(packets[(k[0], "Target H")])
                    ax[5].text(
                        1.02,
                        0.98,
                        text,
                        fontsize=12,
                        va="top",
                        ha="left",
                        transform=ax[5].transAxes,
                        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                    )
                    # ax[5].set_xlim([0, 126.5])

                else:
                    ax[4].plot(self.lstime, packets[(k[0], "Target H")])
                    ax[4].set_ylabel("Target HG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = _build_stats_text(packets[(k[0], "Target H")])
                    ax[4].text(
                        1.02,
                        0.98,
                        text,
                        fontsize=12,
                        va="top",
                        ha="left",
                        transform=ax[4].transAxes,
                        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                    )
                    # ax[4].set_xlim([0, 126.5])
                    
                    ax[5].plot(self.lstime, packets[(k[0], "Target L")])
                    ax[5].set_ylabel("Target LG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = _build_stats_text(packets[(k[0], "Target L")])
                    ax[5].text(
                        1.02,
                        0.98,
                        text,
                        fontsize=12,
                        va="top",
                        ha="left",
                        transform=ax[5].transAxes,
                        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                    )
                    # ax[5].set_xlim([0, 126.5])


                plt.savefig(os.path.join(PlotFolder, f"{fname}_Event_{k[0]}.png"), dpi=100)
                # plt.show()
                plt.close()
                del fig, ax
                fig, ax = plt.subplots(nrows=6)  # Make this general
                fig.set_size_inches(17.5, 10.5)

                
    
    # ||
    # ||
    # || Write the waveform data
    # || to an HDF5 file
    def write_to_hdf5(self, waveforms: dict, filename: str):
        """Persist decoded waveforms, instrument settings, and analysis products."""

        self.analysis_flags = defaultdict(lambda: {
            'failed_fits': [],
            'saturated_channels': [],
            'notes': []
        })

        output_path = Path(filename)
        if not output_path.is_absolute():
            output_dir = Path('HDF5')
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / output_path.name
        else:
            output_path.parent.mkdir(parents=True, exist_ok=True)

        if output_path.exists():
            output_path.unlink()

        def _adc_limit(channel_name: str) -> int:
            return 1023 if channel_name.startswith('TOF') else 4095

        def _as_array(value):
            if isinstance(value, str):
                return value
            if np.isscalar(value):
                return np.asarray(value)
            return np.asarray(value)

        conversion_factors = {
            'TOF H': 2.89e-4,
            'TOF M': 1.13e-2,
            'TOF L': 5.14e-4,
            'Ion Grid': 7.46e-4,
            'Target H': 1.63e-1,
            'Target L': 1.58e1
        }

        channel_units = {
            'TOF H': 'V',
            'TOF M': 'V',
            'TOF L': 'V',
            'Ion Grid': 'V',
            'Target H': 'V',
            'Target L': 'V'
        }

        mass_line_dtype = np.dtype([
            ('mass', np.float64),
            ('peak_index', np.int32),
            ('mu', np.float64),
            ('sigma', np.float64),
            ('lambda', np.float64),
            ('amplitude', np.float64),
            ('signal_amplitude', np.float64),
            ('abundance', np.float64),
            ('fit_success', np.bool_)
        ])

        with h5py.File(output_path, 'w') as h:
            for (evtnum, key), value in self.header.items():
                event_group = h.require_group(str(evtnum))
                self.analysis_flags[int(evtnum)]
                group_name = 'InstrumentSettings' if key.startswith('IDX__') else 'Metadata'
                target_group = event_group.require_group(group_name)
                dataset_value = _as_array(value)
                if isinstance(value, str):
                    dtype = h5py.string_dtype(encoding='utf-8')
                    dataset = create_or_replace_dataset(target_group, key, data=value, dtype=dtype)
                else:
                    dataset = create_or_replace_dataset(target_group, key, data=dataset_value)
                units = self.metadata_units.get((evtnum, key))
                if units is not None:
                    dataset.attrs['units'] = units
                if key == 'Timestamp':
                    create_or_replace_dataset(target_group, 'Epoch', data=np.asarray([value]))

            for (event_id, channel_name), samples in waveforms.items():
                event_group = h.require_group(str(event_id))
                waveform_group = event_group.require_group('Waveforms')
                raw_group = event_group.require_group('WaveformsDN')
                analysis_group = event_group.require_group('Analysis')
                time_group = event_group.require_group('Time')

                waveform_values = np.asarray(samples, dtype=float)
                raw_values = np.asarray(samples, dtype=np.int32)
                converted_values = waveform_values
                dataset_units = 'DN'

                if channel_name in conversion_factors:
                    converted_values = waveform_values * conversion_factors[channel_name]
                    dataset_units = channel_units[channel_name]

                dataset = create_or_replace_dataset(waveform_group, channel_name, data=converted_values)
                dataset.attrs['units'] = dataset_units
                if channel_name in conversion_factors:
                    dataset.attrs['conversion_factor'] = conversion_factors[channel_name]

                create_or_replace_dataset(raw_group, channel_name, data=raw_values)

                # -- Legacy compatibility -------------------------------------------------
                # Older quicklook tooling expects waveforms to live directly under the
                # event group (e.g. ``/1/TOF H``) and for the sampling axes to be stored
                # using the historic dataset names.  Preserve those layouts alongside the
                # structured hierarchy to keep both generations of software working.
                legacy_dataset = create_or_replace_dataset(event_group, channel_name, data=converted_values)
                legacy_dataset.attrs['units'] = dataset_units
                if channel_name in conversion_factors:
                    legacy_dataset.attrs['conversion_factor'] = conversion_factors[channel_name]
                if channel_name.startswith('TOF'):
                    create_dataset_if_not_exists(time_group, 'HighSampling', data=np.asarray(self.hstime, dtype=float))
                    create_dataset_if_not_exists(
                        event_group,
                        'Time (high sampling)',
                        data=np.asarray(self.hstime, dtype=float),
                    )
                else:
                    create_dataset_if_not_exists(time_group, 'LowSampling', data=np.asarray(self.lstime, dtype=float))
                    create_dataset_if_not_exists(
                        event_group,
                        'Time (low sampling)',
                        data=np.asarray(self.lstime, dtype=float),
                    )

                event_flags = self.analysis_flags[int(event_id)]
                if raw_values.size and raw_values.max() >= _adc_limit(channel_name):
                    event_flags['saturated_channels'].append(channel_name)

                channel_analysis_group = analysis_group.require_group(channel_name)

                if channel_name == 'TOF H':
                    _stretch, _shift, mass_scale = time2mass(converted_values, self.hstime)
                    create_dataset_if_not_exists(event_group, 'Mass', data=np.asarray(mass_scale))
                    peaks, _ = find_peaks(converted_values, prominence=.01)
                    if not peaks.size:
                        event_flags['notes'].append(f"{channel_name} no peaks identified")
                    mask = np.logical_and(self.hstime >= -7, self.hstime <= -5)
                    if mask.any():
                        baseline = np.mean(converted_values[mask])
                        sigma = np.std(converted_values[mask])
                    else:
                        baseline = np.nan
                        sigma = np.nan
                    peak_max = float(np.max(converted_values) - baseline) if np.isfinite(baseline) else float(np.max(converted_values))
                    snr = peak_max / sigma if np.isfinite(sigma) and sigma != 0.0 else np.nan
                    kappa = (float(np.mean([mass_scale[peak] - np.round(mass_scale[peak], 1) for peak in peaks]))
                             if peaks.size else np.nan)
                    create_dataset_if_not_exists(analysis_group, 'kappa', data=np.asarray([kappa]))
                    create_dataset_if_not_exists(analysis_group, 'SNR', data=np.asarray([snr]))

                    fit_results = []
                    mass_line_records = []
                    for peak in peaks:
                        start = max(0, peak - 5)
                        end = min(len(converted_values), peak + 6)
                        x_slice = self.hstime[start:end]
                        y_slice = converted_values[start:end]
                        param, param_cov, sig_amp, fitted_curve = FitEMG(x_slice, y_slice)
                        if isinstance(param, (np.ndarray, list, tuple)):
                            param = np.asarray(param, dtype=float)
                            area = calculate_area_under_emg(x_slice, param)
                            fit_results.append((param, sig_amp, x_slice, fitted_curve, area))
                            amplitude_param, mu_param, sigma_param, lam_param = param
                            mass_line_records.append({
                                'mass': float(mass_scale[peak]),
                                'peak_index': int(peak),
                                'mu': float(mu_param),
                                'sigma': float(sigma_param),
                                'lambda': float(lam_param),
                                'amplitude': float(amplitude_param),
                                'signal_amplitude': float(sig_amp),
                                'abundance': float(area),
                                'fit_success': True
                            })
                        else:
                            event_flags['failed_fits'].append(f"{channel_name} peak {int(peak)}")
                            fit_results.append(None)
                            mass_line_records.append({
                                'mass': float(mass_scale[peak]),
                                'peak_index': int(peak),
                                'mu': np.nan,
                                'sigma': np.nan,
                                'lambda': np.nan,
                                'amplitude': np.nan,
                                'signal_amplitude': np.nan,
                                'abundance': np.nan,
                                'fit_success': False
                            })

                    masses_group = channel_analysis_group.require_group('Masses')
                    for result in fit_results:
                        if result is None:
                            continue
                        param, sig_amp, x_slice, fitted_curve, area = result
                        amplitude_param, mu_param, sigma_param, lam_param = param
                        create_dataset_if_not_exists(
                            masses_group,
                            f"{mu_param}FitParams",
                            data=np.asarray([amplitude_param, mu_param, sigma_param, lam_param])
                        )
                        create_dataset_if_not_exists(
                            masses_group,
                            f"{mu_param}AreaUnderFit",
                            data=np.asarray([area])
                        )

                    mass_line_dataset = (np.asarray(mass_line_records, dtype=mass_line_dtype)
                                         if mass_line_records
                                         else np.empty(0, dtype=mass_line_dtype))
                    create_or_replace_dataset(channel_analysis_group, 'MassLines', data=mass_line_dataset)

                if channel_name in ['Target L', 'Target H', 'Ion Grid']:
                    try:
                        param, param_cov, sig_amp, fit_time, fit_signal, fit_result = FitTargetSignal(self.lstime, samples)
                    except Exception as exc:
                        event_flags['failed_fits'].append(f"{channel_name} target fit error: {exc}")
                        param = np.full(5, np.nan)
                        sig_amp = np.nan
                        fit_time = None
                        fit_signal = None
                        fit_result = None

                    if param is not None:
                        param_array = np.asarray(param, dtype=float)
                        create_dataset_if_not_exists(channel_analysis_group, 'FitParams', data=param_array)
                        create_dataset_if_not_exists(
                            analysis_group,
                            f"{channel_name}FitParams",
                            data=param_array,
                        )
                    create_dataset_if_not_exists(channel_analysis_group, 'MassEstimate', data=np.asarray([sig_amp]))
                    create_dataset_if_not_exists(channel_analysis_group, 'ImpactCharge', data=np.asarray([sig_amp]))
                    create_dataset_if_not_exists(
                        analysis_group,
                        f"{channel_name}MassEstimate",
                        data=np.asarray([sig_amp]),
                    )
                    create_dataset_if_not_exists(
                        analysis_group,
                        f"{channel_name}ImpactCharge",
                        data=np.asarray([sig_amp]),
                    )
                    if fit_time is not None and fit_signal is not None:
                        fit_time_array = np.asarray(fit_time, dtype=float)
                        fit_signal_array = np.asarray(fit_signal, dtype=float)
                        create_dataset_if_not_exists(channel_analysis_group, 'FitTime', data=fit_time_array)
                        create_dataset_if_not_exists(channel_analysis_group, 'FitData', data=fit_signal_array)
                        create_dataset_if_not_exists(
                            analysis_group,
                            f"{channel_name}FitTime",
                            data=fit_time_array,
                        )
                        create_dataset_if_not_exists(
                            analysis_group,
                            f"{channel_name}FitData",
                            data=fit_signal_array,
                        )
                    if fit_result is not None:
                        fit_result_array = np.asarray(fit_result, dtype=float)
                        create_dataset_if_not_exists(channel_analysis_group, 'FitResult', data=fit_result_array)
                        create_dataset_if_not_exists(
                            analysis_group,
                            f"{channel_name}FitResult",
                            data=fit_result_array,
                        )

            string_dtype = h5py.string_dtype(encoding='utf-8')
            for event_id, flags in self.analysis_flags.items():
                analysis_group = h.require_group(str(event_id)).require_group('Analysis')
                flags_group = analysis_group.require_group('Flags')
                failed = np.array(sorted(set(flags['failed_fits'])), dtype=string_dtype)
                saturated = np.array(sorted(set(flags['saturated_channels'])), dtype=string_dtype)
                notes = np.array(sorted(set(flags['notes'])), dtype=string_dtype)
                create_or_replace_dataset(flags_group, 'FailedFits', data=failed)
                create_or_replace_dataset(flags_group, 'SaturatedChannels', data=saturated)
                create_or_replace_dataset(flags_group, 'Notes', data=notes)

# ||
# ||
# || Parse the high sampling rate data, this
# || should be 10-bit blocks
def _get_array_module():
    """Return the array module (NumPy or CuPy) for vectorised bit decoding."""

    if cp is not None:
        try:  # pragma: no cover - only exercised when a GPU is present
            cp.cuda.runtime.getDevice()
        except Exception:
            return np
        else:
            return cp
    return np


def _decode_bitpacked_waveform(
    waveform_raw: str,
    skip_bits: int,
    bits_per_sample: int,
    samples_per_frame: int,
    trim_tail: int = 0,
):
    """Decode bit-packed waveforms using vectorised NumPy/CuPy operations.

    Parameters
    ----------
    waveform_raw
        Binary string extracted from the packet payload.
    skip_bits
        Number of pad bits to discard from the beginning of each frame.
    bits_per_sample
        Width of each encoded sample in bits.
    samples_per_frame
        Number of samples stored in each frame after the pad bits.
    trim_tail
        Optional number of decoded samples to drop from the tail.
    """

    if not waveform_raw:
        return []

    xp = _get_array_module()

    ascii_buffer = np.frombuffer(waveform_raw.encode("ascii"), dtype=np.uint8) - 48
    if xp is cp:
        bits = xp.asarray(ascii_buffer)
    else:
        bits = ascii_buffer

    frame_width = skip_bits + bits_per_sample * samples_per_frame
    if frame_width == 0:
        return []

    frame_count = bits.size // frame_width
    if frame_count == 0:
        return []

    trimmed = bits[: frame_count * frame_width]
    frames = trimmed.reshape(frame_count, frame_width)
    data_bits = frames[:, skip_bits:]
    samples = data_bits.reshape(-1, bits_per_sample)

    weights = xp.power(2, xp.arange(bits_per_sample - 1, -1, -1, dtype=xp.int64))
    decoded = samples.dot(weights)

    raw_length = int(decoded.size)
    if trim_tail:
        decoded = decoded[:-trim_tail]

    if xp is cp:  # pragma: no cover - exercised only when a GPU is present
        decoded = cp.asnumpy(decoded)

    print(raw_length)
    return decoded.astype(np.int32, copy=False).tolist()


def parse_hs_waveform(waveform_raw: str):
    """Parse a binary string representing a high sampling-rate waveform."""

    return _decode_bitpacked_waveform(
        waveform_raw,
        skip_bits=2,
        bits_per_sample=10,
        samples_per_frame=3,
        trim_tail=4,
    )

# ||
# ||
# || Parse the low sampling rate data, this
# || should be 12-bit blocks
def parse_ls_waveform(waveform_raw: str):
    """Parse a binary string representing a low sampling-rate waveform."""

    return _decode_bitpacked_waveform(
        waveform_raw,
        skip_bits=8,
        bits_per_sample=12,
        samples_per_frame=2,
    )

# ||
# ||
# || Use the SciType flag to determine the sampling rate of
# || the data we are trying to parse
def parse_waveform_data(waveform: str, scitype: int):
    """Parse a waveform string using the SCI type to select the decoder.

    Parameters
    ----------
    waveform
        Binary waveform payload.
    scitype
        SCI type identifier describing which decoding scheme to apply.

    Returns
    -------
    list[int]
        Decoded waveform data numbers.
    """
    print(f'Parsing waveform for scitype={scitype}')
    if scitype in (2, 4, 8):
        return parse_hs_waveform(waveform)
    else:
        return parse_ls_waveform(waveform)

# ||
# ||
# || Write the waveform data 
# || to CDF files
def write_to_cdf(packets):
    """Write decoded packets to a CDF file using the master template.

    Parameters
    ----------
    packets
        Initialized :class:`IDEXEvent` containing waveform arrays and
        metadata. The function mirrors the template structure contained in
        ``imap_idex_l0-raw_0000000_v01.cdf``.
    """

    cdf_master = cdfread.CDF('imap_idex_l0-raw_0000000_v01.cdf')
    if (cdf_master.file != None):
    # Get the cdf's specification
        info=cdf_master.cdf_info()
        cdf_file=cdfwrite.CDF('./IDEX_SSIM.cdf',cdf_spec=info,delete=True)
    # if (cdf_file.file == None):
    #     cdf_master.close()
    #     raise OSError('Problem writing file.... Stop')

    # Get the global attributes
    globalaAttrs=cdf_master.globalattsget(expand=True)
    # Write the global attributes
    cdf_file.write_globalattrs(globalaAttrs)
    zvars=info['zVariables']
    print('no of zvars=',len(zvars))
    # Loop thru all the zVariables --> What are zvars vs rvars?
    for x in range (0, len(zvars)):
        # Get the variable's specification
        varinfo=cdf_master.varinq(zvars[x])
        print('Z =============>',x,': ', varinfo['Variable'])


# Z =============> 0 :  Epoch
# Z =============> 1 :  IDEX_Trigger
# Z =============> 2 :  TOF_Low
# Z =============> 3 :  TOF_Mid
# Z =============> 4 :  TOF_High
# Z =============> 5 :  Time_Low_SR
# Z =============> 6 :  Time_High_SR
# Z =============> 7 :  Target_Low
# Z =============> 8 :  Target_High
# Z =============> 9 :  Ion_Grid

        if(varinfo['Variable']=="Epoch"):
            vardata = None
        if(varinfo['Variable']=="IDEX_Trigger"):
            vardata = packets.header[(1,"Timestamp")]
        if(varinfo['Variable']=="TOF_Low"):
            print(len(np.array(packets.data[(1,"TOF L")])))
            vardata = np.array(packets.data[(1,"TOF L")], float)
        if(varinfo['Variable']=="TOF_Mid"):
            vardata = np.array(packets.data[(1,"TOF M")])
        if(varinfo['Variable']=="TOF_High"):
            vardata = np.array(packets.data[(1,"TOF H")])
        if(varinfo['Variable']=="Time_Low_SR"):
            vardata = np.linspace(0, len(packets.data[(1,"Target L")]), len(len(packets.data[(1,"Target L")])))
        if(varinfo['Variable']=="Time_High_SR"):
            vardata = np.linspace(0, len(packets.data[(1,"TOF L")]), len(len(packets.data[(1,"Target L")])))
        if(varinfo['Variable']=="Target_Low"):
            vardata = np.array(packets.data[(1,"Target L")])
        if(varinfo['Variable']=="Target_High"):
            vardata = np.array(packets.data[(1,"Target H")])
        if(varinfo['Variable']=="Ion_Grid"):
            vardata = np.array(packets.data[(1,"Ion Grid")])
        # Get the variable's attributes
        varattrs=cdf_master.varattsget(zvars[x], expand=True)
        if (varinfo['Sparse'].lower() == 'no_sparse'):
            # A variable with no sparse records... get the variable data
            # vardata= None
            # Create the zVariable, write out the attributes and data
            cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
        else:
            # A variable with sparse records...
            # data is in this form [physical_record_numbers, data_values]
            # physical_record_numbers (0-based) contains the real record
            # numbers. For example, a variable has only 3 physical records
            # at [0, 5, 10]:
            varrecs=[0,5,10]

            # vardata=None  # np.asarray([.,.,.,..])
            # Create the zVariable, and optionally write out the attributes
            # and data
            cdf_file.write_var(varinfo, var_attrs=varattrs,
                        var_data=[varrecs,vardata])
    rvars=info['rVariables']
    print('no of rvars=',len(rvars))
    # Loop thru all the rVariables
    for x in range (0, len(rvars)):
        varinfo=cdf_master.varinq(rvars[x])
        print('R =============>',x,': ', varinfo['Variable'])
        varattrs=cdf_master.varattsget(rvars[x], expand=True)
        if (varinfo['Sparse'].lower() == 'no_sparse'):
            vardata=None
            # Create the rVariable, write out the attributes and data
            cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
        else:
            varrecs= None  # [.,.,.,..]
            vardata= None  # np.asarray([.,.,.,..])
            cdf_file.write_var(varinfo, var_attrs=varattrs,
                        var_data=[varrecs,vardata])
    cdf_master.close()
    cdf_file.close()

# || Test code: Import file and write the relevant data to an hdf5 file
if __name__ == "__main__":
    # Initalize parsing object to pass filename
    aparser = argparse.ArgumentParser(
        description=(
            "Decode an IDEX L0 binary file from the Data/ directory and "
            "produce higher level quicklook outputs."
        )
    )
    aparser.add_argument(
        "--file",
        "-f",
        type=str,
        required=True,
        help="Path to the raw L0 binary file (typically within Data/).",
    )
    args = aparser.parse_args()

    packets = IDEXEvent(args.file)
    # print(packets.data.keys())
    try:
        packets.plot_all_data(packets.data, args.file)
    except Exception as e:
        print(e)
    packets.write_to_hdf5(packets.data, args.file+'.h5')
    # write_to_cdf(packets)

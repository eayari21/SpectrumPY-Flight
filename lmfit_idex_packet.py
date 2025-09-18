#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

"""
A Python object to store IDEX packets.
__author__      = Ethan Ayari & Gavin Medley, 
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
"""

# || Python libraries
import argparse
import os
import socket
import bitstring
import h5py
import shutil
import lmfit
import struct
import matplotlib.pyplot as plt
import pandas as pd
# try:
#     plt.style.use("seaborn-pastel")
# except:
#     plt.style.use("seaborn-v0_8-pastel")
plt.style.use("seaborn-pastel")
import numpy as np

from datetime import datetime, timedelta, timezone

from scipy.optimize import curve_fit
from scipy.signal import detrend, butter, filtfilt, find_peaks
from scipy.integrate import quad
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
    return P1 + np.heaviside(x-P0, 0) * ( P4 * (1.0 - np.exp(-(x-P0)/P5)) * np.exp( -(x-P0)/P6))

# Define the EMG function
def EMG(x, mu, sigma, lam):
    prefactor = lam / 2
    exponent = np.exp((lam / 2) * (2 * mu + lam * sigma**2 - 2 * x))
    erfc_part = erfc((mu + lam * sigma**2 - x) / (np.sqrt(2) * sigma))
    return prefactor * exponent * erfc_part

# Function to calculate the area under the EMG fit curve
def calculate_area_under_emg(x_slice, param):
    if(type(param) is not int):
        # Extract EMG fit parameters: mu, sigma, lam
        mu, sigma, lam = param

        # Perform numerical integration using scipy.integrate.quad
        area, error = quad(EMG, x_slice[0], x_slice[-1], args=(mu, sigma, lam))
        
        return area
    else:
        return 0.0

# Helper function to apply the polynomial transformation
def apply_polynomial(coeffs, X):
    # Compute the value using the polynomial formula
    return sum(coeffs[i] * (X ** i) for i in range(len(coeffs)))

# Helper function to create dataset if it doesn't exist
def create_dataset_if_not_exists(hdf5_file, dataset_path, data):
    if dataset_path in hdf5_file:
        print(f"Dataset '{dataset_path}' already exists. Skipping creation.")
    else:
        hdf5_file.create_dataset(dataset_path, data=data)

# Fit routine for EMG
def FitEMG(time, amplitude):
    x = np.asarray(time)
    y = np.asarray(amplitude)

    # || Initial Guess for the parameters of the EMG
    mu_guess = x[np.argmax(y)]  # Initial guess for the mean
    sigma_guess = np.std(x) / 10  # Initial guess for standard deviation
    lam_guess = 1 / (x[-1] - x[0])  # Initial guess for decay rate

    p0 = [mu_guess, sigma_guess, lam_guess]  # Initial parameter guesses

    # Fit the data using curve_fit
    try:
        param, param_cov = curve_fit(EMG, x, y, p0=p0, maxfev=100_000)

        # Generate the fitted curve
        result = EMG(x, *param)
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
    x = np.asarray(time)
    y = np.asarray(targetAmp)

    # || Select only raw noise (where we know the signal is not)
    mask = np.logical_and(time >= -7, time <= -5)
    print(mask)

    baselineraw = y[mask]
    baselinedomain = time[mask]

    # || Remove Linear Background
    try:
        # Linear background fit using lmfit
        slopeguess = 0
        def linear_model(params, x):
            slope = params['slope']
            intercept = params['intercept']
            return slope * x + intercept

        linear_params = lmfit.Parameters()
        linear_params.add('slope', value=slopeguess)
        linear_params.add('intercept', value=0)

        linear_model_fn = lambda params, x: linear_model(params, x)
        linear_fit = lmfit.minimize(lambda params: linear_model_fn(params, baselinedomain) - baselineraw, linear_params)

        linparam = linear_fit.params
        linearbase = linear_model(linparam, time)
        y = detrend(y)
    except:
        print(f"Linear background not found.")
        linearbase = None

    # || Remove Sine Wave Background
    try:
        baselinedelined = y[mask]
        def sine_model(params, x):
            amplitude = params['amplitude']
            frequency = params['frequency']
            phase = params['phase']
            return amplitude * np.sin(2 * np.pi * frequency * x + phase)

        sine_params = lmfit.Parameters()
        sine_params.add('amplitude', value=max(baselinedelined))
        sine_params.add('frequency', value=7000)
        sine_params.add('phase', value=45)

        sine_model_fn = lambda params, x: sine_model(params, x)
        sine_fit = lmfit.minimize(lambda params: sine_model_fn(params, baselinedomain) - baselinedelined, sine_params)

        sineparam = sine_fit.params
        sinebase = sine_model(sineparam, time)
        y -= sinebase
        y = butter_lowpass_filter(y, time)

    except:
        print(f"Sinusoidal background not found.")
        sinebase = None

    # print(f"idx = {idx}")
    mask = np.logical_and(time >= -7, time <= -5)

    x = time
    y = y
    pre = -2.0 # Before image charge

    # Assuming x and y are numpy arrays and `pre` is the threshold
    yBaseline = np.where(x < pre, y, np.nan)  # Set elements of y to nan where x >= pre
    yImage = y[(x >= pre) & (x < 0.0)]
    ionError = np.std(yBaseline)
    ionErrorVector = pd.DataFrame([np.nan] * len(yBaseline))
    ionMean = yBaseline.mean()
    yBaseline -= ionMean
    parameters = []

    ionTime = np.array([float(x) for x in x])
    ionAmp = np.array([float(y) for y in y])

    print(f"Calculating Target Fit.")

    # || Initial Guess for the parameters of the ion grid signal
    t0 = 0.0  # P[0] time of impact
    c = 0.0   # P[1] Constant offset
    A = max(ionAmp)  # P[2] amplitude (v)
    t1 = 3.71  # P[3] rise time (s)
    t2 = 37.1  # P[4] discharge time (s)

    # Create the model function for the ion grid
    def ion_grid_model(params, x):
        t0 = params['t0']
        c = params['c']
        A = params['A']
        t1 = params['t1']
        t2 = params['t2']
        return c + A * (1 - np.exp(-(x - t0) / t1)) * np.exp(-(x - t0) / t2)

    # Define parameters for lmfit
    ion_params = lmfit.Parameters()
    ion_params.add('t0', value=t0)
    ion_params.add('c', value=c)
    ion_params.add('A', value=A)
    ion_params.add('t1', value=t1)
    ion_params.add('t2', value=t2)

    try:
        # Fit the ion grid model to the data
        ion_fit = lmfit.minimize(lambda params: ion_grid_model(params, ionTime) - ionAmp, ion_params)

        # Get the fit report, which includes chi-squared and reduced chi-squared
        fit_report = lmfit.fit_report(ion_fit)

        # Print the fit report
        print(fit_report)

        # Get the fit parameters and covariance matrix
        param = ion_fit.params
        param_cov = ion_fit.covar

        # Get the chi-squared and reduced chi-squared from the fit result
        chi_squared = ion_fit.chisqr
        reduced_chi_squared = ion_fit.redchi

        # Calculate the resulting signal amplitude
        result = ion_grid_model(ion_fit.params, ionTime)

        sig_amp = max(result) - yBaseline.mean()

        return param, param_cov, sig_amp, chi_squared, reduced_chi_squared, result, ionTime 
    except Exception as e:
        print(f"Target signal failed due to exception {e}")
        return [0,0,0,0,0], 0, 0, 0, 0

# ||
# ||
# || Generator object from LASP packets
# || to read in the data
class IDEXEvent:
    def __init__(self, filename: str):
        """Test parsing a real XTCE document"""
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
        self.header={}
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






                    # Define the coefficients that are the same for every variable
                    coefficients = ['CO', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']

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

                    # Create a reverse dictionary
                    reverse_mapping_dict = {value: key for key, value in mapping_dict.items()}

                    # Read in Scott K's instrument settings conversions
                    settings_df = pd.read_excel("IDEX CDF Variable Definitions.xlsx")
                    # Normalize Var_notes by converting curly quotes to straight quotes and removing NaN values
                    # settings_df['Var_notes'] = settings_df['Var_notes'].replace(np.nan, '', regex=True)  # Replace NaN with empty strings
                    # settings_df['Var_notes'] = settings_df['Var_notes'].str.replace('“', '"').str.replace('”', '"')  # Replace curly quotes with straight quotes
                    # settings_df['Var_notes'] = settings_df['Var_notes'].str.strip()  # Remove any leading/trailing spaces

                    # print("Var notes: ", [note for note in settings_df["Var_notes"]])

                    # print("Coefficients: ", settings_df.columns)

                    # Step 1: Create the mapping from Var_notes to row indices
                    var_to_row = {var_note: idx for idx, var_note in enumerate(settings_df['Var_notes'])}

                    print(f"var_to_row = {var_to_row}")

                    print(f"\n \n ***** var_name being converted for each instrument setting ***** \n \n")

                    for var_name, row_idx in var_to_row.items():
                        try:
                            print("Matching row ", len(settings_df.iloc[[row_idx]].columns))  # Print the entire row for the problematic variable

                            # Extract the polynomial coefficients from the spreadsheet
                            coeffs = settings_df.iloc[row_idx][coefficients].values  # Access the row by index and then columns by labels
                            print(f"coeffs for {var_name}= {coeffs}")
                            target_value = var_name
                            print(f"var_name = {var_name}")


                            # Get the corresponding key
                            var_name = reverse_mapping_dict.get(target_value)
                            print(f"var_name = {var_name}")
                            
                            # Get the raw DN value for this variable (from your script)
                            X = self.header[(evtnum, var_name)]

                            print(f"Header info = {X}")
                            
                            # Apply the polynomial transformation
                            transformed_value = apply_polynomial(coeffs, X)
                            self.header[(evtnum, var_name)] = transformed_value
                            
                            print(f"Transformed value for {var_name} = {transformed_value}")
                        
                        except KeyError as e:
                            print(f"KeyError: {e} - Could not find entry for {var_name}, {row_idx}. Please check the Var_notes or variable mapping.")
                        except Exception as e:
                            print(f"An error occurred: {e}")


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
        fname = os.path.split(fname)[-1]
        # Create a folder to store the plots
        PlotFolder = os.path.join(os.getcwd(), f"Plots/{fname}")
        if os.path.exists(PlotFolder):  # If it exists, remove it
            shutil.rmtree(PlotFolder)
        os.makedirs(PlotFolder)

        # print("Number of packet items = ", len(packets.items()))
        fig, ax = plt.subplots(nrows=6)  # Make this general
        fig.set_size_inches(18.5, 10.5)
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
                
            plt.subplots_adjust(bottom=0.2)


            plt.suptitle(f"{fname} Event {k[0]}", font="Times New Roman", fontsize=30, fontweight='bold')
            # plt.tight_layout()

            if i==5:  #  End of the event, lets free up some memory
                ax[0].plot(self.hstime, packets[(k[0], "TOF L")])
                ax[0].set_ylabel("TOF L", font="Times New Roman", fontsize=15, fontweight='bold')
                text = f'Min = {min(packets[(k[0], "TOF L")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "TOF L")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "TOF L")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "TOF L")])} [dN]'
                ax[0].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[0].transAxes)
                # ax[0].set_xlim([0, 31.5])
                
                ax[1].plot(self.hstime, packets[(k[0], "TOF M")])
                ax[1].set_ylabel("TOF M", font="Times New Roman", fontsize=15, fontweight='bold')
                text = f'Min = {min(packets[(k[0], "TOF M")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "TOF M")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "TOF M")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "TOF M")])} [dN]'
                ax[1].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[1].transAxes)
                # ax[1].set_xlim([0, 31.5])
                
                ax[2].plot(self.hstime, packets[(k[0], "TOF H")])
                ax[2].set_ylabel("TOF H", font="Times New Roman", fontsize=15, fontweight='bold')
                text = f'Min = {min(packets[(k[0], "TOF H")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "TOF H")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "TOF H")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "TOF H")])} [dN]'
                ax[2].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[2].transAxes)
                # ax[2].set_xlim([0, 31.5])

                ax[3].plot(self.lstime, packets[(k[0], "Ion Grid")])
                ax[3].set_ylabel("Ion Grid", font="Times New Roman", fontsize=15, fontweight='bold')
                text = f'Min = {min(packets[(k[0], "Ion Grid")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Ion Grid")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Ion Grid")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Ion Grid")])} [dN]'
                ax[3].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[3].transAxes)
                # ax[3].set_xlim([0, 126.5])
                
                if(self.header[(k[0], 'Timestamp')] < 494_733_600):  # If we are before September 27th, 2023 then we use the old definitions
                
                    ax[4].plot(self.lstime, packets[(k[0], "Target L")])
                    ax[4].set_ylabel("Target LG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = f'Min = {min(packets[(k[0], "Target L")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Target L")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Target L")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Target L")])} [dN]'
                    ax[4].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[4].transAxes)
                    # ax[4].set_xlim([0, 126.5])
                    
                    ax[5].plot(self.lstime, packets[(k[0], "Target H")])
                    ax[5].set_ylabel("Target HG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = f'Min = {min(packets[(k[0], "Target H")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Target H")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Target H")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Target H")])} [dN]'
                    ax[5].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[5].transAxes)
                    # ax[5].set_xlim([0, 126.5])

                else:
                    ax[4].plot(self.lstime, packets[(k[0], "Target H")])
                    ax[4].set_ylabel("Target HG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = f'Min = {min(packets[(k[0], "Target H")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Target H")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Target H")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Target H")])} [dN]'
                    ax[4].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[4].transAxes)
                    # ax[4].set_xlim([0, 126.5])
                    
                    ax[5].plot(self.lstime, packets[(k[0], "Target L")])
                    ax[5].set_ylabel("Target LG", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = f'Min = {min(packets[(k[0], "Target L")])} [dN]'+ '\n'+ f'Avg={np.mean(packets[(k[0], "Target L")]): 4.2f} [dN]'+ '\n' + f'Std={np.std(packets[(k[0], "Target L")]): 4.2f} [dN]'+ '\n' + f'Max = {max(packets[(k[0], "Target L")])} [dN]'
                    ax[5].text(1.125, 0.85, text, fontsize=12, va="top", ha="right", transform=ax[5].transAxes)
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
        os.chdir('./HDF5/')
        filename = os.path.split(filename)[-1]  # Just get the name of the file
        # Prepend HDF5 folder to filename

        # print(waveforms.keys())
        # print(waveforms.values())

        if os.path.exists(filename):
            os.remove(filename)
        h = h5py.File(filename,'w')

        for (evtnum, key), value in self.header.items():
            # Skip fields starting with "IDX__" for waveforms
            if "IDX__" in key:
                print(f"Skipping header field {key}")
                continue
            # Create the dataset for each item in self.header
            dataset_path = f"/{evtnum}/Metadata/{key}"
            
            # Check if the dataset exists before creating it
            if dataset_path not in h:
                # Check if the value is a string and handle accordingly
                if isinstance(value, str):
                    # Use string_dtype to handle string data
                    dtype = h5py.string_dtype(encoding='utf-8')
                    h.create_dataset(dataset_path, data=value, dtype=dtype)
                else:
                    # Convert value to numpy array and write normally for non-strings
                    h.create_dataset(dataset_path, data=np.array(value))
                
                print(f"Created dataset: {dataset_path} with value: {value}")

        # Conversion factors based on channel type
        conversion_factors = {
            'TOF H': 2.89e-4,
            'TOF M': 1.13e-2,
            'TOF L': 5.14e-4,
            'Ion Grid': 7.46e-4,
            'Target H': 1.63e-1,
            'Target L': 1.58e1
        }

        for k, v in waveforms.items():
            # Convert degrees to radians for RA and Dec
            ra_values = np.deg2rad(np.random.uniform(0, 15, size=1))  # e.g., 100 sample points
            dec_values = np.deg2rad(np.random.uniform(-5, 5, size=1))


            # Generating random Euler angles (Roll, Pitch, Yaw) for Attitude in radians
            # Roll, Pitch, Yaw in an L1 orbit: small deviations typically within +/- 0.5 degrees
            # Converting these to radians for precise control

            roll_values = np.deg2rad(np.random.uniform(-0.5, 0.5, size=1))
            pitch_values = np.deg2rad(np.random.uniform(-0.5, 0.5, size=1))
            yaw_values = np.deg2rad(np.random.uniform(-0.5, 0.5, size=1))

            # Generating random position and velocity values for Ephemeris in an L1 orbit
            # Position in km around L1: minor deviations from the L1 point in the Sun-Earth system
            # Velocity in km/s around L1: typical small velocities in each direction

            position_x = np.random.uniform(1.5e6, 1.52e6, size=1)  # in km, slight offsets from exact L1
            position_y = np.random.uniform(-2000, 2000, size=1)     # minor deviation in km
            position_z = np.random.uniform(-2000, 2000, size=1)     # minor deviation in km

            velocity_x = np.random.uniform(-1.0, 1.0, size=1)       # in km/s, small orbital velocity
            velocity_y = np.random.uniform(-1.0, 1.0, size=1)       # in km/s, minor variations
            velocity_z = np.random.uniform(-0.5, 0.5, size=1)       # in km/s, typically smaller z-velocity

            # Adding these values under SpiceData
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/RightAscension", ra_values)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Declination", dec_values)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Attitude/Roll", roll_values)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Attitude/Pitch", pitch_values)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Attitude/Yaw", yaw_values)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Ephemeris/PositionX", position_x)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Ephemeris/PositionY", position_y)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Ephemeris/PositionZ", position_z)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Ephemeris/VelocityX", velocity_x)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Ephemeris/VelocityY", velocity_y)
            create_dataset_if_not_exists(h, f"/{k[0]}/SpiceData/Ephemeris/VelocityZ", velocity_z)

            # print(np.array(v))
            # h.create_dataset(k, data=np.array(v, dtype=np.int8))
            # print(f"Time = {self.header}")
            if(f"/{k[0]}/Metadata/Epoch" not in h):
                h.create_dataset(f"/{k[0]}/Metadata/Epoch", data = self.header[(int(k[0]), 'Timestamp')])
            # Apply transformation based on k[1] (waveform name)
            if k[1] in conversion_factors:
                transformed_data = np.array(v) * conversion_factors[k[1]]
                h.create_dataset(f"/{k[0]}/{k[1]}", data=transformed_data)

                if k[1] in ['TOF H']:  # Calculate mass scale from TOF H arrays (best quality of the 3 gain stages)
                    stretch, shift, mass_scale = time2mass(transformed_data, self.hstime)
                    create_dataset_if_not_exists(h, f"/{k[0]}/Mass", data=np.array(mass_scale))
                    peaks, _ = find_peaks(transformed_data, prominence=.01)
                    print(f"peaks = {peaks}")

                    # Plot the data
                    plt.figure(figsize=(10, 6))
                    plt.plot(self.hstime, transformed_data, label="Transformed Data", color="blue")
                    plt.scatter(self.hstime[peaks], transformed_data[peaks], color="red", label="Peaks", marker="x", s=100)

                    # Add labels and legend
                    plt.xlabel("Time (hstime)")
                    plt.ylabel("Transformed Data")
                    plt.title("Transformed Data with Detected Peaks")
                    plt.legend()
                    plt.grid(True)

                    # Show the plot
                    # plt.show()

                    kappa = np.mean([mass_scale[peak]-np.round(mass_scale[peak], 1) for peak in peaks])
                    print(f"Kappa = {kappa}")
                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/kappa", data=np.array([kappa]))

                    # Find indices where Time (High Sampling) is between -7 and -5
                    mask = np.logical_and(self.hstime >= -7, self.hstime <= -5)

                    # Append the corresponding TOF H values to the baseline array
                    TOFmax = max(transformed_data) - np.mean(transformed_data[mask])
                    TOFsigma = np.std(transformed_data[mask])
                    SNR = TOFmax/TOFsigma
                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/SNR", data=np.array([SNR]))

                    fit_results = []
                    for peak in peaks:
                        print("Analyzing peak {peak}")
                        # Take a slice of 5 samples on either side of the peak
                        start = max(0, peak - 5)
                        end = min(len(transformed_data), peak + 6)  # 6 because slicing is exclusive of the end
                        x_slice = self.hstime[start:end]
                        y_slice = transformed_data[start:end]
                        print(f"x_slice = {x_slice}, y_slice = {y_slice}")
                        
                        # Fit the EMG to the slice
                        print(f"Calculating mass fit for peak {mass_scale[peak]}")
                        param, param_cov, sig_amp, fitted_curve = FitEMG(x_slice, y_slice)
                        if param is not None:
                            area = calculate_area_under_emg(x_slice, param)
                            print(f"Area under the EMG fit for peak {mass_scale[peak]}: {area}")
                            fit_results.append((param, sig_amp, x_slice, fitted_curve, area))
                        else:
                            fit_results.append(None)

                        # Overlay EMG fits
                        for result in fit_results:
                            if result is not None:
                                print(f"result = {result}")
                                param, sig_amp, x_slice, fitted_curve, area = result
                                if(type(param) is not int):
                                    param = param.tolist()
                                    print(f"Param = {param}")
                                    print(f"writing fit results for mass {param[0]}")
                                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/{k[1]}/Masses/{param[0]}FitParams", data=np.array(param))
                                    print(f"Area under fit for mass {param[0]}: {area}")
                                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/{k[1]}/Masses/{param[0]}AreaUnderFit", data=np.array([area]))



                if k[1] in ['Target L', 'Target H', 'Ion Grid']:  # Fit target and ion grid signals
                    param, param_cov, sig_amp, chi_s, r_chi_s, ion_result, ion_time = FitTargetSignal(self.lstime, v)
                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/{k[1]}FitParams", data=np.array(param))
                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/{k[1]}MassEstimate", data=sig_amp)
                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/{k[1]}ImpactCharge", data=sig_amp)
                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/{k[1]}ChiSquared", data=chi_s)
                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/{k[1]}ReducedChiSquared", r_chi_s)
                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/{k[1]}FitResult", ion_result)
                    create_dataset_if_not_exists(h, f"/{k[0]}/Analysis/{k[1]}FitTime", ion_time)


            else:
                h.create_dataset(f"/{k[0]}/{k[1]}", data=np.array(v))

            if(k[1]=='TOF L'):
                h.create_dataset(f"/{k[0]}/Time (high sampling)", data=self.hstime)
            if(k[1]=='Ion Grid'):
                h.create_dataset(f"/{k[0]}/Time (low sampling)", data=self.lstime)
        os.chdir('../')
        # h.create_dataset("Time since ")

# ||
# ||
# || Parse the high sampling rate data, this
# || should be 10-bit blocks
def parse_hs_waveform(waveform_raw: str):
    """Parse a binary string representing a high gain waveform"""
    w = bitstring.ConstBitStream(bin=waveform_raw)
    ints = []
    while w.pos < len(w):
        w.read('pad:2')  # skip 2
        ints += w.readlist(['uint:10']*3)
    print(len(ints))
    return ints[:-4]

# ||
# ||
# || Parse the low sampling rate data, this
# || should be 12-bit blocks
def parse_ls_waveform(waveform_raw: str):
    """Parse a binary string representing a low gain waveform"""
    w = bitstring.ConstBitStream(bin=waveform_raw)
    ints = []
    while w.pos < len(w):
        w.read('pad:8')  # skip 2
        ints += w.readlist(['uint:12']*2)
    print(len(ints))
    return ints

# ||
# ||
# || Use the SciType flag to determine the sampling rate of
# || the data we are trying to parse
def parse_waveform_data(waveform: str, scitype: int):
    """Parse the binary string that represents a waveform"""
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
    aparser = argparse.ArgumentParser()
    aparser.add_argument("--file", "-f", type=str, required=True)
    args = aparser.parse_args()

    packets = IDEXEvent(args.file)
    # print(packets.data.keys())
    try:
        packets.plot_all_data(packets.data, args.file)
    except Exception as e:
        print(e)
    packets.write_to_hdf5(packets.data, args.file+'.h5')
    # write_to_cdf(packets)

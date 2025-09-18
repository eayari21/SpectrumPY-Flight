#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================
A series of objects to store dust accelerator data effectively.
__author__      = Ethan Ayari,
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
===============================================================
"""

import os
import datetime
import time
import h5py
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from readTrc import Trc
from csv import writer
from scipy.optimize import curve_fit
from scipy.signal import detrend, butter, filtfilt
from sqlalchemy import create_engine
from lmfit import Model
from concurrent.futures import ThreadPoolExecutor


# %% GENERAL LINEAR FUNCTION DEFINITION
# || Used to subtract an overall linear background of noise
# || from the "baseline" of our signal

def LinearFit(Time, a, b):
    return a*Time + b

# %% GENERAL SINUSOIDAL FUNCTION DEFINITION
# || Used to subtract an overall beat pattern of noise
# || from the "baseline" of our signal

def SineFit(Time, c, d, e):
    return c*np.sin(d*Time + e)

# %% DUST ACCELERATOR PICKUP TUBE FUNCTION DEFINITION
# || Mihály's fit, derived from the charged particle dynamics
# || 

def QDFit(Time, t0, q, v):
    # |\Fit pickup tube (PUT) signals using charge (Q) and speed (v), and time of entering it t0
    # ||Accelerator U_acc = 2.2MV
    # ||PUT gain = .05 +/- something pc/V /*TO-DO*/ 
    # ||PUT length = 20 cm 
    # ||PUT diameter = 0.9 cm 
    # ||PUT RC time constant tau=0.25e-6
    # ||Particle mass = 2*Q*U_acc/v^2  to be calculated from best fit
    # |/Scope sampling rate dt= 0.8e-9

    # Input parameters and set up empty signal
    signal = np.zeros(len(Time))
    t0, q, v = t0, q, v  # To be fitted
    PUT_l = 20   # Length of the entire tube in cm
    PUT_d = 0.9  # Diameter of the PUT in cm
    tau = 0.16e-3  # Decay time, different than above??

    # Delineate first diameter vs. main body flight times
    dt1    = PUT_d/v  # Travel time for 1 diameter of the PUT
    dt2   = (PUT_l - 2*PUT_d)/v  # Travel time in the main body of PUT
    # print(f"Travel time for 1 diameter of the PUT = {dt1}s, and Travel time in the main body of PUT = {dt2}s")
    # print(f"")
    # print(f"Time goes from {min(Time)}s to {max(Time)}s")

    # Take a slice where the particle is in the first diameter

    index  = np.where((Time > t0) & (Time < t0+dt1))  # Time spend in 1st diameter
    index = [a for a in index]
    # print(f"Index 1 = {index}")

    signal[index] = q/dt1*(Time[index]-t0)  # signal from the first 1 diameter

    # Take a slice where the particle is in the pickup tube main body
    index = np.where((Time > t0+dt1) & (Time < t0+dt1+dt2))
    index = [a for a in index]
    # print(f"Index 2 = {index}")

    signal[index] = q*np.exp(-(Time[index]-(t0+dt1))/tau)  # signal from the main body of PUT

    ii = len(index)
    # print(ii)
    # print(f"Last index = {max(index[ii-1])}")
    qq = signal[max(index[ii-1])]  # last value in the main body of PUT
    # print(f"qq = {qq}")

    # Take a slice where the particle is in the last diameter
    index = np.where((Time > t0+dt1+dt2) & (Time < t0+2*dt1+dt2))  # signal from the last 1 diameter
    index = [a for a in index]

    # print(f"Index 3 = {index}")
    # print(index)
    signal[index] = qq-q / dt1*(Time[index] - (t0+dt1+dt2))

    ii = len(index)
    qq = signal[max(index[ii-1])]  # value after leaving the tube

    # Take a slice of everything after the particle leaves the pickup tube
    index = np.where(Time > t0+2*dt1+dt2)

    signal[index] = qq*np.exp(-1*(Time[index]-(t0+2*dt1+dt2))/tau) 

    return signal


# %%IDEX ION GRID FUNCTION DEFINITON
def IDEXIonGrid(x, P0, P1, P4, P5, P6):
    return P1 + np.heaviside(x-P0, 0) * ( P4 * (1.0 - np.exp(-(x-P0)/P5)) * np.exp( -(x-P0)/P6))


# %%OBJECT TO STORE WAVEFORMS, METADATA, ETC. FOR A GIVEN EXPERIMENT
class ImpactBook():
    def __init__(self, ChannelNames, trcdir="", ExperimentName = "", ChannelUnits=[]):

        self.ExperimentName = os.path.basename(os.path.normpath(trcdir))

        times, amps, metas = generalReadTRC(trcdir)
        self.NumEvents = len(times)
        ImpactEventList = []

        IonFolder = os.path.join(os.getcwd(), f"IonGridFits/{ExperimentName}_IonGridFits")
        os.makedirs(IonFolder)
        TargetFolder = os.path.join(os.getcwd(), f"TargetFits/{ExperimentName}_TargetFits")
        os.makedirs(TargetFolder)
        QDFolder = os.path.join(os.getcwd(), f"QDFits/{ExperimentName}_QDFits")
        os.makedirs(QDFolder)

        h5file = "./HDF5/" + f"{ExperimentName}.h5"

        for tracenum in range(self.NumEvents - 1):
            ImpactEventList.append(ImpactEvent(tracenum, ChannelNames, times[tracenum], amps[tracenum], metas[tracenum], ExperimentName, IonFolder, TargetFolder, QDFolder, h5file, ChannelUnits))


# %%OBJECT TO STORE WAVEFORMS FOR A SPECIFIC TRACE
class ImpactEvent():
    def __init__(self, Number, ChannelNames, times, amps, metas, ExperimentName, IonFolder, TargetFolder, QDFolder, h5file, ChannelUnits=[]):
        # print(f"There are {len(amps)} channels")
        # print(f"Times shape is {np.asarray(times).shape}, Metas shape is {np.asarray(metas).shape}, Amps shape is {np.asarray(amps).shape}")
        self.h5file = h5file
        self.TraceNumber = Number
        # self.ChannelUnits = ChannelUnits
        self.MetaData = metas
        # self.DataDict = dict(zip(ChannelNames, amps))
        # self.DataDict.update({'Time': times[0]})
        self.IonFolder = IonFolder
        self.TargetFolder = TargetFolder
        self.QDFolder = QDFolder

        # for k, v in self.DataDict.items():
        #    print(k, " : ", len(v), " items")

        # print(f"Shape of amps = {np.array(amps[0], dtype=object).shape}. {len(times)}")

        # Convert to a DataFrame
        self.Waveforms = pd.DataFrame(data=np.array(amps).T, dtype=float)  # Transpose the array to make rows correspond to N samples
        self.Waveforms.columns = ChannelNames       # Set column names
        self.Waveforms["Time"] = times[0]
        # print(self.Waveforms)


        with h5py.File(h5file, mode="a") as h:
            for col in self.Waveforms.columns:
                h.create_dataset(f"/{Number}/{col}", data=np.array(self.Waveforms[col]))

        """
        ||=================================================||
        Time to gather our information:
        ||=================================================||
        Remember, what's important for each IDEX
        dataset is:
        1.  ||Particle charge   (QD)
        2.  ||Particle velocity (QD Fit)
        3.  ||Trigger time (QD)
        4.  ||QD mass
        5.  ||Particle mass (Accelerator)
        6.  ||Particle velocity (Accelerator)
        7.  ||Particle charge (Accelerator)

        ————------------————^^^ dust ^^^^---------------------

        8.  ||Target rise time
        9.  ||Target amplitude (Peak-Baseline)
        10. ||Impact time estimate (5 % ? of the max charge)
        11. ||Ion grid rise time
        12. ||Ion grid amplitude (Peak)
        13. ||Ion grid signal start time (5% of peak?)
        ||=================================================||
        """

        try:
            self.IonGridFit, self.IonChiSquared, self.IonRedChiSquared, self.IonGridFitResult, self.IonFitTime = self.FitIonGridSignal(self.Waveforms["Ion Grid"])
            self.TargetFit, self.TargetChiSquared, self.TargetRedChiSquared, self.TargetFitResult, self.TargetFitTime, self.TargetFitData = self.FitTargetSignal(self.Waveforms["Target H"])

            print("Fitting accelerator data")
            self.QDFit, self.QDChiSquared, self.QDRedChiSquared, self.QDFitResult, self.QDFitTime = self.FitQDSignal(self.Waveforms["QD High"])
            self.FittedTrigger = self.QDFit[0]
            self.FittedCharge = self.QDFit[1]
            self.FittedVelocity = self.QDFit[2]
            self.SQLData = self.RetrieveSQL(self.QDFit[2])

            self.IonGridAmplitude = 0.519602*self.IonGridFit[2]
            self.IonGridRiseTime = self.IonGridFit[3]
            self.IonGridTrigger = self.IonGridFit[0]

            self.TargetAmplitude = 2.735834*self.TargetFit[2]
            self.TargetRiseTime = self.TargetFit[3]
            self.TargetTrigger = self.TargetFit[0]

            self.AssignedVelocity = self.SQLData["Velocity (km/s)"].values[0]
            self.AssignedMass = self.SQLData["Mass (kg)"].values[0]
            self.AssignedCharge = self.SQLData["Charge (C)"].values[0]*1e16
            self.AssignedRadius = self.SQLData["Radius (m)"].values[0]



            # Labels and values
            data_dict = {

                "Trace Number": Number,
                "Target Rise Time (s)": self.TargetRiseTime,
                "Target Amplitude (pC)": self.TargetAmplitude,
                "Ion Grid Rise Time (s)": self.IonGridRiseTime,
                "Ion Grid Amplitude (pC)": self.IonGridAmplitude,
                "Ion Grid Trigger (s)": self.IonGridTrigger,
                "Target Trigger (s)": self.TargetTrigger,
                "Ion Chi Squared": self.IonChiSquared,
                "Ion Reduced Chi Squared": self.IonRedChiSquared,
                "Ion Grid Fit Result": self.IonGridFitResult,
                "Ion Grid Fit Time": self.IonFitTime,

                "Target Chi Squared": self.TargetChiSquared,
                "Target Reduced Chi Squared": self.TargetRedChiSquared,
                "Target Fit Result": self.TargetFitResult,
                "Target Fit Time": self.TargetFitTime,
                "Target Fit Data": self.TargetFitData,

                "Particle Charge (pC, Accelerator)" : self.AssignedCharge, 
                "Particle Velocity (kmps, Accelerator)" : self.AssignedVelocity, 
                "Particle Mass (kg, Accelerator)" :  self.AssignedMass, 
                "Particle Radius (m, Accelerator)" : self.AssignedRadius, 
                "Particle Charge (pC, QD Fit)" : self.FittedCharge, 
                "Particle Velocity (kmps, QD Fit)" : self.FittedVelocity,
                "QD Chi Squared": self.QDChiSquared,
                "QD Reduced Chi Squared": self.QDRedChiSquared,
                "QD Fit Result": self.QDFitResult,
                "QD Fit Time": self.QDFitTime

            }

            # Write to HDF5
            with h5py.File(self.h5file, "a") as h:  # "a" mode to append data
                # Write each label-value pair as a separate dataset
                for label, value in data_dict.items():
                    # If the value is a scalar, write it directly
                    if np.isscalar(value):
                        print(f"Writing {label}")
                        h.create_dataset(f"{Number}/Analysis/{label}", data=np.array(value))
                    # If the value is an array, write it as a dataset
                    else:
                        print(f"Writing {label}")
                        h.create_dataset(f"{Number}/Analysis/{label}", data=np.array(value))

        except Exception as e:
            print(f"Calculations failed for shot number {Number} due to {e}")
            pass

        # try:

        #     if (Number==0):

        #         # with open(f"{ExperimentName}_data.csv",'w') as f:

        #         with open(os.path.join(os.getcwd(), f"{ExperimentName}.csv"),'w') as f:
        #             csv_writer = writer(f)
        #             csv_writer.writerow(["Trace Number",  "Target Rise Time (s)", "Target Amplutude (pC)", "Ion Grid Rise Time (s)", "Ion Grid Amplitude (pC)", "Ion Grid Trigger (s)", "Target Trigger (s)", "Ion Fit Quality", "Target Fit Quality"])
        #             csv_writer.writerow([Number, self.TargetRiseTime, self.TargetAmplitude, self.IonGridRiseTime, self.IonGridAmplitude, self.IonGridTrigger, self.TargetTrigger, np.linalg.det(self.IonChiSquared), np.linalg.det(self.TargetChiSquared)])
        #             # csv_writer.writerow(["Target Rise Time (s)", "Target Amplutude (pC)", "Ion Grid Rise Time (s)", "Ion Grid Amplitude (pC)", "Ion Grid Trigger (s)", "Target Trigger (s)", "Ion Fit Quality", "Target Fit Quality"])
        #             # csv_writer.writerow([self.TargetRiseTime, self.TargetAmplitude, self.IonGridRiseTime, self.IonGridAmplitude, self.IonGridTrigger, self.TargetTrigger, self.IonChiSquared, self.TargetChiSquared])
        #     else:
        #         # with open(f"{ExperimentName}_data.csv",'a') as f:

        #         with open(os.path.join(os.getcwd(), f"{ExperimentName}.csv"),'a') as f:
        #             csv_writer = writer(f)
        #             csv_writer.writerow([Number, self.AssignedCharge, self.AssignedVelocity, self.AssignedMass, self.AssignedRadius, self.FittedCharge, self.FittedVelocity, self.FittedTrigger, self.TargetRiseTime, self.TargetAmplitude, self.IonGridRiseTime, self.IonGridAmplitude, self.IonGridTrigger, self.TargetTrigger, np.linalg.det(self.IonChiSquared), np.linalg.det(self.TargetChiSquared), np.linalg.det(self.QDQual)])
        #             # csv_writer.writerow([self.TargetRiseTime, self.TargetAmplitude, self.IonGridRiseTime, self.IonGridAmplitude, self.IonGridTrigger, self.TargetTrigger, self.IonChiSquared, self.TargetChiSquared])

        # except:
        #     print(f"Csv dump failed for shot number {Number}")
        #     pass

    # %%Target Signal Fitting Routine %% #

    # || Very noisy due to "microphonics", so we will:
    # || 1) Remove a linear baseline (y = a*x + b), and 
    # || 2) Remove a sinusoidal background (y = c*sin(d*x + e)

    def FitTargetSignal(self, targetAmp, Plot=True):
        self.x, self.y = self.Waveforms["Time"], targetAmp
        # Use the "idx" variable to delineate where our target signal should be
        idx = self.Waveforms.index[(self.Waveforms["Time"] >= -5e-5) & (self.Waveforms["Time"] <= .0004)]

        # || Select only raw noise (where we know the signal is not)

        basedex = self.Waveforms.index[self.Waveforms["Time"] <= -1e-5]
        baselineraw = self.y[basedex]
        baselinedomain = self.Waveforms["Time"][basedex]

        # || Remove Linear Background

        try:
            # slopeguess = (baselineraw[len(baselineraw)-1] - baselineraw[0]) / (baselinedomain[len(baselinedomain-1)] - baselinedomain[0])
            slopeguess = 0
            linparam, lin_cov = curve_fit(LinearFit, baselinedomain, baselineraw, p0=[slopeguess,0], maxfev=100_000)
            linearbase = LinearFit(self.Waveforms["Time"], linparam[0], linparam[1])
            # self.y -= linearbase
            self.y = detrend(self.y)
        except:
            print(f"Linear background not found for trace number {self.TraceNumber}.")
            linearbase = None


        # || Remove Sine Wave Background

        try:
            baselinedelined = self.y[basedex]
            sineparam, sine_cov = curve_fit(SineFit, baselinedomain, baselinedelined, p0=[max(baselinedelined), 7000, 45], maxfev=100_000)
            sinebase = SineFit(self.Waveforms["Time"], sineparam[0], sineparam[1], sineparam[2])
            self.y -= sinebase
            self.y = butter_lowpass_filter(self.y, self.Waveforms["Time"])


        except:
            print(f"Sinusoidal background not found for trace number {self.TraceNumber}.")
            sinebase=None


        # print(f"idx = {idx}")
        self.x = self.Waveforms["Time"][idx]
        self.y = self.y[idx]
        self.traceNum = self.TraceNumber
        self.pre = -3.0e-5  # Before image charge
        self.yBaseline = self.y.where(self.x < self.pre)
        self.yImage = self.y[(self.x.values >= self.pre) & (self.x.values < 0.0)]
        self.ionError = self.yBaseline.std()
        print(f"Ion error = {self.ionError}")
        self.ionErrorVector = pd.DataFrame([np.nan] * len(self.yBaseline))
        self.ionMean = self.yBaseline.mean()
        self.yBaseline -= self.ionMean
        self.parameters = []

        ionTime = np.array(self.x.apply(lambda x: float(x)))
        ionAmp = np.array(self.y.apply(lambda x: float(x)))

        print(f"Calculating Target Fit number {self.TraceNumber}")
        # || Initial Guess for the parameters of the ion grid signal

        t0 = 0.0                         # P[0] time of impact
        c = 0.                           # P[1] Constant offset
        # b = np.abs(min(self.yImage))     # P[2] Image amplitude
        # b = .01

        # s = 4.e-6                        # P[3] Image pulse width
        A  = max(ionAmp)  # np.abs(min(ionAmp) - max(ionAmp))     # P[4] amplitude (v)
        # A = .05

        t1 = 3.71e-5                      # P[5] rise  time (s)
        t2 = 3.71e-4                      # P[6] discharge time (s)
        c = 0

        P0, P1, P4, P5, P6 = t0, c, A, t1, t2

        model = Model(IDEXIonGrid)
        param = model.make_params(P0=P0, P1=P1, P4=P4, P5=P5, P6=P6)
        param['P5'].set(value=P5, min=1e-6, max=1e-3)
        # Fit the data
        result = model.fit(ionAmp, param, x=ionTime)

        # Calculate chi-squared values
        chi_squared = result.chisqr
        reduced_chi_squared = result.redchi

        self.result = result.best_fit
        param['P4'].value = max(self.result) - self.yBaseline.mean()


        if Plot:

            plt.style.use('seaborn-pastel')

            # || First, plot the original signal and the noise that has been removed from it.
            plt.cla()
            plt.clf()

            plt.ylabel("Voltage (V)", fontsize=15)
            plt.xlabel("Time (s)", fontsize=15)
            plt.title(f"Original Target Signal number {self.TraceNumber}", fontweight='bold', fontsize=20)
            plt.plot(self.Waveforms["Time"][idx], targetAmp[idx], label="Original target signal (V)")
            plt.plot(self.Waveforms["Time"][idx], self.y, label="Filtered Signal (V)")


            if linearbase is not None:
                plt.plot(self.Waveforms["Time"], linearbase, label="Linear noise background (V)")

            if sinebase is not None:
                plt.plot(self.Waveforms["Time"], sinebase, label="Sine wave noise background (V)")

            plt.legend(loc='best')
            plt.savefig(os.path.join(self.TargetFolder, f"TargetSignal{self.TraceNumber}.png"))
    
            plt.cla()
            plt.clf()

            plt.ylabel("Voltage (V)", fontsize=15)
            plt.xlabel("Time (s)", fontsize=15)

            plt.plot(self.x, self.y, label="Target (V)")
            plt.plot(self.x, self.result, label = "Fit (V)")

            plt.figtext(0.5, 0.2,
                r"Amplitude: %.2e C, $\tau_{Rise}$: %.2e s, $\tau_{Discharge}$: %.2e s"%(2.735834054*max(self.result), param['P5'], 3.71e-4),
                horizontalalignment ="center", 
                verticalalignment ="center", 
                wrap = True, fontsize = 14, 
                color ="orange")
            plt.title(f"Target Signal Fit number {self.TraceNumber}", fontweight='bold', fontsize=20)
            plt.legend(loc='best')
            plt.savefig(os.path.join(self.TargetFolder, f"TargetSignalFit{self.TraceNumber}.png"))

        return([param[key].value for key in param.keys()], chi_squared, reduced_chi_squared, self.result, self.x, self.y)

    # %% Ion Grid Signal Fit
    # || Identical to the Target signal fit with less noise reduction

    def FitIonGridSignal(self, ionGridAmp, Plot=True):
        self.x, self.y = self.Waveforms["Time"], ionGridAmp
        idx = self.Waveforms.index[(self.Waveforms["Time"] >= -5e-5) & (self.Waveforms["Time"] <= .001)]
        basedex = self.Waveforms.index[self.Waveforms["Time"] <= -2e-5]
        basemean = ionGridAmp[basedex].mean()

        self.y -= basemean
        # print(f"idx = {idx}")
        self.x = self.Waveforms["Time"][idx]
        self.y = ionGridAmp[idx]
        self.traceNum = self.TraceNumber
        self.pre = -3.0e-5  # Before image charge
        self.yBaseline = self.y.where(self.x < self.pre)
        self.yImage = self.y[(self.x.values >= self.pre) & (self.x.values < 0.0)]
        self.ionError = self.yBaseline.std()
        print(f"Ion error = {self.ionError}")
        self.ionErrorVector = pd.DataFrame([np.nan] * len(self.yBaseline))
        self.ionMean = self.yBaseline.mean()
        self.yBaseline -= self.ionMean
        self.parameters = []

        ionTime = np.array(self.x.apply(lambda x: float(x)))
        ionAmp = np.array(self.y.apply(lambda x: float(x)))

        # || Apply a low pass filter to improve our fits

        originalAmp = ionAmp
        self.y = butter_lowpass_filter(ionAmp, ionTime)

        print(f"Calculating Ion Fit number {self.TraceNumber}")
        # %% Initial Guess for the parameters of the ion grid signal
        t0 = 0.0                         # P[0] time of impact
        c = 0.                           # P[1] Constant offset
        # b = np.abs(min(self.yImage))     # P[2] Image amplitude
        # b = .01
        # s = 4.e-6                        # P[3] Image pulse width
        A  = max(ionGridAmp) # np.abs(min(self.y) - max(self.y))     # P[4] amplitude (v)
        # A = .05
        t1 = 3.71e-5                      # P[5] rise  time (s)
        t2 = 3.71e-4                      # P[6] discharge time (s)
        c = 0

        P0, P1, P4, P5, P6 = t0, c, A, t1, t2

        model = Model(IDEXIonGrid)
        param = model.make_params(P0=P0, P1=P1, P4=P4, P5=P5, P6=P6)
        # Fit the data
        result = model.fit(ionAmp, param, x=ionTime)

        # Calculate chi-squared values
        chi_squared = result.chisqr
        reduced_chi_squared = result.redchi

        self.result = result.best_fit
        param['P4'].value = max(self.result) - self.yBaseline.mean()




        if Plot:

            plt.style.use('seaborn-pastel')
            plt.cla()
            plt.clf()
            plt.xlabel("Time (s)", fontsize=15)
            plt.ylabel("Voltage (V)", fontsize=15)

            plt.plot(self.x, originalAmp, lw = .5, label="Ion Grid (V)")
            plt.plot(self.x, self.y, label="Filtered Signal (V)")
            plt.plot(self.x, self.result, label = "Fit (V)")

            plt.axhline(max(self.result), color='black', linestyle='-.')  # Black dash-dot line
            plt.axhline(min(self.result), color='black', linestyle='--')  # Black dashed line

            plt.figtext(0.5, 0.2,
                r"Amplitude: %.2e C, $\tau_{Rise}$: %.2e s, $\tau_{Discharge}$: %.2e s"%(.519602*max(self.result), param['P5'], 3.71e-4),
                horizontalalignment ="center", 
                verticalalignment ="center", 
                wrap = True, fontsize = 14, 
                color ="orange")

            plt.title(f"Ion Grid Signal Fit number {self.TraceNumber}", fontweight='bold', fontsize=20)
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(os.path.join(self.IonFolder, f"IonGrid{self.TraceNumber}.png"))

        return([param[key].value for key in param.keys()], chi_squared, reduced_chi_squared, self.result, self.x)
        
    # %%
    # || Employ Mihály's function (above) to assign a tentative velocity

    def FitQDSignal(self,QDAmplitude, Plot=True):
        print(f"Calculating QD Fit number {self.TraceNumber}")

        idx = self.Waveforms.index[self.Waveforms["Time"] <= 0.0]
        domain = self.Waveforms["Time"][idx]
        # print(f"Timestep = {domain[1] - domain[0]}")

        range = QDAmplitude[idx]

        range2 = butter_lowpass_filter(range.to_numpy(), domain.to_numpy())

        # print("Running Calcs...")

        # Assume the trigger time is where the QD waveform hits its minimum
        mindex = range.idxmin().astype(int)
        # print(mindex)
        range2 -= range2[0:mindex-200_000].mean()
        range -= range2[0:mindex-200_000].mean()
        # t_trig = domain[mindex-15_000]
        first_qdl = range2[0:mindex]
        second_qdl = range2[mindex:len(range2)]
        midpoint = 2*min(range2)/3.0
        # print(f"Midpoint = {midpoint}")
        t_a = domain[find_nearest(first_qdl, midpoint)]
        t_b = domain[find_nearest(second_qdl, midpoint)+mindex]
        t_trig = domain[find_nearest(first_qdl, midpoint)]
        # print(f"t_a: {t_a} t_b: {t_b}")
        # range = pd.Series((min(range)/min(range2))*range2)

        # ||PU Tube length - both diameters used for distance

        v_init = 18.2/(t_b - t_a)
        model = Model(QDFit)
        params = model.make_params(t0=t_trig, q=min(range2), v=v_init)

        # Fit the data
        result = model.fit(range2, params, Time=domain.to_numpy())

        # Extract fit results
        t0_fit = result.params['t0'].value
        q_fit = result.params['q'].value
        v_fit = result.params['v'].value

        # Calculate chi-squared values
        chi_squared = result.chisqr
        reduced_chi_squared = result.redchi

        if Plot:
            plt.style.use('seaborn-pastel')
            plt.cla()
            plt.clf()

            plt.ylabel("Voltage (V)", fontsize=15)
            plt.plot(domain, range, label="Original Signal (V)")
            plt.plot(domain, range2, label="Filtered Signal (V)")
            plt.plot(domain, result.best_fit, label="Fit (V)")
            plt.figtext(
                0.5, 0.2,
                r"$t_{0}$: %.2e s, q: %.2e pC, v: %.2e $\frac{km}{s}$" % (t0_fit, 0.05 * q_fit, v_fit / 100_000),
                horizontalalignment="center",
                verticalalignment="center",
                wrap=True,
                fontsize=14,
                color="orange",
            )
            plt.title(f"QD Pickup Tube Fit number {self.TraceNumber}", fontweight='bold', fontsize=20)
            plt.legend(loc='best')
            plt.savefig(os.path.join(self.QDFolder, f"QDFit{self.TraceNumber}.png"))


        return [t0_fit, 0.05 * q_fit, v_fit / 100_000], chi_squared, reduced_chi_squared, result.best_fit, domain

    # %%
    # || Fetch all relevant metadata from the accelerator's SQL server

    def RetrieveSQL(self, velocity):
        timeArr = self.MetaData[4]['TRIGGER_TIME']
        date_time = datetime.datetime(*timeArr)
        date_time = date_time.timestamp()

        timestamp = 1000*int(date_time)
        db_connection_str = 'mysql+pymysql://admin:1jlwbXqCVkNdt91KCijx@impactlab.cbhrpdddmhcs.us-east-1.rds.amazonaws.com/CCLDAS_PRODUCTION'
        db_connection = create_engine(db_connection_str)

        time = str(timestamp)
        # velocity = str(velocity)

        self.SQL_df = pd.read_sql("select estimate_quality, integer_timestamp, velocity, mass, charge, radius from dust_event where mass != -1 and integer_timestamp < {} and velocity > {} and velocity < {} order by integer_timestamp desc limit 1;".format(timestamp, str(velocity*1000-100), str(velocity*1000+100)), con=db_connection)

        self.SQL_df.columns = ["Estimate Quality", "Time", "Velocity (km/s)", "Mass (kg)", "Charge (C)", "Radius (m)"] 
        self.SQL_df["Velocity (km/s)"] /= 1000.0
        self.SQL_df = self.SQL_df.sort_values(by=["Time"], ascending=True)

        self.SQL_df["Time"] = pd.to_datetime(self.SQL_df['Time'], unit='ms')
        self.SQL_df["Time"] = self.SQL_df["Time"].dt.tz_localize('utc').dt.tz_convert('US/Mountain')

        self.SQL_df['Trace Number'] = range(1, len(self.SQL_df) + 1)
        self.SQL_df = self.SQL_df.reset_index(drop=True)
        return self.SQL_df


# %%GENERAL TRACE READER
def generalReadTRC(dataDir, verbose=False):
    """A function to read in the binary trace files from the oscilloscopes. It
    is general for different instruments and scope configurations. The number
    of channels and shots is automatically detected and sorted.

        Args:
           dataDir (str): The folder containing the trace files. It doesn't
           matter if there are other files in that directory. If you are
           manually entering a path and encounter issues, try os.listdir() to
           ensure things look right.'

        Kwargs:
           verbose (bool): A boolean representing the user's choice of whether
           or not they want real-time print statements of all relevant
           quantities.

        Returns:
           times (list): A 2 dimensional list containing the time arrays for
           each shot.

           amps (list): A 2 dimensional list containing the amplitude arrays
           for each shot.

           metas (list): A 2 dimensional list containing the metadata arrays
           for each shot. **TODO** Figure out this format and make it a dict.

        Exceptions:
           FileNotFoundError: The cue to stop parsing the user-provuded
           directory for any more trace files.

            ZeroDivisionError: This will be reaised if there are no files
            detected at all.
           """
    # __author__= Ethan Ayari
    # Iterate through all files in the folder and attempt to sort them
    parse = True

    global times
    global amps
    global metas
    global nChannels
    global traceList
    times = []
    amps = []
    metas = []
    shotKey = "00000"
    nChannels = 0
    while parse:
        # Shot number must be given as a 5 digit number
        # Initialize empty scope object
        trc = Trc()
        tmpTime = []
        tmpAmp = []
        tmpMeta = []
        # Iterate through all files in passed directory
        parsetwo = False

        directoryList = sorted(os.scandir(dataDir), key=lambda e: e.name)
        for filename in directoryList:
            # print(str(filename))
            # Gather all channels for a given
            # shotKey (tracenumber)
            if(str(filename).find(shotKey) != -1):
                parsetwo = True
                if verbose:
                    print(str(filename.path))
                try:
                    t, y, meta = trc.open(str(filename.path))
                    # t, y, meta = trc.open(str(filename))
                    tmpTime.append(t)
                    tmpAmp.append(y)
                    tmpMeta.append(meta)
                    nChannels += 1
                except FileNotFoundError:
                    # We've found everything
                    parse = False
        if not parsetwo:
            parse = False

        # print(f"tmpTime has length {len(tmpTime)}")
        times.append(tmpTime)
        amps.append(tmpAmp)
        metas.append(tmpMeta)
        i = int(shotKey) + 1
        shotKey = f'{i:05d}'
    shotKey = int(shotKey)-1  # Get rid of extra shot count
    try:
        nChannels = int(int(nChannels) / int(shotKey))
    except ZeroDivisionError:
        print("No data was detected.")

    """
        |\1. Laser Photodiode
        ||2. TOF L
        ||3. TOF M
        ||4. TOF H
        ||5. Ion Grid (Flight = 1.499326 pc/V, GSE = 0.519602 pc/V)
        ||6. Target Low Range (Flight = 319.2277 pc/V, GSE = 0.232115 pc/V)
        |/7. Target High Range (Flight = 31936.65 pc/V, GSE = Nonexistent)
   
    for count, trace in enumerate(amps):
        print(f"Analyzing trace file {count-1}")
        amps[int(count-1)][4] *= 1.499326
        amps[int(count-1)][5] *= 319.2277
        amps[int(count-1)][6] *= 31936.65
    """

    if verbose:
        print(f"{nChannels} channels detected.")
        print(f"{shotKey} shots detected per each channel.")
    return times, amps, metas


# %%
def butter_lowpass_filter(data, time):
    # Filter requirements.
    T = time[1] - time[0]       # |\Sample Period (s)
    fs = (time[-1] - time[0])/T # ||sample rate, Hz
    cutoff = 100                # ||desired cutoff frequency of the filter, Hz
    nyq = 0.5 * fs              # ||Nyquist Frequency
    order = 2                   # ||sine wave can be approx represented as quadratic
    n = len(data)               # |/total number of samples
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

# %%
def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return idx


# %%EXECUTABLE CODE BELOW
if __name__ == "__main__":
    """
        |\1. Laser Photodiode
        ||2. TOF L
        ||3. TOF M
        ||4. TOF H
        ||5. Ion Grid (Flight = 1.499326 pc/V, GSE = 0.519602 pc/V)
        ||6. Target Low Range (Flight = 319.2277 pc/V, GSE = 0.232115 pc/V)
        |/7. Target High Range (Flight = 31936.65 pc/V, GSE = Nonexistent)
    """

    # ChannelNames = ["Laser Photodiode", "TOF L", "TOF M", "TOF H", "Ion Grid", "Target Low Range", "Target High Range"]
    # ChannelNames = ["QD Low", "QD High", "TOF L", "TOF M", "TOF H", "Ion Grid", "Target L", "Target H"]
    # ChannelNames = ["QD Low", "QD High", "TOF L", "TOF M", "TOF H", "Ion Grid", "Target H"]

    # ChannelNames = ["Ion Grid", "TOF H Gain", "TOF L Gain", "TOF M Gain", "Target CSA", "Velocity Grid CSA"]
    # # LaserPrimer = ImpactBook(ChannelNames, trcdir="/Users/impact/Dropbox/IDEX Pipeline/python/SpectrumPY/IDEX EM1 Laser Testing/Finding Focus")

    ChannelNames = ["QD Low", "QD High", "TOF L", "TOF M", "TOF H", "Ion Grid", "Target L", "Target H"]
    # ChannelNames = ["QD High", "Ion Grid", "Target H", "TOF L", "TOF M", "TOF H"]

    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Run ImpactBook with specified directories and experiment names.")
    parser.add_argument("--trcdir", type=str, required=True, help="Directory containing TRC files.")
    parser.add_argument("--experimentname", type=str, required=True, help="Name of the experiment.")

    # Parse the arguments
    args = parser.parse_args()
    trcdir = args.trcdir
    experimentname = args.experimentname

    # Call the ImpactBook function with the provided arguments
    ImpactBook(ChannelNames, trcdir=trcdir, ExperimentName=experimentname)

    # ImpactBook(ChannelNames, trcdir = "vbin_5-8", ExperimentName="Grid vel 5-8")



    # ImpactBook(ChannelNames, trcdir = "/Users/impact/Dropbox/IDEX Pipeline/3_Codes/3_SpectrumPY/0_SpectrumPY Dust Edition/0_IDEX_EM1_Dust_Testing/IDEX_Grid_Characterization/vbin_5-8", ExperimentName="Grid vel 5-8")
        # DustPrimer = ImpactBook(ChannelNames, trcdir = "~/Dropbox/IDEX Pipeline/3_Codes/3_SpectrumPY/0_SpectrumPY Dust Edition/0_IDEX_EM1_Dust_Testing/Dust_09.20.22/Target_Location_1", ExperimentName="Target Location 1")

    # os.chdir("0_IDEX_EM1_Dust_Testing/IDEX_Grid_Characterization")
    # print(os.listdir())

    # ImpactBook(ChannelNames, trcdir = "vbin_5-8", ExperimentName="Grid vel 5-8")
        # DustPrimer = ImpactBook(ChannelNames, trcdir = "~/Dropbox/IDEX Pipeline/3_Codes/3_SpectrumPY/0_SpectrumPY Dust Edition/0_IDEX_EM1_Dust_Testing/Dust_09.20.22/Target_Location_1", ExperimentName="Target Location 1")



    # ImpactBook(ChannelNames, trcdir = "vbin_12-20", ExperimentName="Grid vel 12-20")
        # DustPrimer = ImpactBook(ChannelNames, trcdir = "~/Dropbox/IDEX Pipeline/3_Codes/3_SpectrumPY/0_SpectrumPY Dust Edition/0_IDEX_EM1_Dust_Testing/Dust_09.20.22/Target_Location_1", ExperimentName="Target Location 1")




    # # try:
    # # ImpactBook(ChannelNames, trcdir = "/Users/impact/Dropbox/IDEX Pipeline/python/SpectrumPY/IDEX EM1 Laser Testing/Finding Focus", ExperimentName="Finding Focus Dataset")
    #     # ImpactBook(ChannelNames, trcdir = "/Users/impact/Dropbox/Mac (2)/Documents/IDEX_dust_09.28.22/Target_Location_Detector", ExperimentName="Target Location Detector v2")
    #     # DustPrimer = ImpactBook(ChannelNames, trcdir = "/Users/impact/Dropbox/IDEX Pipeline/python/SpectrumPY/Dust_09.20.22/Target_Location_1", ExperimentName="Target Location 1")
    # # except:
    # #    pass

    # try:
    #     # ImpactBook(ChannelNames, trcdir = "/Users/impact/Dropbox/IDEX Pipeline/python/SpectrumPY/IDEX EM1 Laser Testing/Finding Focus", ExperimentName="Finding Focus Dataset")

    #     # ImpactBook(ChannelNames, trcdir = "/Users/impact/Dropbox/Mac (2)/Documents/IDEX_dust_09.20.22/Target_Location_Detector", ExperimentName="Target Location Detector")
    #     DustPrimer = ImpactBook(ChannelNames, trcdir = "/Users/impact/Dropbox/IDEX Pipeline/3_Codes/2_Level 0->1/Python_Packet_Crushing/lasp_packets/tests/integration/test_xtce_based_parsing/sciData_2022_130_17_41_53.h5", ExperimentName="SciData53")
    # except:
    #     pass

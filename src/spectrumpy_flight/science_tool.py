#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

title = """

//  ================================================================================
//  ||                                                                            ||
//  ||              science_tool                                                  ||
//  ||              ------------------------------------------------------        ||
//  ||                           S C I E N C E   T O O L                          ||
//  ||              ------------------------------------------------------        ||
//  ||                                                                            ||
//  ||                __author__      = Ethan Ayari & Gavin Medley,               ||
//  ||                IMPACT/LASP, CU Boulder                                     ||
//  ||                                                                            ||
//  ||                For: IDEX Flight Model Integration and Test, L0-L1A         ||
//  ||                                                                            ||
//  ||                2023                                                        ||
//  ||                                                                            ||
//  ||                Dubbed by Sam Oberg, tested with Delphine Sherman           ||
//  ||                Works with Python 3.10.4                                    ||
//  ||                                                                            ||
//  ================================================================================


Science Tool: An all-in-one Python script to read in and 
visualize IDEX data from the GS-OASIS server.

"""


print(title)
# %%
# || Python libraries
import argparse
import os
import socket
from pathlib import Path
import bitstring
import h5py
import shutil
import datetime
# import signal
import matplotlib.pyplot as plt

from .plot_style import apply_plot_style
from .paths import default_hdf5_dir
from .spacecraft_clock import (
    SPACECRAFT_EPOCH,
    combine_coarse_fine_seconds,
    spacecraft_seconds_to_datetime,
)

apply_plot_style()
import numpy as np

from datetime import datetime, timedelta, timezone

# || LASP software
from lasp_packets import xtcedef  # Gavin Medley's xtce UML implementation
from lasp_packets import parser  # Gavin Medley's constant bitstream implementation
from .rice_decode import idex_rice_Decode
# import cdflib.cdfwrite as cdfwrite
# import cdflib.cdfread as cdfread

# ||
# ||
# || Generator object from LASP packets
# || to read in the data
class IDEXEvent:
    def __init__(self, filename: str):
        """Test parsing a real XTCE document"""
        # TODO: Change location of xml definition
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
                    # print(pkt.data)
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


                     # Account for HS trigger delay
                    self.TOFdelay = pkt.data['IDX__TXHDRSAMPDELAY'].raw_value >> 2  # First two bits are padding

                    # Mask to extract 10-bit values
                    mask = 0b1111111111

                    self.hgdelay = (self.TOFdelay) & mask
                    self.mgdelay = (self.TOFdelay >> 10) & mask
                    self.lgdelay = (self.TOFdelay >> 20) & mask
                    print(f"High gain delay = {self.hgdelay} samples.")
                    print(f"Mid gain delay = {self.mgdelay} samples.")
                    print(f"Low gain delay = {self.lgdelay} samples.")

                    if(pkt.data['IDX__TXHDRLSTRIGMODE'].derived_value!=0):  # If this was a LS (Target Low Gain) trigger
                        self.Triggerorigin = 'LS' 
                        print("Low sampling trigger mode enabled.")
                        # Check the first 25th-bit integer for a coincidence trigger
                        # coincidence = (trigmode >> 24) &  0b1
                        # if(coincidence==1):
                            # self.Triggermode = 'Coincidence'
                        # else:
                            # self.Triggermode = 'Threshold'
                    print("Packet trigger mode = ", pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value, pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value, pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value)
                    if(pkt.data['IDX__TXHDRLGTRIGMODE'].derived_value!=0):
                        print("Low gain TOF trigger mode enabled.")
                        self.Triggerorigin = 'LG'
                    if(pkt.data['IDX__TXHDRMGTRIGMODE'].derived_value!=0):
                        print("Mid gain TOF trigger mode enabled.")
                        self.Triggerorigin = 'MG'
                    if(pkt.data['IDX__TXHDRHGTRIGMODE'].derived_value!=0):
                        print("High gain trigger mode enabled.")
                        self.Triggerorigin = 'HG'

                    print(f"AID = {pkt.data['IDX__SCI0AID'].derived_value}")  # Instrument event number
                    print(f"Event number = {pkt.data['IDX__SCI0EVTNUM'].raw_value}")  # Event number out of how many events constitute the file
                    # print(f"Time = {pkt.data['IDX__SCI0TIME32'].derived_value}")  # Time in 20 ns intervals


                    print(f"Rice compression enabled = {bool(pkt.data['IDX__SCI0COMP'].raw_value)}")
                    compressed = bool(pkt.data['IDX__SCI0COMP'].raw_value)  # If we need to decompress the data


                    # self.header[evtnum][f"TimeIntervals"] = pkt.data['IDX__SCI0TIME32'].derived_value  # Store the number of 20 us intervals in the respective CDF "Time" variables
                    seconds_since_spacecraft_epoch = combine_coarse_fine_seconds(
                        pkt.data['SHCOARSE'].derived_value,
                        pkt.data['SHFINE'].derived_value,
                    )
                    print(
                        "Timestamp ="
                        f" {seconds_since_spacecraft_epoch} seconds since epoch (Midnight January 1st,"
                        f" {SPACECRAFT_EPOCH.year})"
                    )

                    utc_time = spacecraft_seconds_to_datetime(seconds_since_spacecraft_epoch)
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


        # Parse the waveforms according to the scitype present (high gain and low gain channels encode waveform data differently).
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
                        # decompressor = RiceGolombDecompressor(waveform)
                        
                        if(scitype[1] < 12):  # LS
                            nsamples = 8*(self.lspretrigblocks + self.lsposttrigblocks)
                            # copy = gpt_rice_Decode(waveform, True, nsamples)
                            # copy = rice_Decode(compressedFile, f"test.txt", False, nsamples)
                            # copy = decompressor.decompress(10)
                            # waveform = copy
                        else:  # HS
                            nsamples = 512*(self.hspretrigblocks + self.hsposttrigblocks) # - pkt.data['IDX__TXHDRSAMPDELAY']
                            # copy = gpt_rice_Decode(waveform, True, nsamples)
                            # copy = rice_Decode(compressedFile, f"test.txt", True, nsamples)
                            # copy = decompressor.decompress(12)
                            # waveform = copy

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
        conversion_factors = {
            "TOF H": 2.89e-4,
            "TOF M": 1.13e-2,
            "TOF L": 5.14e-4,
            "Ion Grid": 7.46e-4,
            "Target H": 1.63e-1,
            "Target L": 1.58e1,
        }
        unit_labels = {
            "TOF L": "pC/Δt",
            "TOF M": "pC/Δt",
            "TOF H": "pC/Δt",
            "Ion Grid": "pC",
            "Target L": "pC",
            "Target H": "pC",
        }

        def _build_stats_text(values, unit):
            return (
                f"Min = {np.min(values):.3g} [{unit}]\n"
                f"Avg = {np.mean(values):.3g} [{unit}]\n"
                f"Std = {np.std(values):.3g} [{unit}]\n"
                f"Max = {np.max(values):.3g} [{unit}]"
            )

        fig, ax = plt.subplots(nrows=6)  # Make this general
        fig.set_size_inches(18.5, 10.5)
        fig.subplots_adjust(top=0.9, bottom=0.08, left=0.08, right=0.78, hspace=0.4)
        for i, (k, v) in enumerate(packets.items()):  # k[0] = Event number, k[1] = channel name, v=waveform data

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

            
            print(f"Length of channel {k[1]} = {len(v)}")
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
                tof_l = np.asarray(packets[(k[0], "TOF L")], dtype=float) * conversion_factors["TOF L"]
                ax[0].plot(self.hstime, tof_l)
                ax[0].set_ylabel("TOF L [pC/Δt]", font="Times New Roman", fontsize=15, fontweight='bold')
                text = _build_stats_text(tof_l, unit_labels["TOF L"])
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
                
                tof_m = np.asarray(packets[(k[0], "TOF M")], dtype=float) * conversion_factors["TOF M"]
                ax[1].plot(self.hstime, tof_m)
                ax[1].set_ylabel("TOF M [pC/Δt]", font="Times New Roman", fontsize=15, fontweight='bold')
                text = _build_stats_text(tof_m, unit_labels["TOF M"])
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
                
                tof_h = np.asarray(packets[(k[0], "TOF H")], dtype=float) * conversion_factors["TOF H"]
                ax[2].plot(self.hstime, tof_h)
                ax[2].set_ylabel("TOF H [pC/Δt]", font="Times New Roman", fontsize=15, fontweight='bold')
                text = _build_stats_text(tof_h, unit_labels["TOF H"])
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

                ion_grid = np.asarray(packets[(k[0], "Ion Grid")], dtype=float) * conversion_factors["Ion Grid"]
                ax[3].plot(self.lstime, ion_grid)
                ax[3].set_ylabel("Ion Grid [pC]", font="Times New Roman", fontsize=15, fontweight='bold')
                text = _build_stats_text(ion_grid, unit_labels["Ion Grid"])
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
                
                    target_l = np.asarray(packets[(k[0], "Target L")], dtype=float) * conversion_factors["Target L"]
                    ax[4].plot(self.lstime, target_l)
                    ax[4].set_ylabel("Target LG [pC]", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = _build_stats_text(target_l, unit_labels["Target L"])
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
                    
                    target_h = np.asarray(packets[(k[0], "Target H")], dtype=float) * conversion_factors["Target H"]
                    ax[5].plot(self.lstime, target_h)
                    ax[5].set_ylabel("Target HG [pC]", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = _build_stats_text(target_h, unit_labels["Target H"])
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
                    target_h = np.asarray(packets[(k[0], "Target H")], dtype=float) * conversion_factors["Target H"]
                    ax[4].plot(self.lstime, target_h)
                    ax[4].set_ylabel("Target HG [pC]", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = _build_stats_text(target_h, unit_labels["Target H"])
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
                    
                    target_l = np.asarray(packets[(k[0], "Target L")], dtype=float) * conversion_factors["Target L"]
                    ax[5].plot(self.lstime, target_l)
                    ax[5].set_ylabel("Target LG [pC]", font="Times New Roman", fontsize=15, fontweight='bold')
                    text = _build_stats_text(target_l, unit_labels["Target L"])
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
                plt.show()
                plt.close()
                del fig, ax
                fig, ax = plt.subplots(nrows=6)  # Make this general
                fig.set_size_inches(17.5, 10.5)
                
    
    # ||
    # ||
    # || Write the waveform data 
    # || to an HDF5 file
    def write_to_hdf5(self, waveforms: dict, filename: str):
        output_dir = default_hdf5_dir(create=True)
        output_path = output_dir / Path(filename).name

        if output_path.exists():
            output_path.unlink()

        h = h5py.File(output_path, 'w')
        for k, v in waveforms.items():
            # print(np.array(v))
            # h.create_dataset(k, data=np.array(v, dtype=np.int8))
            # print(f"Time = {self.header}")
            if(f"/{k[0]}/Metadata/Epoch" not in h):
                h.create_dataset(f"/{k[0]}/Metadata/Epoch", data = self.header[(int(k[0]), 'Timestamp')])
            h.create_dataset(f"/{k[0]}/{k[1]}", data=np.array(v))
            if(k[1]=='TOF L'):
                h.create_dataset(f"/{k[0]}/Time (high sampling)", data=self.hstime)
            if(k[1]=='Ion Grid'):
                h.create_dataset(f"/{k[0]}/Time (low sampling)", data=self.lstime)
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
    return ints[:-1]

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
# def write_to_cdf(packets):
    
#     cdf_master = cdfread.CDF('imap_idex_l0-raw_0000000_v01.cdf')
#     if (cdf_master.file != None):
#     # Get the cdf's specification
#         info=cdf_master.cdf_info()
#         cdf_file=cdfwrite.CDF('./IDEX_SSIM.cdf',cdf_spec=info,delete=True)
#     # if (cdf_file.file == None):
#     #     cdf_master.close()
#     #     raise OSError('Problem writing file.... Stop')

#     # Get the global attributes
#     globalaAttrs=cdf_master.globalattsget(expand=True)
#     # Write the global attributes
#     cdf_file.write_globalattrs(globalaAttrs)
#     zvars=info['zVariables']
#     print('no of zvars=',len(zvars))
#     # Loop thru all the zVariables --> What are zvars vs rvars?
#     for x in range (0, len(zvars)):
#         # Get the variable's specification
#         varinfo=cdf_master.varinq(zvars[x])
#         print('Z =============>',x,': ', varinfo['Variable'])


# # Z =============> 0 :  Epoch
# # Z =============> 1 :  IDEX_Trigger
# # Z =============> 2 :  TOF_Low
# # Z =============> 3 :  TOF_Mid
# # Z =============> 4 :  TOF_High
# # Z =============> 5 :  Time_Low_SR
# # Z =============> 6 :  Time_High_SR
# # Z =============> 7 :  Target_Low
# # Z =============> 8 :  Target_High
# # Z =============> 9 :  Ion_Grid

#         if(varinfo['Variable']=="Epoch"):
#             vardata = None
#         if(varinfo['Variable']=="IDEX_Trigger"):
#             vardata = packets.header[(1,"Timestamp")]
#         if(varinfo['Variable']=="TOF_Low"):
#             print(len(np.array(packets.data[(1,"TOF L")])))
#             vardata = np.array(packets.data[(1,"TOF L")], float)
#         if(varinfo['Variable']=="TOF_Mid"):
#             vardata = np.array(packets.data[(1,"TOF M")])
#         if(varinfo['Variable']=="TOF_High"):
#             vardata = np.array(packets.data[(1,"TOF H")])
#         if(varinfo['Variable']=="Time_Low_SR"):
#             vardata = np.linspace(0, len(packets.data[(1,"Target L")]), len(len(packets.data[(1,"Target L")])))
#         if(varinfo['Variable']=="Time_High_SR"):
#             vardata = np.linspace(0, len(packets.data[(1,"TOF L")]), len(len(packets.data[(1,"Target L")])))
#         if(varinfo['Variable']=="Target_Low"):
#             vardata = np.array(packets.data[(1,"Target L")])
#         if(varinfo['Variable']=="Target_High"):
#             vardata = np.array(packets.data[(1,"Target H")])
#         if(varinfo['Variable']=="Ion_Grid"):
#             vardata = np.array(packets.data[(1,"Ion Grid")])
#         # Get the variable's attributes
#         varattrs=cdf_master.varattsget(zvars[x], expand=True)
#         if (varinfo['Sparse'].lower() == 'no_sparse'):
#             # A variable with no sparse records... get the variable data
#             # vardata= None
#             # Create the zVariable, write out the attributes and data
#             cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
#         else:
#             # A variable with sparse records...
#             # data is in this form [physical_record_numbers, data_values]
#             # physical_record_numbers (0-based) contains the real record
#             # numbers. For example, a variable has only 3 physical records
#             # at [0, 5, 10]:
#             varrecs=[0,5,10]

#             # vardata=None  # np.asarray([.,.,.,..])
#             # Create the zVariable, and optionally write out the attributes
#             # and data
#             cdf_file.write_var(varinfo, var_attrs=varattrs,
#                         var_data=[varrecs,vardata])
#     rvars=info['rVariables']
#     print('no of rvars=',len(rvars))
#     # Loop thru all the rVariables
#     for x in range (0, len(rvars)):
#         varinfo=cdf_master.varinq(rvars[x])
#         print('R =============>',x,': ', varinfo['Variable'])
#         varattrs=cdf_master.varattsget(rvars[x], expand=True)
#         if (varinfo['Sparse'].lower() == 'no_sparse'):
#             vardata=None
#             # Create the rVariable, write out the attributes and data
#             cdf_file.write_var(varinfo, var_attrs=varattrs, var_data=vardata)
#         else:
#             varrecs= None  # [.,.,.,..]
#             vardata= None  # np.asarray([.,.,.,..])
#             cdf_file.write_var(varinfo, var_attrs=varattrs,
#                         var_data=[varrecs,vardata])
#     cdf_master.close()
#     cdf_file.close()
    
# %%
# || Executable code: Connect to the OIS server, 
# || Inherited from Julie Barnum's SUDA code.
if __name__ == "__main__":
    try:
 
        data = []
        data_file = ''

        # SSIM1 was chosen ---> We will only have one SSIM for IDEX
        print('Connecting to the rack')
        HOST = ''  # Set to an empty string when on the rack
        # If not on the rack, set it to the rack's IP.
    

        PORT = 7514  # The port used by the server

        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((HOST, PORT))
        print(f'Connected to host {HOST} and port {PORT}')

        # signal.signal(signal.SIGINT, keyboard_interrupt_handler)


        try:
            # Going go read in data from the OIS socket server and append it to a data list
            while True:
                d = s.recv(1024)
                if d:
                   print(f'Appending data: \n{repr(d)}')
                   data.append(d)
                else:
                   print('No data received')
        except KeyboardInterrupt:
            print('\nFinished taking data, writing data to file now.')
            current_time = datetime.now().strftime("%m%d%Y_%H%M%S")  # Adding timestamp to filename
            # if not (os.path.exists(os.path.join(os.getcwd(), 'output_from_ois'))):  # If the "output_from_ois" folder doesn't exist,
            #     os.makedirs(os.path.join(os.getcwd(), 'output_from_ois'))  # Create it.
            data_file = os.path.join(os.getcwd(), 'output_from_ois', f'ois_output_{current_time}')  # Name the output file "ois_output_<timestamp>"
            with open(data_file, 'wb') as f:
                f.write(b''.join(data))
            print('Finished writing data to file.')

    except Exception as e:  # Read off server error. If connection is refused, make sure OASIS is running.
        error_st = e
        print(f'Reading from OASIS encountered a critical error: {error_st}')

    # Pass generated data file to the IDEXEvent struct
    packets = IDEXEvent(data_file)
    # print(packets.data.keys())
    packets.plot_all_data(packets.data, data_file)
    packets.write_to_hdf5(packets.data, data_file+'.h5')
    # write_to_cdf(packets)

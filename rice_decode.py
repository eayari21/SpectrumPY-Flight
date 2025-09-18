#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A Python command line tool to decompress IDEX packets.
__author__      = Ethan Ayari & Gavin Medley, 
Institute for Modeling Plasmas, Atmospheres and Cosmic Dust

Works with Python 3.8.10
"""
import argparse
import bitstring
import numpy as np
import matplotlib.pyplot as plt

# || Translation of Corrine's code
def idex_rice_Decode(input, nBit10, sampleCount):
    # 10 bit data has 512 samples per frame
    # 12 bit data has 1024 samples per frame
    # SF_SIZE is the compression block size
    # FRAME_SIZE is the expected amount of data

    # Constants:
    k_bits = 4
    if nBit10:
        n_bits = 10
    else:
        n_bits = 12
    SF_SIZE = 64
    FRAME_SIZE = sampleCount
    SF_PER_FRAME = FRAME_SIZE / SF_SIZE

    # Read ASCII Hex data from file
    indata = []

    for line in input:
        words = line.split()
        for s in words:
            indata.append(int(s, 16))

    # Convert data to an array of bytes
    byte_data = bytearray()
    for n in indata:
        for byte_index in range(3, -1, -1):
            byte_data.append((n >> (byte_index * 8)) & 0xFF)

    print('Read', len(byte_data) // 4, 'words from input file')

    # Create a ConstBitStream from the byte data
    bits = bitstring.ConstBitStream(bytes=byte_data)

    outdata = []
    sample_count = 0
    total_count = 0
    sf_count = 0
    iternum = 0
    while iternum < len(bits) and sf_count < SF_PER_FRAME:
        iternum += 1

        # Next two bits are the predictor select bits
        psel = bits.read(2)

        if psel.int > 1:
            k = bits.read(k_bits)

        sample_count = 0
        sf_data = []
        # Ensure byte alignment before accessing bytepos
        bits.bytealign()
        start_byte = bits.bytepos
        compress_size = 0
        while sample_count < SF_SIZE and iternum < len(bits):
            if sample_count == 0:
                # Read warmup sample
                d1 = bits.read(n_bits)
                sf_data.append(d1)
                sample_count += 1

                # Constant value frame contains 512 copies of the same data value
                if psel == 0:
                    sf_data.extend([d1] * (SF_SIZE - 1))
                    sample_count = SF_SIZE

            elif psel == 1 or (sample_count == 1 and psel == 3):
                d1 = bits.read(n_bits)
                sf_data.append(d1)
                sample_count += 1
            else:
                q = 0
                while True:
                    bit = bits.read(1)
                    if bit:
                        break
                    q += 1

                if q & 0x1:
                    q = int(-((q + 1) / 2))
                else:
                    q = int(q / 2)

                r = bits.read(k + 1)
                d1 = (q << (k + 1)) + r

                if psel == 2:
                    d1 = d1 + sf_data[sample_count - 1]
                elif sample_count > 1 and psel == 3:
                    d1 = d1 + 2 * sf_data[sample_count - 1] - sf_data[sample_count - 2]

                if int(d1.uint) > 2 ** n_bits or int(d1.uint) < -2 ** n_bits:
                    print('Overflow Error at byte offset', bits.bytepos)
                    print('Current byte =', hex(byte_data[bits.bytepos]))
                    print('Current Sample Count =', total_count + sample_count)
                    print('k =', k, 'q =', q, 'r =', r, 'd1 =', d1)
                    print('DataOut =', sf_data)
                    return

                sf_data.append(d1)
                sample_count += 1

        # Align your bits
        length = len(bits)
        alignment_bits = (-length) % 8  # Calculate the number of additional bits needed for alignment
        padding = bitstring.BitStream(bin='0' * alignment_bits)  # Create a new BitStream with padding bits
        aligned_bits = bits + padding  # Concatenate the existing bits with the padding bits
        byte_position = aligned_bits.bytepos

        stop_byte = byte_position % 8
        compressed_size = stop_byte - start_byte
        compression_rate = (compressed_size * 8) / (SF_SIZE * n_bits)
        outdata.extend(sf_data)
        total_count += sample_count
        sf_count += 1

    print('Outdata length =', len(outdata))
    

    print('Decoded', len(outdata), 'samples from compressed data')

    if len(bits) - iternum < FRAME_SIZE - total_count:
        print('**** ERROR: End of file reached before', FRAME_SIZE, 'samples decoded ****')

    return(outdata)


# || Translation of Corrine's code
def rice_Decode(sourceFile, destinationFile, nBit10, sampleCount):
    # 10 bit data has 512 samples per frame
    # 12 bit data has 1024 samples per frame
    # SF_SIZE is the compression block size
    # FRAME_SIZE is the expected amount of data

    # Constants:
    k_bits = 4
    if nBit10:
        n_bits = 10
    else:
        n_bits = 12
    SF_SIZE = 64
    FRAME_SIZE = sampleCount
    SF_PER_FRAME = FRAME_SIZE / SF_SIZE

    # Read ASCII Hex data from file
    indata = []
    with open(sourceFile, 'rb') as infile:
        for line in infile:
            words = line.split()
            for s in words:
                indata.append(int(s, 16))

    # Convert data to an array of bytes
    byte_data = bytearray()
    for n in indata:
        for byte_index in range(3, -1, -1):
            byte_data.append((n >> (byte_index * 8)) & 0xFF)

    print('Read', len(byte_data) // 4, 'words from input file')

    # Create a ConstBitStream from the byte data
    bits = bitstring.ConstBitStream(bytes=byte_data)

    outdata = []
    sample_count = 0
    total_count = 0
    sf_count = 0
    iternum = 0
    while iternum < len(bits) and sf_count < SF_PER_FRAME:
        iternum += 1

        # Next two bits are the predictor select bits
        psel = bits.read(2)

        if psel.int > 1:
            k = bits.read(k_bits)

        sample_count = 0
        sf_data = []
        # Ensure byte alignment before accessing bytepos
        bits.bytealign()
        start_byte = bits.bytepos
        compress_size = 0
        while sample_count < SF_SIZE and iternum < len(bits):
            if sample_count == 0:
                # Read warmup sample
                d1 = bits.read(n_bits)
                sf_data.append(d1)
                sample_count += 1

                # Constant value frame contains 512 copies of the same data value
                if psel == 0:
                    sf_data.extend([d1] * (SF_SIZE - 1))
                    sample_count = SF_SIZE

            elif psel == 1 or (sample_count == 1 and psel == 3):
                d1 = bits.read(n_bits)
                sf_data.append(d1)
                sample_count += 1
            else:
                q = 0
                while True:
                    bit = bits.read(1)
                    if bit:
                        break
                    q += 1

                if q & 0x1:
                    q = int(-((q + 1) / 2))
                else:
                    q = int(q / 2)

                r = bits.read(k + 1)
                d1 = (q << (k + 1)) + r

                if psel == 2:
                    d1 = d1 + sf_data[sample_count - 1]
                elif sample_count > 1 and psel == 3:
                    d1 = d1 + 2 * sf_data[sample_count - 1] - sf_data[sample_count - 2]

                if int(d1.uint) > 2 ** n_bits or int(d1.uint) < -2 ** n_bits:
                    print('Overflow Error at byte offset', bits.bytepos)
                    print('Current byte =', hex(byte_data[bits.bytepos]))
                    print('Current Sample Count =', total_count + sample_count)
                    print('k =', k, 'q =', q, 'r =', r, 'd1 =', d1)
                    print('DataOut =', sf_data)
                    return

                sf_data.append(d1)
                sample_count += 1

        # Align your bits
        length = len(bits)
        alignment_bits = (-length) % 8  # Calculate the number of additional bits needed for alignment
        padding = bitstring.BitStream(bin='0' * alignment_bits)  # Create a new BitStream with padding bits
        aligned_bits = bits + padding  # Concatenate the existing bits with the padding bits
        byte_position = aligned_bits.bytepos

        stop_byte = byte_position % 8
        compressed_size = stop_byte - start_byte
        compression_rate = (compressed_size * 8) / (SF_SIZE * n_bits)
        outdata.extend(sf_data)
        total_count += sample_count
        sf_count += 1

    print('Outdata length =', len(outdata))
    with open(destinationFile, 'ab') as fpHex:
        for data in outdata:
            fpHex.write(data.to_bytes(2, byteorder='big'))

    print('Decoded', len(outdata), 'samples from compressed data')

    if len(bits) - iternum < FRAME_SIZE - total_count:
        print('**** ERROR: End of file reached before', FRAME_SIZE, 'samples decoded ****')


if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument("--sourcefile", "-s", type=str, required=True)
    aparser.add_argument("--targetfile", "-t", type=str, required=True)
    args = aparser.parse_args()

    compressedFile = args.sourcefile
    destinationFile = args.targetfile
    sampleCount = 500
    rice_Decode(compressedFile, destinationFile, True, sampleCount)
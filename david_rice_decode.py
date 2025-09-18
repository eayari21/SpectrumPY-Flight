import numpy as np
import matplotlib.pyplot as plt


class bitstream():
    def __init__(self, b):
        # Stream to read from or write to. It probably only makes sense to
        # either read or write, but never both, on the same stream
        self.b = b
        self.buf = 0
        self.size = len(b)
        self.byteptr = 0
        self.bitptr = 7
        self.neverread = True
        self.neverwritten = True
        self.eos_flag = False

    def read(self, nbits):
        # Do it the stupid way first, shift one a bit at a time
        result = 0
        if self.neverread:
            self.neverread = False
            self.buf = self.b[self.byteptr]
            self.bitptr = 7

        for _ in range(nbits):
            if self.bitptr < 0:
                self.bitptr = 7
                if self.byteptr < self.size - 1:
                    self.byteptr += 1
                    self.buf = self.b[self.byteptr]
                else:
                    self.eos_flag = True
            bufbit = (self.buf >> self.bitptr) & 1
            result = (result << 1) | bufbit
            self.bitptr = self.bitptr - 1
        return result


    def read_byte(self):
        result = 0

        self.neverread = False
        self.bitptr = 7
        result = self.b[self.byteptr]
        if self.byteptr < self.size - 1:
            self.byteptr += 1
            self.buf = self.b[self.byteptr]
        return result


    def read_signed(self, nbits):
        # Do it the stupid way first, shift on a bit at a time
        result = 0
        if self.neverread:
            self.neverread = False
            self.buf = self.b[self.byteptr]
            self.bitptr = 7
        for _ in range(nbits):
            if self.bitptr < 0:
                self.bitptr = 7
                self.byteptr += 1
                self.buf = self.b[self.byteptr]
            bufbit = (self.buf >> self.bitptr) & 1
            result = (result << 1) | bufbit
            self.bitptr = self.bitptr - 1
            if (result > 2**(nbits-1)):
                result = result - 2**nbits
        return result

    def eos(self):
        if self.eos_flag:
            return 1
        else:
            return 0

    def get_byte_ptr(self):
        return self.byteptr

    def get_bit_ptr(self):
        return self.bitptr



    # def write(self, val, nbits, meaning=""):
    #     # Do it the stupid way first, shift on a bit at a time
    #     if verbose and nbits>0: print("Appending value %d, %d bits"%(val,nbits))
    #     self.neverwritten = False
    #     for i in range(nbits):
    #         valbit = (val >> (nbits - 1 - i)) & 1
    #         #if verbose: print(valbit, end='', flush=True)
    #         self.buf |= valbit << self.ptr
    #         self.ptr = self.ptr - 1
    #         if self.ptr < 0:
    #             trace_write(self.f,bytes([self.buf]),meaning="Writing bit data")
    #             if self.buf != 0xff:
    #                 self.ptr = 7
    #             else:
    #                 # Do a zero-byte stuffing thing so as to follow the JPEG
    #                 # specification
    #                 self.ptr = 6
    #             self.buf = 0
#        if meaning is None:
#            if verbose: print(" ", end='', flush=True)
#        else:
#            if verbose: print(" "+meaning)

    def flush(self):
        if self.bitptr != 0 and not self.neverwritten:
            trace_write(self.f,bytes([self.buf]),"Flushing")
        self.buf = 0
        self.bitptr = 7
        self.byteptr = 7
        self.neverread = True
        self.neverwritten = True

    def close(self):
        self.flush()
        return self.f.close()



def rice_Decode():

    #10 bit data has 512 samples per frame
    #12 bit data has 1024 samples per frame

    #Constants:
    k_bits = 4
    n_bits = 10
    SF_SIZE = 64
    FRAME_SIZE = 1024
    SF_PER_FRAME = FRAME_SIZE / SF_SIZE


    #Read ASCII Hex data from file

    #filename="C:/LASP_SVN/projects/SUDA/trunk/ProcBoard/FPGA/SUDAProc/sim/TestData/Compress Test/output_07-28-2020-13-49-19_sciData_compressedDT5/event_0_000001e0478d/HS-ADC-0-Channel-I-High-Gain/HS-ADC-0-Channel-I-High-Gain.compressed"

    filename="C:/LASP_SVN/projects/SUDA/trunk/ProcBoard/FPGA/SUDAProc/sim/TestData/output_08-06-2020-14-31-23_sciData_2020_219_10_54_55_compressed/LS-ADC-2-Channel-iCSA.compressed"
    filename = "C:/LASP_SVN/projects/SUDA/trunk/ProcBoard/FPGA/SUDAProc/sim/TestData/output_08-06-2020-14-31-23_sciData_2020_219_10_54_55_compressed/LS-ADC-2-Channel-iCSA.compressed"

    filename = "C:/LASP_SVN/projects/SUDA/trunk/ProcBoard/FPGA/SUDAProc/sim/results/rice_compression_output_10bit_hex.txt"



    #filename2 = "C:/LASP_SVN/projects/SUDA/trunk/ProcBoard/FPGA/SUDAProc/sim/TestData/Compress Test/output_07-28-2020-13-49-19_sciData_compressedDT5/event_0_000001e0478d/HS-ADC-0-Channel-I-High-Gain/HS-ADC-0-Channel-I-High-Gain.compressed.decode.txt"
    filename2="C:/LASP_SVN/projects/SUDA/trunk/ProcBoard/FPGA/SUDAProc/sim/TestData/output_08-06-2020-14-31-23_sciData_2020_219_10_54_55_compressed/event_0_000001e0478d/HS-ADC-0-Channel-I-High-Gain/HS-ADC-0-Channel-I-High-Gain.compressed_decode.txt"

    filename2 = "C:/LASP_SVN/projects/SUDA/trunk/ProcBoard/FPGA/SUDAProc/sim/TestData/output_08-06-2020-14-31-23_sciData_2020_219_10_54_55_compressed/event_0_000001e0478d/LS-ADC-0-Channel-vCSA/LS-ADC-0-Channel-vCSA.compressed_decode.txt"


    LOGFILE = open(filename2, 'w')

    #Read in data from VHDL testbench as 32 bit ASCII hex values
    indata = []
    with open(filename, 'r') as infile:
      for line in infile:
          words = line.split()
          for s in words:
              try:
                  indata.append(int(s, 16))
              except ValueError:
                  print(s)
              print(s)


    #convert data to an array of bytes
    byte_data = []
    byte_index = 3
    for n in indata:
        for byte_index in range(4):
            byte_data.append((n>>((3 - byte_index)*8)) & 0xFF)

    print('Read ', len(byte_data)/4, ' words from input file')
    bits = bitstream(byte_data)


    # print(byte_data)
    # for i in range(20):
    #     print(bits.read(7))





    outdata = []

    sample_count = 0





    total_count = 0
    sf_count = 0
    #while(sample_count < SF_SIZE) and (not bits.eos()):
    while not bits.eos() and (sf_count < SF_PER_FRAME):

        print(sf_count)

        #Next two bits are the predictor select bits
        psel = bits.read(2)
        print('Predictor Selector = ', psel, ' byte offset ', bits.get_byte_ptr() )
        #if psel == 0:
         #   print(' Unexpected')

        if psel > 1:
            k = bits.read(k_bits)
            print('k = ', k)


        sample_count = 0
        sf_data = []
        start_byte =  bits.get_byte_ptr();
        compress_size = 0
        while (sample_count < SF_SIZE) and (not bits.eos()):
            #print('sample count = ', sample_count)
            #print(len(outdata))
            #print(outdata)
            if sample_count == 0:
                #Read warmup sample
                d1 = bits.read(n_bits)
                sf_data.append(d1)
                sample_count += 1

                #Constant value frame contains 512 copies of the same data value
                if psel==0:
                    print(' constant frame of value ', hex(d1))
                    for i in range(SF_SIZE-1):
                        sf_data.append(d1)

                    sample_count  = SF_SIZE
                else:
                    pass
                    #print('Warmup sample #1 ', hex(d1))


            elif (psel == 1) or ((sample_count == 1) and (psel == 3)):
                d1 = bits.read(n_bits)
                sf_data.append(d1)
                #print('Warmup sample #2 ', hex(d1))
                sample_count += 1
            else:

                q = 0
                while(bits.read(1) == 0):
                    q += 1
                #print('q = ', q)

                if q == 47:
                    d1=bits.read_signed(n_bits+2)
                    #d1 = bits.read_signed(14)
                else:
                    if ((q & 0x1)):
                        q = int(-((q+1)/2))
                    else:
                        q = int(q/2)

                    r = bits.read(k+1)
                    #print('r = ', r)
                    d1 = (q << (k+1))+r

                if psel == 2:
                    d1 = d1 + sf_data[sample_count-1]
                elif (sample_count > 1) and (psel == 3):
                    #print('sf_data[sample_count-2]=',hex(sf_data[(sample_count-2)]), '  sf_data[sample_count-1]',hex(sf_data[(sample_count-1)]))
                    #print(sample_count,' q = ', q, '  r = ', r, ' d1 = ', d1)
                    d1 = d1 + 2*sf_data[(sample_count-1)] - sf_data[(sample_count-2)]
                    #print('                        d1 = ', hex(d1))
                if (d1 > 2**n_bits) or (d1 < -2**n_bits):
                    print('Overflow Error at byte offset ',bits.get_byte_ptr() )
                    print('Current byte = ', hex(byte_data[bits.get_byte_ptr()]) )
                    print('Current Sample Count =  ', total_count + sample_count)
                    print('k = ', k, ' q= ', q, ' r = ', r, ' d1 = ', d1)
                    print('DataOut = ', sf_data)
                    plt.plot(outdata)
                    plt.show()
                    exit()
                sf_data.append(d1)
                #print(sf_data)
                sample_count += 1

        stop_byte =  bits.get_byte_ptr();
        compressed_size = stop_byte - start_byte
        compression_rate = (compressed_size*8)/(SF_SIZE*n_bits)
        print('Subframe compression rate = ', compression_rate)
        outdata.extend(sf_data)
        total_count += sample_count

        # if (sf_count < SF_PER_FRAME-1):
        #     sf_count += 1
        # else:
        #     sf_count = 0
        sf_count += 1

    # Write data samples to a text file as base16 integers
    for j in range(len(outdata)):
        LOGFILE.write('%03X\n' % (outdata[j]))


    print('Decoded ', len(outdata), ' samples from compressed data')

    if bits.eos() and (len(outdata) < FRAME_SIZE):
        print('**** ERROR:  End of file reached before' , FRAME_SIZE, 'samples decoded ****')

    print(outdata[-16:])
    plt.plot(outdata)
    plt.show()

    LOGFILE.close()



if __name__ == "__main__":
    rice_Decode()
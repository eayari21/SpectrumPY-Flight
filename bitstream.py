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
            if (result > 2 ** (nbits - 1)):
                result = result - 2 ** nbits
        return result

    def eos(self):
        if self.eos_flag:
            return 1
        else:
            return 0
        # if self.byteptr < self.size - 2:
        #     return 0
        # else:
        #     return 1

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
            trace_write(self.f, bytes([self.buf]), "Flushing")
        self.buf = 0
        self.bitptr = 7
        self.byteptr = 7
        self.neverread = True
        self.neverwritten = True

    def close(self):
        self.flush()
        return self.f.close()


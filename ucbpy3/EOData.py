
#########
# Imports
import os
import struct
import numpy as np
from sys import stderr

# EOD Header structure size (bytes) and format
HEADER_STRUCT_SIZE = 48
HEADER_STRUCT_FORMAT = 'i15cxffiiff4x'

# Format prefix for endian-ness
ENDIAN_PREFIX = { 'l': '<', 'b': '>' }

# Padding size at the beginning of the eoh file
EOH_PAD_SIZE = 16

#
# ---
#

class EOHeader:
    """
    Handler class for .eoh waveform data header files - currently read-only
    """
    def __init__(self, data_path, cmt_name, endian = 'l'):
        self.data_path = data_path
        self.cmt_name = cmt_name
        self.endian = endian
    def open_for_reading(self):
        file_name = self.data_path + '/' + self.cmt_name + '.eoh'
        if not os.path.exists(file_name):
            stderr.write('EOHeader.open_for_reading(): %s does not exist!\n' % (file_name))
            return False
        else:
            self.f_eoh = open(file_name, 'rb')
            self.f_eoh.seek(EOH_PAD_SIZE, os.SEEK_SET)
            return True
    def close(self):
        self.f_eoh.close()
    def __load_next_header__(self):
        raw_struct_str = self.f_eoh.read(HEADER_STRUCT_SIZE)
        if raw_struct_str:
            header_vars = struct.unpack(ENDIAN_PREFIX[self.endian] + HEADER_STRUCT_FORMAT, raw_struct_str)
            return { 'rcno': header_vars[0],
                     'netwk': ''.join(header_vars[1:5]),
                     'chnnl': ''.join(header_vars[5:9]),
                     'stn': ''.join(header_vars[9:13]),
                     'compnt': header_vars[13],
                     'dttype': header_vars[14],
                     'reftime': header_vars[15],
                     'starttm': header_vars[16],
                     'smplintv': header_vars[17],
                     'ndata': header_vars[18],
                     'locatn': header_vars[19],
                     'slat': header_vars[20],
                     'slon': header_vars[21] }
        else:
            return None
    def __iter__(self):
        return self
    def next(self):
        header = self.__load_next_header__()
        if header:
            return header
        else:
            self.reset()
            raise StopIteration
    def reset(self):
        self.f_eoh.seek(EOH_PAD_SIZE, os.SEEK_SET)

class EOData:
    """
    Handler class for .eod waveoform data trace files - currently read-only
    """
    def __init__(self, data_path, cmt_name, endian = 'l'):
        self.data_path = data_path
        self.cmt_name = cmt_name
        self.endian = endian
    def open_for_reading(self):
        file_name = self.data_path + '/' + self.cmt_name + '.eod'
        if not os.path.exists(file_name):
            stderr.write('EOData.open_for_reading(): %s does not exist!\n' % (file_name))
            return False
        else:
            self.f_eod = open(file_name, 'rb')
            return True
    def close(self):
        self.f_eod.close()
    def fetch_data(self, header):
        self.f_eod.seek(header['locatn'], os.SEEK_SET)
        data = np.fromfile(self.f_eod, dtype = ENDIAN_PREFIX[self.endian] + 'f4', count = header['ndata'], sep = '')
        return data


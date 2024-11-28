#########
# Imports
import os
import struct
import numpy as np
from sys import exit, stderr

# Format prefix for endian-ness
ENDIAN_PREFIX = { 'l': '<', 'b': '>' }

####
# First, we define a series of handlers for converting struct.unpack()
# output tuples to dicts...

#
# handler for unpacking title_st
def handler_title_st(struct_tuple):
    moment = { 'rr': struct_tuple[11],
               'tt': struct_tuple[12],
               'pp': struct_tuple[13],
               'rt': struct_tuple[14],
               'rp': struct_tuple[15],
               'tp': struct_tuple[16] }
    return { 'event': ''.join(struct_tuple[:8]),
             'theta': struct_tuple[8],
             'phi': struct_tuple[9],
             'depth': struct_tuple[10],
             'moment': moment,
             'dt': struct_tuple[17] }

#
# handler for unpacking tracehdr_st
def handler_tracehdr_st(struct_tuple):
    return { 'stn': ''.join(struct_tuple[:4]),
             'locatn': struct_tuple[4],
             'theta': struct_tuple[5],
             'phi': struct_tuple[6],
             'delta': struct_tuple[7],
             'az': struct_tuple[8],
             'dip': struct_tuple[9],
             'smplintv': struct_tuple[10],
             'w1': struct_tuple[11],
             'w2': struct_tuple[12],
             'w3': struct_tuple[13],
             'w4': struct_tuple[14],
             'id': struct_tuple[15],
             'reftime': struct_tuple[16],
             'comp': struct_tuple[17],
             'netwk': ''.join(struct_tuple[18:22]),
             'chnnl': ''.join(struct_tuple[22:26]),
             'extr': ''.join(struct_tuple[26:30]) }
#
# handler for unpacking tracehdrH_st
def handler_tracehdrH_st(struct_tuple):
    return { 'stn': ''.join(struct_tuple[:4]),
             'locatn': struct_tuple[4],
             'theta': struct_tuple[5],
             'phi': struct_tuple[6],
             'delta': struct_tuple[7],
             'az': struct_tuple[8],
             'dip': struct_tuple[9],
             'smplintv': struct_tuple[10],
             'w1': struct_tuple[11],
             'w2': struct_tuple[12],
             'w3': struct_tuple[13],
             'w4': struct_tuple[14],
             'id': struct_tuple[15],
             'reftime': struct_tuple[16],
             'comp': struct_tuple[17],
             'netwk': ''.join(struct_tuple[18:22]),
             'chnnl': ''.join(struct_tuple[22:26]),
             'extr': ''.join(struct_tuple[26:30]),
             'locatnA': struct_tuple[30] }


def handler_packhdrbr_st(struct_tuple):
    return { 'phase': struct_tuple[0],
             'ndata': struct_tuple[1],
             'id': struct_tuple[2],
             't0': struct_tuple[3],
             'gv1': struct_tuple[4],
             'gv2': struct_tuple[5],
             'pv1': struct_tuple[6],
             'pv2': struct_tuple[7],
             'rmsd': struct_tuple[8],
             'rmsr': struct_tuple[9],
             'rmss': struct_tuple[10],
             'weight': struct_tuple[11] }

def handler_packhdrmp_st(struct_tuple):
    print 'TBD'
    exit(1)

def handler_packhdrcr_st(struct_tuple):
    print 'TBD'
    exit(1)


####
# Now, we assemble the structure sizes, formats, and handlers into dicts

#
# wph/H title
TITLE_STRUCT_FORMAT  = { 'title_st': '8c3f6ff' }
TITLE_STRUCT_HANDLER = { 'title_st': handler_title_st }

#
# wph/H trace header
TRACEHDR_STRUCT_FORMAT  = { 'tracehdr_st': '4ci10fh10c4c',
                            'tracehdrH_st': '4ci10fh10c4ci' }
TRACEHDR_STRUCT_HANDLER = { 'tracehdr_st': handler_tracehdr_st,
                            'tracehdrH_st': handler_tracehdrH_st }

#
# wpd packet header
PACKHDR_STRUCT_FORMAT  = { 'packhdrbr_st': '3h2x9f',
                           'packhdrmp_st': '3hcx10f',
                           'packhdrcr_st': '3h2x12x' }
PACKHDR_STRUCT_HANDLER = { 'packhdrbr_st': handler_packhdrbr_st,
                           'packhdrmp_st': handler_packhdrmp_st,
                           'packhdrcr_st': handler_packhdrcr_st }

####
# Now we define classes to wrap wph and wpd file instances

class WPHeader:
    def __init__(self, data_path, cmt_name, extension = 'wph', endian = 'l'):
        #
        self.data_path = data_path
        self.cmt_name = cmt_name
        self.extension = extension
        self.endian = endian
        #
        title_type = 'title_st'
        tracehdr_type = 'tracehdr_st'
        #
        self.__init_format__(title_type, tracehdr_type)
    def __init_format__(self, title_type, tracehdr_type):
        self.title_format = ENDIAN_PREFIX[self.endian] + TITLE_STRUCT_FORMAT[title_type]
        self.title_size = struct.calcsize(self.title_format)
        self.title_handler = TITLE_STRUCT_HANDLER[title_type]
        #
        if self.extension[-1] == 'H':
            tracehdr_type = 'tracehdrH_st'
        self.tracehdr_format = ENDIAN_PREFIX[self.endian] + TRACEHDR_STRUCT_FORMAT[tracehdr_type]
        self.tracehdr_size = struct.calcsize(self.tracehdr_format)
        self.tracehdr_handler = TRACEHDR_STRUCT_HANDLER[tracehdr_type]
    def open_for_reading(self):
        file_name = self.data_path + '/' + self.cmt_name + '.' + self.extension
        if not os.path.exists(file_name):
            stderr.write('WPHeader.open_for_reading(): %s does not exist!\n' %
                             (file_name))
            return False
        else:
            self.f_wph = open(file_name, 'rb')
            self.__load_title__()
            return True
    def close(self):
        self.f_wph.close()
    def get_title(self):
        return self.title
    def __load_title__(self):
        raw_struct_str = self.f_wph.read(self.title_size)
        self.title = self.title_handler(struct.unpack(self.title_format, raw_struct_str))
    def __load_next_trace_header__(self):
        raw_struct_str = self.f_wph.read(self.tracehdr_size)
        if len(raw_struct_str) == self.tracehdr_size:
            tracehdr = self.tracehdr_handler(struct.unpack(self.tracehdr_format, raw_struct_str))
        else:
            tracehdr = None
        return tracehdr
    def __iter__(self):
        return self
    def next(self):
        tracehdr = self.__load_next_trace_header__()
        if tracehdr:
            return tracehdr
        else:
            self.reset()
            raise StopIteration
    def reset(self):
        self.f_wph.seek(self.title_size, os.SEEK_SET)
        
class WPData:
    def __init__(self, data_path, cmt_name, extension = 'wpd', endian = 'l'):
        self.data_path = data_path
        self.cmt_name = cmt_name
        self.extension = extension
        self.endian = endian
        #
        packhdr_type = 'packhdrbr_st'
        #
        self.__init_format__(packhdr_type)
    def __init_format__(self, packhdr_type):
        self.packhdr_format = ENDIAN_PREFIX[self.endian] + PACKHDR_STRUCT_FORMAT[packhdr_type]
        self.packhdr_size = struct.calcsize(self.packhdr_format)
        self.packhdr_handler = PACKHDR_STRUCT_HANDLER[packhdr_type]
    def open_for_reading(self):
        file_name = self.data_path + '/' + self.cmt_name + '.' + self.extension
        if not os.path.exists(file_name):
            stderr.write('WPData.open_for_reading(): %s does not exist!\n' %
                             (file_name))
            return False
        else:
            self.f_wpd = open(file_name, 'rb')
            return True
    def close(self):
        self.f_wpd.close()
    def load_wavepackets(self, tracehdr):
        self.f_wpd.seek(tracehdr['locatn'], os.SEEK_SET)
        raw_struct_str = self.f_wpd.read(self.packhdr_size)
        wavepackets = []
        while raw_struct_str:
            packhdr = self.packhdr_handler(struct.unpack(self.packhdr_format, raw_struct_str))
            if packhdr['id'] == tracehdr['id']:
                wavepackets.append((packhdr,
                                    np.fromfile(self.f_wpd,
                                                dtype = ENDIAN_PREFIX[self.endian] + 'f4',
                                                count = packhdr['ndata'], sep = '')))
                raw_struct_str = self.f_wpd.read(self.packhdr_size)
            else:
                break
        return wavepackets

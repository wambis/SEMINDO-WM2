#!/usr/bin/env python

####
# Imports
import struct
import numpy as np

# Format prefix for endian-ness
_ENDIAN_PREFIX = { 'l': '<', 'b': '>' }

# Padding size at the beginning of the eoh file
_EOH_PAD_SIZE = 16

# eod trace header format
_HEADER_STRUCT_FORMAT = 'i16cffiiff4c'

encoding = 'ISO-8859-1'#utf-8
####
# handlers for packing / unpacking eoh header structs
def _handler_header_struct_r(header_vars):
    return { 'rcno':     header_vars[0]
           , 'netwk':    ''.join([elem.decode(encoding) for elem in header_vars[1:5]])
           , 'chnnl':    ''.join([elem.decode(encoding) for elem in header_vars[5:9]])
           , 'stn':      ''.join([elem.decode(encoding) for elem in header_vars[9:13]])
           , 'compnt':   header_vars[13].decode(encoding)
           , 'dttype':   header_vars[14].decode(encoding)
           , 'reftime':  header_vars[15].decode(encoding)
           , 'extr0':    header_vars[16].decode('latin')
#           , 'extr0':    header_vars[16].decode(encoding)
           , 'starttm':  header_vars[17]
           , 'smplintv': header_vars[18]
           , 'ndata':    header_vars[19]
           , 'locatn':   header_vars[20]
           , 'slat':     header_vars[21]
           , 'slon':     header_vars[22]
           , 'extr1':    "".join([elem.decode('latin') for elem in header_vars[23:]])}
#           , 'extr1':    "".join([elem.decode(encoding) for elem in header_vars[23:]])}

def _handler_header_struct_w(struct_dict):
    stn   = [c for c in struct_dict['stn']]
    netwk = [c for c in struct_dict['netwk']]
    chnnl = [c for c in struct_dict['chnnl']]
    extr1 = [c for c in struct_dict['extr1']]
    order = ['compnt', 'dttype', 'reftime', 'extr0', 'starttm', 'smplintv', 'ndata', 'locatn', 'slat', 'slon']
    return tuple([struct_dict['rcno']] + netwk + chnnl + stn + [struct_dict[k] for k in order] + extr1)
_HEADER_STRUCT_HANDLER = {'r': _handler_header_struct_r, 'w': _handler_header_struct_w}


####
# Finally, the actual eod trace data class implementation

class EOData:

    def __init__(self, headers_file, data_file, endian_r = 'l', endian_w = 'l'):
        # named 
        self._headers_file = headers_file
        self._data_file = data_file
        # keyword args
        self._endian_r = endian_r
        self._endian_w = endian_w
        # init
        self._open = False

    def open(self, mode='r'):
        # check / store the mode
        assert mode in ['r','w']
        self._mode = mode
        # attempt to open both the header and data files
        try:
            self._f_headers = open(self._headers_file, self._mode + 'b')
        except IOError as e:
            print( "Note: [I/O error] cannot open %s - %s" % (self._headers_file, e.strerror))
            return False
        try:
            self._f_data = open(self._data_file, self._mode + 'b')
        except IOError as e:
            print( "Note: [I/O error] cannot open %s - %s" % (self._data_file, e.strerror))
            return False
        # success
        self._open = True
        return True

    def close(self):
        assert self._open == True
        self._f_headers.close()
        self._f_data.close()
        self._open = False

    def read_headers(self):
        """
        Returns all eod trace headers from the target eoh file
        """
        # check
        assert self._open
        assert self._mode == 'r'
        # seek to start of headers
        self._f_headers.seek(_EOH_PAD_SIZE)
        # read header structs
        fmt = _ENDIAN_PREFIX[self._endian_r] + _HEADER_STRUCT_FORMAT
        struct_size = struct.calcsize(fmt)
        raw = self._f_headers.read(struct_size)
        headers = []
        while raw:
            headers.append(_HEADER_STRUCT_HANDLER['r'](struct.unpack(fmt, raw)))
            raw = self._f_headers.read(struct_size)
        return headers

    def read_data(self, header):
        """
        read_data(self, header)

        Returns eod trace associated with the supplied trace header
        """
        # check
        assert self._open
        assert self._mode == 'r'
        # jump to data start
        self._f_data.seek(header['locatn'])
        # read and return trace
        return np.fromfile(self._f_data, dtype=_ENDIAN_PREFIX[self._endian_r] + 'f4', count=header['ndata'])

    def write_headers_and_data(self, traces, verbose=True):
        """
        write_headers_and_data(self, traces, verbose=True)

        Builds a eoh / eod file pair from the supplied list of traces.
        
          Trace list is of the form [(tracehdr1,data1),...]

          Consistent header 'locatn' parameters are ensured
        """
        # check
        assert self._open
        assert self._mode == 'w'
        # will fill in locatn entries as needed ...
        # make sure we are rewound
        self._f_headers.seek(_EOH_PAD_SIZE)
        self._f_data.seek(0)
        # set up header format
        hfmt = _ENDIAN_PREFIX[self._endian_w] + _HEADER_STRUCT_FORMAT
        # write headers and data
        for header, data in traces:
            # write trace header with properly sync'd locatn field
            header['locatn'] = self._f_data.tell()
            self._f_headers.write(struct.pack(hfmt, *_HEADER_STRUCT_HANDLER['w'](header)))
            # now the data
            dfmt = _ENDIAN_PREFIX[self._endian_w] + '%if' % (header['ndata'])
            self._f_data.write(struct.pack(dfmt, *data.astype(np.float32).tolist()))
        # done - report what we wrote
        if verbose:
            print( 'wrote: (%s, %i bytes) (%s, %i bytes)' % (
                self._headers_file, self._f_headers.tell(), self._data_file, self._f_data.tell()) )


####
# TEST TEST TEST TEST
if __name__ == '__main__':
    ep = EOData('C122203C.eoh', 'C122203C.eod')
    ep2 = EOData('C122203C.eoh2', 'C122203C.eod2')
    ep.open(mode='r')
    traces = ep.read_headers()
    ep2.open(mode='w')
    data_traces = []
    for trace in traces:
        data_traces.append((trace,ep.read_data(trace)))
    ep2.write_headers_and_data(data_traces)
    ep.close()
    ep2.close()


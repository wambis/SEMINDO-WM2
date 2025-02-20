#!/usr/bin/env python

####
# Imports
import struct
import numpy as np

# Format prefix for endian-ness
_ENDIAN_PREFIX = { 'l': '<', 'b': '>' }

encoding = 'ISO-8859-1'#utf-8
#encoding = 'utf-8'

####
# First, we define a series of handlers for converting struct.unpack()
# output tuples to dicts and reverse

#
# handlers for title_st
def _handler_title_st_r(struct_tuple):
    moment = { 'rr': struct_tuple[11]
             , 'tt': struct_tuple[12]
             , 'pp': struct_tuple[13]
             , 'rt': struct_tuple[14]
             , 'rp': struct_tuple[15]
             , 'tp': struct_tuple[16]}
    return { 'event':  "".join([elem.decode(encoding) for elem in struct_tuple[:8]])
           , 'theta':  struct_tuple[8]
           , 'phi':    struct_tuple[9]
           , 'depth':  struct_tuple[10]
           , 'moment': moment
           , 'dt':     struct_tuple[17]}
def _handler_title_st_w(struct_dict):
    event  = [c.encode(encoding) for c in struct_dict['event']]
    moment = [struct_dict['moment'][k]
        for k in ['rr','tt','pp','rt','rp','tp']]
    return tuple( event
                + [struct_dict[k] for k in ['theta', 'phi', 'depth']]
                + moment
                + [struct_dict['dt']])

#
# handlers for tracehdr_st
def _handler_tracehdr_st_r(struct_tuple):
    return { 'stn':      "".join([elem.decode(encoding) for elem in struct_tuple[:4]])
           , 'locatn':   struct_tuple[4]
           , 'theta':    struct_tuple[5]
           , 'phi':      struct_tuple[6]
           , 'delta':    struct_tuple[7]
           , 'az':       struct_tuple[8]
           , 'dip':      struct_tuple[9]
           , 'smplintv': struct_tuple[10]
           , 'w1':       struct_tuple[11]
           , 'w2':       struct_tuple[12]
           , 'w3':       struct_tuple[13]
           , 'w4':       struct_tuple[14]
           , 'id':       struct_tuple[15]
           , 'reftime':  struct_tuple[16].decode(encoding)
           , 'comp':     struct_tuple[17].decode(encoding)
           , 'netwk':    "".join([elem.decode(encoding) for elem in struct_tuple[18:22]])
           , 'chnnl':    "".join([elem.decode(encoding) for elem in struct_tuple[22:26]])
           , 'extr':     "".join([elem.decode(encoding) for elem in struct_tuple[26:30]])}
def _handler_tracehdr_st_w(struct_dict):
    stn   = [c.encode(encoding) for c in struct_dict['stn']]
    netwk = [c.encode(encoding) for c in struct_dict['netwk']]
    chnnl = [c.encode(encoding) for c in struct_dict['chnnl']]
    extr  = [c.encode(encoding) for c in struct_dict['extr']]

    if isinstance(struct_dict['reftime'], str):
        struct_dict['reftime'] = struct_dict['reftime'].encode(encoding)
    if isinstance(struct_dict['comp'], str): 
        struct_dict['comp'] = struct_dict['comp'].encode(encoding)
    
    order = ['locatn', 'theta', 'phi', 'delta', 'az', 'dip', 'smplintv',
        'w1', 'w2', 'w3', 'w4', 'id', 'reftime', 'comp']
    return tuple(stn + [struct_dict[k] for k in order] + netwk + chnnl + extr)

#
# handlers for tracehdrH_st
def _handler_tracehdrH_st_r(struct_tuple):
    return { 'stn':      ''.join(struct_tuple[:4])
           , 'locatn':   struct_tuple[4]
           , 'theta':    struct_tuple[5]
           , 'phi':      struct_tuple[6]
           , 'delta':    struct_tuple[7]
           , 'az':       struct_tuple[8]
           , 'dip':      struct_tuple[9]
           , 'smplintv': struct_tuple[10]
           , 'w1':       struct_tuple[11]
           , 'w2':       struct_tuple[12]
           , 'w3':       struct_tuple[13]
           , 'w4':       struct_tuple[14]
           , 'id':       struct_tuple[15]
           , 'reftime':  struct_tuple[16]
           , 'comp':     struct_tuple[17]
           , 'netwk':    ''.join(struct_tuple[18:22])
           , 'chnnl':    ''.join(struct_tuple[22:26])
           , 'extr':     ''.join(struct_tuple[26:30])
           , 'locatnA':  struct_tuple[30]}
def _handler_tracehdrH_st_w(struct_dict):
    stn   = [c.encode(encoding) for c in struct_dict['stn']]
    netwk = [c.encode(encoding) for c in struct_dict['netwk']]
    chnnl = [c.encode(encoding) for c in struct_dict['chnnl']]
    extr  = [c.encode(encoding) for c in struct_dict['extr']]
    order = ['locatn', 'theta', 'phi', 'delta', 'az', 'dip', 'smplintv',
        'w1', 'w2', 'w3', 'w4', 'id', 'reftime', 'comp']
    return tuple(stn + [struct_dict[k] for k in order] + netwk + chnnl + extr
        + [struct_dict['locatnA']])

#
# handlers for packhdrbr_st
def _handler_packhdrbr_st_r(struct_tuple):
    return { 'phase':  struct_tuple[0]
           , 'ndata':  struct_tuple[1]
           , 'id':     struct_tuple[2]
           , 't0':     struct_tuple[3]
           , 'gv1':    struct_tuple[4]
           , 'gv2':    struct_tuple[5]
           , 'pv1':    struct_tuple[6]
           , 'pv2':    struct_tuple[7]
           , 'rmsd':   struct_tuple[8]
           , 'rmsr':   struct_tuple[9]
           , 'rmss':   struct_tuple[10]
           , 'weight': struct_tuple[11]}
def _handler_packhdrbr_st_w(struct_dict):
    order = ['phase', 'ndata', 'id', 't0', 'gv1', 'gv2', 'pv1', 'pv2', 'rmsd',
        'rmsr', 'rmss', 'weight']
    return tuple([struct_dict[k] for k in order])

#
# handlers for packhdrbr_st
def _handler_packhdrmp_st_r(struct_tuple):
    return { 'phase':  struct_tuple[0]
           , 'ndata':  struct_tuple[1]
           , 'id':     struct_tuple[2]
           , 'keep':   struct_tuple[3]
           , 't0':     struct_tuple[4]
           , 'gv1':    struct_tuple[5]
           , 'gv2':    struct_tuple[6]
           , 'pv1':    struct_tuple[7]
           , 'pv2':    struct_tuple[8]
           , 'rmsd':   struct_tuple[9]
           , 'rmsr':   struct_tuple[10]
           , 'rmss':   struct_tuple[11]
           , 'lag':    struct_tuple[12]
           , 'weight': struct_tuple[13]}
def _handler_packhdrmp_st_w(struct_dict):
    order = ['phase', 'ndata', 'id', 'keep', 't0', 'gv1', 'gv2', 'pv1', 'pv2',
        'rmsd', 'rmsr', 'rmss', 'lag', 'weight']
    return tuple([struct_dict[k] for k in order])

####
# Now, we assemble the struct formats and handlers into dicts

#
# wph/H title
_TITLE_STRUCT_FORMAT  = {'title_st': '8c3f6ff'}
_TITLE_STRUCT_HANDLER = {
    'title_st': {'r': _handler_title_st_r, 'w': _handler_title_st_w}}

#
# wph/H trace header
_TRACEHDR_STRUCT_FORMAT  = {
    'tracehdr_st':  '4ci10fh10c4c',
    'tracehdrH_st': '4ci10fh10c4ci'}
_TRACEHDR_STRUCT_HANDLER = {
    'tracehdr_st':  {'r': _handler_tracehdr_st_r,  'w': _handler_tracehdr_st_w},
    'tracehdrH_st': {'r': _handler_tracehdrH_st_r, 'w': _handler_tracehdrH_st_w}}

#
# wpd packet header
_PACKHDR_STRUCT_FORMAT  = {'packhdrbr_st': '3h2x9f', 'packhdrmp_st': '3hcx10f'}
_PACKHDR_STRUCT_HANDLER = {
    'packhdrbr_st': {'r': _handler_packhdrbr_st_r,
                     'w': _handler_packhdrbr_st_w},
    'packhdrmp_st': {'r': _handler_packhdrmp_st_r,
                     'w': _handler_packhdrmp_st_w}}

####
# Finally, the actual wavepacket data class implementation

class WPData:

    def __init__(self, headers_file, data_file,
            endian_r = 'l',
            endian_w = 'l',
            title_type = 'title_st',
            tracehdr_type = 'tracehdr_st',
            packhdr_type = 'packhdrbr_st'):
        # names
        self._headers_file = headers_file
        self._data_file = data_file
        # keyword args
        self._endian_r = endian_r
        self._endian_w = endian_w
        self._title_type = title_type
        self._tracehdr_type = tracehdr_type
        self._packhdr_type = packhdr_type
        # init
        self._open = False

    def get_available_title_types(self):
        return _TITLE_STRUCT_FORMAT.keys()

    def get_available_tracehdr_types(self):
        return _TRACEHDR_STRUCT_FORMAT.keys()

    def get_available_packhdr_types(self):
        return _PACKHDR_STRUCT_FORMAT.keys()

    def open(self, mode='r', quiet=True):
        """
        open(mode = 'r', quiet = True)

        Open the header and data files associated with this object in the
        specified mode. On failure, the resulting IOError is caught and False
        is returned (True is returned on success).

        Parameters
        ----------
        mode : char, optional
            mode must be one of ether 'r' and 'w' - defaults to 'r'
        quiet : boolean, optional
            whether to suppress warning messages on failure to open the header /
            data files under the specified mode - defaults to True

        """
        # check / store the mode
        assert mode in ['r','w']
        self._mode = mode
        # attempt to open both the wp[hH] and wpd files
        try:
            self._f_headers = open(self._headers_file, self._mode + 'b')
        except IOError as e:
            if not quiet:
                print( "Note: [I/O error] cannot open %s - %s" % (self._headers_file, e.strerror) )
            return False
        try:
            self._f_data = open(self._data_file, self._mode + 'b')
        except IOError as e:
            if not quiet:
                print( "Note: [I/O error] cannot open %s - %s" % (self._data_file, e.strerror)) 
            return False
        # success
        self._open = True
        return True

    def close(self):
        """
        close()

        Close the currently-open header and data files, which must currently
        be in an open state, otherwise an AssertionError is thrown.

        Parameters
        ----------
        none

        """
        assert self._open == True
        self._f_headers.close()
        self._f_data.close()
        self._open = False

    def read_headers(self):
        """
        read_headers()

        Returns the title header and a list of all trace headers from the
        target wp[hH] file.

        Parameters
        ----------
        none

        """
        # check
        assert self._open
        assert self._mode == 'r'
        # rewind
        self._f_headers.seek(0)
        # read title struct
        fmt = _ENDIAN_PREFIX[self._endian_r] + _TITLE_STRUCT_FORMAT[self._title_type]
        title = _TITLE_STRUCT_HANDLER[self._title_type]['r'](struct.unpack(fmt, self._f_headers.read(struct.calcsize(fmt))))
        # read header structs
        fmt = _ENDIAN_PREFIX[self._endian_r] + _TRACEHDR_STRUCT_FORMAT[self._tracehdr_type]
        raw = self._f_headers.read(struct.calcsize(fmt))
        headers = []
        while raw:
            headers.append(_TRACEHDR_STRUCT_HANDLER[self._tracehdr_type]['r'](struct.unpack(fmt, raw)))
            raw = self._f_headers.read(struct.calcsize(fmt))
        return title, headers

    def read_data(self, tracehdr):
        """
        read_data(tracehdr)

        Returns all wavepackets (packet headers and waveform data) associated
        with the supplied trace header.

        Parameters
        ----------
        tracehdr : trace header structure
            dict representing a valid wp[hH] trace header

        """
        # check
        assert self._open
        assert self._mode == 'r'
        # jump to data start
        self._f_data.seek(tracehdr['locatn'])
        # read wavepacket header structs and data
        fmt = _ENDIAN_PREFIX[self._endian_r] + _PACKHDR_STRUCT_FORMAT[self._packhdr_type]
        raw = self._f_data.read(struct.calcsize(fmt))
        data = []
        while raw:
            packhdr = _PACKHDR_STRUCT_HANDLER[self._packhdr_type]['r'](struct.unpack(fmt, raw))
            if packhdr['id'] != tracehdr['id']:
                break
            data.append((packhdr, np.fromfile(self._f_data, dtype=_ENDIAN_PREFIX[self._endian_r] + 'f4', count=packhdr['ndata'])))
            raw = self._f_data.read(struct.calcsize(fmt))
        return data

    def write_headers_and_data(self, title, traces, verbose=True):
        """
        write_headers_and_data(title, traces, verbose=True)

        Builds a wp[hH] / wpd file pair from the supplied title header and list
        of traces.

          Trace list is of the form [(tracehdr1,[(packhdr1,data1),...]),...]

          Consistent tracehdr 'locatn' parameters are ensured

        Parameters
        ----------
        title : dict
            dictionary containing the wp[hH] title structure
        traces : list traces
            see above
        verbose : boolean, optional
            display number of bytes written to the wp[hH] and wpd files

        """
        # check
        assert self._open
        assert self._mode == 'w'
        # will fill in locatn entries as needed ...
        # make sure we are rewound
        self._f_headers.seek(0)
        self._f_data.seek(0)
        # write title
        fmt = _ENDIAN_PREFIX[self._endian_w] + _TITLE_STRUCT_FORMAT[self._title_type]
        self._f_headers.write(struct.pack(fmt, *_TITLE_STRUCT_HANDLER[self._title_type]['w'](title)))
        # write wp data and headers
        for tracehdr, data in traces:
            # write trace header with properly sync'd locatn field
            tracehdr['locatn'] = self._f_data.tell()
            fmt = _ENDIAN_PREFIX[self._endian_w] + _TRACEHDR_STRUCT_FORMAT[self._tracehdr_type]
            self._f_headers.write(struct.pack(fmt, *_TRACEHDR_STRUCT_HANDLER[self._tracehdr_type]['w'](tracehdr)))
            # now the wavepackets
            for packhdr, wp in data:
                # packhdr is easy ...
                fmt = _ENDIAN_PREFIX[self._endian_w] + _PACKHDR_STRUCT_FORMAT[self._packhdr_type]
                self._f_data.write(struct.pack(fmt, *_PACKHDR_STRUCT_HANDLER[self._packhdr_type]['w'](packhdr)))
                # however, numpy does not offer a "tofile" function, so we go with struct
                fmt = _ENDIAN_PREFIX[self._endian_w] + '%if' % (packhdr['ndata'])
                self._f_data.write(struct.pack(fmt, *wp.astype(np.float32).tolist()))
        # done - report what we wrote
        if verbose:
            print( 'wrote: (%s, %i bytes) (%s, %i bytes)' % (
                self._headers_file, self._f_headers.tell(), self._data_file, self._f_data.tell()) )


####
# TEST TEST TEST TEST
#if __name__ == '__main__':
#    wp = WPData('C122203C.wph', 'C122203C.wpd')
#    wp2 = WPData('C122203C.wph2', 'C122203C.wpd2')
#    wp.open(mode='r')
#    title, traces = wp.read_headers()
#    wp2.open(mode='w')
#    data_traces = []
#    for trace in traces:
#        data_traces.append((trace,wp.read_data(trace)))
#    wp2.write_headers_and_data(title, data_traces)
#    wp.close()
#    wp2.close()


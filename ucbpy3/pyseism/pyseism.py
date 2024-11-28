
import numpy as np
import pyseism_driver as _ps

_JCOM = {'T': 2, 'S': 3}

class PySeism:
    def __init__(self, source, passband, component, datatype):
        assert component in ['L','T','Z']
        assert datatype in ['a','v','d']
        self._component = component
        self._datatype = datatype
        self._source = source
        self._passband = passband
        _ps._initialize(
            self._source['lon'], self._source['lat'], self._source['dep'],
            self._source['mt_rr'], self._source['mt_tt'], self._source['mt_pp'], 
            self._source['mt_rt'], self._source['mt_rp'], self._source['mt_tp'], 
            self._passband[0], self._passband[1], self._passband[2], self._passband[3], 
            self._component, self._datatype)
    def change_component(self, component):
        assert component in ['L','T','Z']
        self._component = component
        _ps._initialize(
            self._source['lon'], self._source['lat'], self._source['dep'],
            self._source['mt_rr'], self._source['mt_tt'], self._source['mt_pp'], 
            self._source['mt_rt'], self._source['mt_rp'], self._source['mt_tp'], 
            self._passband[0], self._passband[1], self._passband[2], self._passband[3], 
            self._component, self._datatype)
    def change_source(self, source):
        self._source = source
        _ps._initialize(
            self._source['lon'], self._source['lat'], self._source['dep'],
            self._source['mt_rr'], self._source['mt_tt'], self._source['mt_pp'], 
            self._source['mt_rt'], self._source['mt_rp'], self._source['mt_tp'], 
            self._passband[0], self._passband[1], self._passband[2], self._passband[3], 
            self._component, self._datatype)
    def change_passband(self, passband):
        self._passband = passband
        _ps._initialize(
            self._source['lon'], self._source['lat'], self._source['dep'],
            self._source['mt_rr'], self._source['mt_tt'], self._source['mt_pp'], 
            self._source['mt_rt'], self._source['mt_rp'], self._source['mt_tp'], 
            self._passband[0], self._passband[1], self._passband[2], self._passband[3], 
            self._component, self._datatype)
    def calculate_seismogram(self, lon, lat, t0, dt, ndata, ellipticity=True):
        if ellipticity:
            no_ell = 0
        else:
            no_ell = 1
        u = np.zeros((ndata,),'d')
        _ps._calculate_seismogram(lon, lat, no_ell, t0, dt, u)
        return u

#!/usr/bin/env python

import pyseism
import numpy as np
import matplotlib.pyplot as plt

source = {
    'lon': 0.0,
    'lat': 0.0,
    'dep': 15.0,
    'mt_rr': 0.0,
    'mt_tt': 0.0,
    'mt_pp': 0.0,
    'mt_rt': 0.0,
    'mt_rp': 0.0,
    'mt_tp': 1e+18}

bandpass = np.array([400.0,250.0,53.0,40.0])
passband = 2.0 * np.pi / bandpass
datatype = 'a'
component = 'T'

ps = pyseism.PySeism(source, passband, component, datatype)

u = ps.calculate_seismogram(100.0, 0.0, 0.0, 1.0, 5000, ellipticity=False)
plt.plot(u)

u = ps.calculate_seismogram(60.0, 0.0, 0.0, 1.0, 5000, ellipticity=False)
plt.plot(u)

ps.change_component('Z')
u = ps.calculate_seismogram(40.0, 60.0, 0.0, 1.0, 5000, ellipticity=False)
plt.plot(u)



plt.show()

# UCB python modules

Standalone python modules for manipulating various UCB data file formats, model
files, etc.

## Contents

### Maintained

* `EODataRW.py` -- Manipualtion of `eod` trace-data files
* `Model1D.py` -- Manipulation of minos / yannos style tabular 1D models files
* `ModelA3d.py` -- Manipualtion of A3d-format spline-based seismic models
* `PhaseCodes.py` -- Mapping from UCB wavepacket phase codes to human-readable form
* `pyminos` -- Minos wrapper (eigenfunctions / dispersion only) for python; **Requires compilation**
* `pyseism` -- Computation of 1D mode synthetics using Yannos formatted eigenfunctions; **Requires compilation**
* `pyspl` -- Python wrapper for C-based A3d spline evaluation; **Requires compilation**
* `Scalings.py` -- Common scalings between seismologically relevant material quantities
* `Sphere.py` -- Common spherical geometry operations (distance, azimuth, geocentric / geodetic coversion, etc.)
* `UCBColorMaps.py` -- Colormaps (`matplotlib`) from Lekic and Romanowicz (2011)
* `UCBFilter.py` -- Frequency-domain cosine-taper bandpass filter
* `WPDataRW.py` -- Manipulation of `wpd` wavepacket-data files

### Legacy (unmaintained)

* `EOData.py` -- An older, read-only module supporting manipulation of `eod` trace files
* `Splines.py` -- A (slower) pure-python implementation of A3d spline evaluation
* `WPData.py` -- An older, read-only module supporting manipulation of `wpd` wavepacket files

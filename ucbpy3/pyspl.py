#!/usr/bin/env python

import numpy as np
import scipy.sparse as sp
import pyspl_driver as _pd


class CubicBSplines:
    """
    Cubic b-splines on arbitrarily-spaced knots in the style of
    Megnin & Romanowicz (2000).
    """
    # constructor
    def __init__(self, knots):
        self._knots = np.asarray(knots).flatten()
        self._N = self._knots.size
    # stringification
    def __str__(self):
        return 'CubicBSplines: %i spline knots' % (self._N)
    # getters / setters ...
    def get_knots(self):
        return self._knots
    def set_knots(self, knots):
        self._knots = np.asarray(knots).flatten()
        self._N = self._knots.size
    # evaluation
    def evaluate(self, radii, deriv=0, **kwargs):
        for kw in kwargs:
            if not kw in ['k']:
                raise ValueError('Unsupported keyword argument: "%s"' % (kw))
        if 'k' in kwargs:
            ks = np.asarray(kwargs['k']).flatten()
        else:
            ks = np.arange(self._N)
        radii = np.asarray(radii).flatten()
        res = []
        for k in ks:
            res_k = np.zeros(radii.shape)
            nres = _pd._bspl_driver(deriv, self._knots, radii, int(k), res_k)
            assert nres == radii.size
            res.append(sp.csr_matrix(res_k))
        return sp.vstack(res).T.tocsr()


class SphericalSplines:
    """
    Spherical splines on knots derived from an icosahedral tesselation of the
    sphere in the manner of Wang & Dahlen (1995).

    Only tesselations refined in powers of two are supported - i.e. for each
    edge of the parent tesselation, the next "level" contains two.

    Level 1 corresponds to the base icosahedron.
    """
    # mean inter-knot distances by level
    _std_sspl_dist = [0, 63.4, 31.7, 15.8, 7.9, 3.95, 1.98, 0.99]
    # constructor
    def __init__(self, knot_lons, knot_lats, knot_levels):
        self._knot_lons = np.asarray(knot_lons).flatten()
        self._knot_lats = np.asarray(knot_lats).flatten()
        self._knot_levels = np.asarray(knot_levels).flatten().astype(int)
        self._std_dist = np.asarray([self._std_sspl_dist[l] for l in self._knot_levels]).flatten()
        self._N = self._std_dist.size
    # stringification
    def __str__(self):
        return 'SphericalSplines: %i spline knots' % (self._N)
    # getters / setters ...
    def get_knots_and_levels(self):
        return self._knots, self._levels
    def set_knots_and_levels(self, knot_lons, knot_lats, knot_levels):
        self._knot_lons = np.asarray(knot_lons).flatten()
        self._knot_lats = np.asarray(knot_lats).flatten()
        self._knot_levels = np.asarray(knot_levels).flatten().astype(int)
        self._std_dist = np.asarray([self._std_sspl_dist[l] for l in self._knot_levels]).flatten()
        self._N = self._std_dist.size
    # evaluation
    def evaluate(self, lons, lats, deriv=0, **kwargs):
        lons = np.asarray(lons).flatten()
        lats = np.asarray(lats).flatten()
        assert lons.size == lats.size
        if 'nnz_hint' in kwargs:
            nnz_hint = kwargs['nnz_hint']
        else:
            nnz_hint = 20 * lons.size
        resvals = np.zeros(nnz_hint)
        resinds = np.zeros(nnz_hint).astype(np.int32)
        reslocs = np.zeros(nnz_hint).astype(np.int32)
        nres = _pd._sspl_driver(deriv, self._knot_lons, self._knot_lats, self._std_dist, lons, lats, resvals, resinds, reslocs)
        assert nres > 0
        return sp.csr_matrix(
                (resvals[:nres], (reslocs[:nres], resinds[:nres])),
                shape=(lons.size,self._N))

###def test_sspl():
###    knots = np.loadtxt('grid.6', skiprows=1)
###    grid = np.loadtxt('grid.7', skiprows=1)
###    ss = SphericalSplines(knots[:,0], knots[:,1], knots[:,2])
###    print ss
###    import time
###    t0 = time.time()
###    S = ss.evaluate(grid[:10000,0], grid[:10000,1])
###    print time.time() - t0
###    print S.shape
###
###if __name__ == '__main__':
###    test_sspl()

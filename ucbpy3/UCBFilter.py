#########
# Imports
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

class Filter:
    """Implements the UCB-style cosine taper filter"""

    def __init__(self, w):
        """
        Initialize the filter

        Arguments
        ---------
        w : array-like (float)
            list of frequencies (rad/s) [w1,w2,w3,w4] where w1<=w2<=w3<=w4 and
            w1/w4 are lower/upper cutoff frequencies and w2/w3 are lower/upper
            corner frequencies
        """
        self.set_passband(w)

    def get_passband(self):
        """
        Get the filter passband

        Returns
        -------
        The current filter passband (ndarray)
        """
        return self._w

    def set_passband(self, w):
        """
        Set the filter passband

        Arguments
        ---------
        w : array-like (float)
            list of frequencies (rad/s) [w1,w2,w3,w4] where w1<=w2<=w3<=w4 and
            w1/w4 are lower/upper cutoff frequencies and w2/w3 are lower/upper
            corner frequencies
        """
        if len(w) != 4:
            raise ValueError('Frequency list is not the right length (expecting 4)')
        self._w = np.asarray(w)

    def _window(self, x, p):
        W = np.zeros(x.shape)
        ind_low = np.nonzero(np.logical_and(x >= p[0], x <  p[1]))
        ind_mid = np.nonzero(np.logical_and(x >= p[1], x <= p[2]))
        ind_hi  = np.nonzero(np.logical_and(x >  p[2], x <= p[3]))
        if np.any(ind_low):
            W[ind_low] = 0.5 * (1.0 - np.cos(np.pi * (x[ind_low] - p[0]) / (p[1] - p[0])))
        if np.any(ind_mid):
            W[ind_mid] = 1.0
        if np.any(ind_hi):
            W[ind_hi]  = 0.5 * (1.0 + np.cos(np.pi * (x[ind_hi] - p[2]) / (p[3] - p[2])))
        return W

    def __call__(self, *args, **kwargs):
        return self.apply(*args, **kwargs)

    def apply(self, u, dt, taper_len = None, differentiate = 0, integrate = 0, dt_new = None):
        """
        Apply the cosine-taper filter to the supplied time series

        Arguments
        ---------
        u : array-like (float)
            input time series
        dt : float
            sampling interval of input time series
        taper_len : float (optional, default : lowest period in passband)
            taper applied to beginning and end of time series prior to FFT
        differentiate : int (optional, default: 0)
            for int n, differentiate the time series n times
        integrate : int (optional, default: 0)
            for int n, integrate the time series n times
        dt_new : float (optional)
            interpolate filtered time series to new sampling interval

        Returns
        -------
        filtered (possibly differentiated, integrated, or resampled) time series (ndarray)
        """
        u = np.asarray(u)
        N = u.size
        t = dt * np.arange(0, N)
        w_u = 2.0 * np.pi * np.fft.fftfreq(N, dt)
        if not taper_len:
            taper_len = 2.0 * np.pi / self._w[0]
        Ct = self._window(t, [0.0, taper_len, t[-1] - taper_len, t[-1]])
        if np.any(self._w):
            Cw = self._window(np.abs(w_u), self._w)
        else:
            Cw = np.ones(w_u.shape)
        Fu = np.fft.fft(Ct * u)
        if differentiate:
            for _ in range(differentiate):
                Fu *= 1j * w_u
        elif integrate:
            for _ in range(integrate):
                Fu[w_u != 0.0] /= 1j * w_u[w_u != 0.0]
        u_filt = np.fft.ifft(Cw * Fu).real
        if dt_new:
            resample = InterpolatedUnivariateSpline(t, u_filt)
            u_filt = resample(np.arange(0, N * dt, dt_new))
        return u_filt

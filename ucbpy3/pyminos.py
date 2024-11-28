
import numpy as np
import pyminos_driver as _pm

_JCOM = {'T': 2, 'S': 3}

class PyMinos:
    def __init__(self, wmin, wmax, lmin, lmax, nmin, nmax, eps = 1.0e-9, num_modes_hint = 10000, **kwargs):
        self._wmin = wmin
        self._wmax = wmax
        self._lmin = lmin
        self._lmax = lmax
        self._nmin = nmin
        self._nmax = nmax
        self._eps  = eps
        self._num_modes_hint = num_modes_hint
        if 'wgrv' in kwargs:
            self._wgrv = kwargs['wgrv']
        else:
            self._wgrv = self._wmax
    def run_minos(self, model, mode_type, store_eigs = False):
        # check mode_type
        if not mode_type in ['S', 'T']:
            raise ValueError('Bad mode type selection "%s"' % (mode_type))
        # preallocate storage
        n = np.zeros((self._num_modes_hint,),'i')
        l = np.zeros((self._num_modes_hint,),'i')
        w = np.zeros((self._num_modes_hint,),'d')
        c = np.zeros((self._num_modes_hint,),'d')
        U = np.zeros((self._num_modes_hint,),'d')
        q = np.zeros((self._num_modes_hint,),'d')
        if store_eigs:
            if mode_type == 'S':
                neig = 6
                eig_fields = ['U','dU','V','dV','P','dP']
            else:
                neig = 2
                eig_fields = ['W','dW']
            nvec = neig * model.num_layers
            eig = np.zeros((self._num_modes_hint * nvec,),'d')
        # run pyminos
        _pm._configure_minos(self._eps, self._wgrv, _JCOM[mode_type], self._lmin, self._lmax, 
                             self._wmin, self._wmax, self._nmin, self._nmax)
        if store_eigs:
            num_modes = _pm._run_minos_eig(model.name, model.t_ref, model.num_layers, 
                                           model.index_icb + 1, model.index_cmb + 1, model.num_ocean_layers, 
                                           model.table.reshape((model.table.size,)), 
                                           n, l, w, c, U, q, eig)
        else:
            num_modes = _pm._run_minos(model.name, model.t_ref, model.num_layers, 
                                       model.index_icb + 1, model.index_cmb + 1, model.num_ocean_layers, 
                                       model.table.reshape((model.table.size,)), 
                                       n, l, w, c, U, q)
        # sort by overtone branch / degree and return
        ind_sort = np.lexsort((l[:num_modes],n[:num_modes]))
        n = n[:num_modes][ind_sort]
        l = l[:num_modes][ind_sort]
        w = w[:num_modes][ind_sort]
        c = c[:num_modes][ind_sort]
        U = U[:num_modes][ind_sort]
        q = q[:num_modes][ind_sort]
        res = {'n': n, 'l': l, 'w': w, 'c': c, 'U': U, 'q': q}
        # parse eigenfunctions
        if store_eigs:
            eig = eig.reshape((self._num_modes_hint, nvec))
            eig = eig[:num_modes,:]
            res['eig'] = [
                    dict(zip(eig_fields,np.hsplit(eig[ind,:],neig)))
                    for ind in ind_sort]
        # TODO: add normalization consts?
        return res

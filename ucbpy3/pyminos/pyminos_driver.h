
#if ! defined(_PYMINOS_DRIVER_H)

#define _PYMINOS_DRIVER_H

void
_configure_minos( double eps, double wgrav, int jcom, int lmin, int lmax, 
                  double wmin, double wmax, int nmin, int nmax );

int
_run_minos( char *model_name, double tref, int npts, int icb, int cmb, int noc, 
            double *model_vec, int ndata,
            int *n_, int ndim_, 
            int *l_, int ldim_,
            double *w_, int wdim_,
            double *c_, int cdim_,
            double *U_, int Udim_,
            double *q_, int qdim_ );

int
_run_minos_eig( char *model_name, double tref, int npts, int icb, int cmb, int noc, 
                double *model_vec, int ndata,
                int *n_, int ndim_, 
                int *l_, int ldim_,
                double *w_, int wdim_,
                double *c_, int cdim_,
                double *U_, int Udim_,
                double *q_, int qdim_, 
                double *eig_, int eigdim_ );

#endif

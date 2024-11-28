
#if ! defined(_PYSEISM_DRIVER_H)

#define _PYSEISM_DRIVER_H

void _initialize(
	double srclon_, double srclat_, double srcdep_,
	double mt_rr_, double mt_tt_, double mt_pp_, 
	double mt_rt_, double mt_rp_, double mt_tp_, 
	double w1_, double w2_, double w3_, double w4_,
	char compnt_, char datatype_);

int _calculate_seismogram(double stnlon, double stnlat,
	int no_ell, double t0, double dt, double *u_t_, int ndata);

#endif

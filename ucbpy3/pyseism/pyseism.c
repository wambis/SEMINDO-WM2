#include <stdio.h>
#include <stdlib.h>

#define CONV 0.0174533

#define NDATA_MAX 100000

#include "structsynth.h"

void seism2y_ellrot_auto(int, source_st *, recvlctn_st *, syntrace_st *, int *, int, int, float *);

static int configured = 0, reinit = 1;
static double srclon, srclat, srcdep,
              mt_rr, mt_tt, mt_pp, 
              mt_rt, mt_rp, mt_tp, 
              w1, w2, w3, w4;
static char compnt, datatype;

void _initialize(
	double srclon_, double srclat_, double srcdep_,
	double mt_rr_, double mt_tt_, double mt_pp_, 
	double mt_rt_, double mt_rp_, double mt_tp_, 
	double w1_, double w2_, double w3_, double w4_,
	char compnt_, char datatype_)
{
	srclon   = srclon_;
	srclat   = srclat_;
	srcdep   = srcdep_;
	mt_rr    = mt_rr_;
	mt_tt    = mt_tt_;
	mt_pp    = mt_pp_;
	mt_rt    = mt_rt_;
	mt_rp    = mt_rp_;
	mt_tp    = mt_tp_;
	w1       = w1_;
	w2       = w2_;
	w3       = w3_;
	w4       = w4_;
	compnt   = compnt_;
	datatype = datatype_;
	configured = 1;
	reinit = 1;
}

int _calculate_seismogram(double stnlon, double stnlat,
	int no_ell, double t0, double dt, double *u_t_, int ndata)
{
	int i, nmode;
	float u_t[NDATA_MAX];
	source_st src;
	recvlctn_st rcv;
	syntrace_st trc;

	if (ndata > NDATA_MAX) {
		fprintf(stderr, "Error [%s]: ndata > NDATA_MAX\n", __func__);
		return 0;
	}

	src.moment[0] = mt_rr;
	src.moment[1] = mt_tt;
	src.moment[2] = mt_pp;
	src.moment[3] = mt_rt;
	src.moment[4] = mt_rp;
	src.moment[5] = mt_tp;

	src.theta = CONV * (90.0 - srclat);
	src.phi   = CONV * srclon;
	src.dep   = 1000 * srcdep;

	rcv.theta = CONV * (90.0 - stnlat);
	rcv.phi   = CONV * stnlon;

	trc.t0 = t0;
	trc.whghcut = w1;
	trc.whghcorner = w2;
	trc.wlowcorner = w3;
	trc.wlowcut = w4;
	trc.whichbrnch = 'a';
	trc.com = compnt;
	trc.pv1 = trc.gv1 = 0.;
	trc.pv2 = trc.gv2 = 1.;
	if (compnt == 'Z' || compnt == 'L')
		trc.modetype = 'S';
	else
		trc.modetype = 'T';
	trc.datatype = datatype;
	trc.dt = dt;
	trc.ndata = ndata;

	seism2y_ellrot_auto(reinit, &src, &rcv, &trc, &nmode, 0, no_ell, u_t);
	reinit = 0;

	for (i = 0; i < ndata; i++)
		u_t_[i] = (double) u_t[i];

	return ndata;
}
